"""Real-product regressions and in-memory fault injection (no truth in solver)."""
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
from collections import Counter
from statistics import mean
from run_height_bias_audit import dataset

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'gnss_replay/startup_qc_audit'
EXE=OUT/'bin/Release/startup_qc_replay.exe'

def run(day,fixture,limit=40,off=False):
    products=dataset(day)
    name=fixture+('_OFF' if off else '')
    out=OUT/day/name;out.mkdir(parents=True,exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('STARTUP_TEST_')}
    if fixture!='CLEAN':env['STARTUP_TEST_FIXTURE']=fixture
    if off:env['STARTUP_TEST_OFF']='1'
    cmd=[str(EXE),'--source','rinex']
    for key,path in products.items():cmd+=['--'+key,path]
    cmd+=['--result-file',str(out/'solution.csv'),'--state-dump',str(out/'state.stat'),
          '--trace',str(out/'ppp.trace'),'--trace-level','3']
    if limit:cmd+=['--max-epochs',str(limit)]
    manifest={'products':products,'command':cmd,'fixture':fixture,'disabled':off,
              'disabled_features':['STARTQC','SIGQC'] if off else [],
              'source_sha256':{f:hashlib.sha256((ROOT/'src'/f).read_bytes()).hexdigest()
                 for f in ('rtkpos.c','pntpos.c','ppp.c')}}
    (out/'inputs.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    with (out/'run.log').open('wb') as log:subprocess.run(cmd,env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(day,name,'complete',flush=True)
    return out

def analyze():
    records=[]
    for day in ('0317','0321'):
        clean=list(csv.DictReader((OUT/day/'CLEAN/solution.csv').open()))
        old=list(csv.DictReader((ROOT/f'gnss_replay/height_bias_audit/{day}/BASE/solution.csv').open()))
        fields=[k for k in clean[0] if k!='processing_ms']
        assert len(clean)==len(old)
        off_rows=list(csv.DictReader((OUT/day/'CLEAN_OFF/solution.csv').open()))
        assert all(all(a[k]==b[k] for k in fields) for a,b in zip(off_rows,old)),f'{day}: disabled-control numerical regression'
        assert (OUT/day/'CLEAN_OFF/state.stat').read_bytes()==(ROOT/f'gnss_replay/height_bias_audit/{day}/BASE/state.stat').read_bytes(),f'{day}: disabled-control filter state regression'
        for path in sorted((OUT/day).iterdir()):
            if not (path/'solution.csv').exists():continue
            rows=list(csv.DictReader((path/'solution.csv').open()))
            trace=(path/'ppp.trace').read_text(errors='replace')
            ms=sorted(float(r['processing_ms']) for r in rows)
            control=off_rows if path.name.endswith('_OFF') else clean
            diff=[math.sqrt(sum((float(r[k])-float(control[i][k]))**2 for k in ('x_m','y_m','z_m')))
                  for i,r in enumerate(rows) if r['Q']!='0' and control[i]['Q']!='0']
            rec={'day':day,'fixture':path.name,'epochs':len(rows),'Q':dict(Counter(r['Q'] for r in rows)),
                 'mean_ms':mean(ms),'p95_ms':ms[math.ceil(.95*len(ms))-1],'max_ms':max(ms),
                 'deadline_processing_gt_1000ms':sum(v>1000 for v in ms),
                 'max_3D_diff_clean_m':max(diff),'final_3D_diff_clean_m':diff[-1],
                 'SPP_REJECT':trace.count('$START_SPP_REJECT,'),'SPATIAL_REJECT':trace.count('$START_CODE_REJECT,'),
                 'INPUT_REJECT':trace.count('$PRE_INPUT_REJECT,'),'SPP_FAIL':trace.count('$START_SPP_FAIL,')}
            rec['control']='CLEAN_OFF' if path.name.endswith('_OFF') else 'CLEAN'
            rec['matched_valid_epochs']=len(diff)
            rec['DOP_SLIP']=trace.count('$PRE_DOP_SLIP,')
            rec['MW_SLIP']=trace.count('$PRE_MW_SLIP,')
            rec['CLOCK_RESET']=trace.count('$PRE_CLOCK_RESET,')
            if path.name in ('CLEAN','CLEAN_OFF'):
                ref=(30.759928,103.984162,482.253) if day=='0317' else (30.759991,103.984034,482.202)
                lat,lon=map(math.radians,ref[:2]);s=math.sin(lat);c=math.cos(lat)
                rn=6378137/math.sqrt(1-6.6943799901413165e-3*s*s)
                xyz=((rn+ref[2])*c*math.cos(lon),(rn+ref[2])*c*math.sin(lon),(rn*(1-6.6943799901413165e-3)+ref[2])*s)
                def error(r):
                    dx,dy,dz=[float(r[k])-v for k,v in zip(('x_m','y_m','z_m'),xyz)]
                    return (-math.sin(lon)*dx+math.cos(lon)*dy,-s*math.cos(lon)*dx-s*math.sin(lon)*dy+c*dz,
                            c*math.cos(lon)*dx+c*math.sin(lon)*dy+s*dz)
                for label,segment in [('first30',rows[:30]),('tail300',rows[-300:]),('all',rows)]:
                    errors=[error(r) for r in segment if r['Q']=='6']
                    rec[label+'_U_RMS']=math.sqrt(mean(e[2]**2 for e in errors))
                    rec[label+'_3D_RMS']=math.sqrt(mean(sum(v*v for v in e) for e in errors))
                rec['reference_provisional']=ref
            if path.name=='CLEAN':assert rec['INPUT_REJECT']==0,'unexpected invalid-value rejection on original data'
            if path.name in ('CODE_NAN','CODE_NEGATIVE','PHASE_NAN','PHASE_HUGE','DOPPLER_NAN'):assert rec['INPUT_REJECT']>0
            if path.name in ('PRIMARY_SPIKE','PERSISTENT_PRIMARY','TWO_SPIKES'):assert rec['SPP_REJECT']>0
            if path.name=='SECONDARY_SPIKE':assert rec['SPATIAL_REJECT']>0
            if path.name=='PHASE_SPIKE':assert '$PRE_DOP_SLIP,sat=1,F1,' in trace
            if path.name=='DELAYED_START':
                assert '$START_TIMER_RESTART,reason=FIRST_USABLE_SPP' in trace
                assert all(r['Q']=='0' for r in rows[:40])
                assert any(r['Q']=='6' for r in rows[40:])
            if path.name in ('CLOCK_OFFSET','SIGNAL_OFFSET'):
                # Original March17 has a natural E12 fault. Do not disguise
                # its rejection as a false positive caused by the fixture.
                assert not re.search(r'\$START_(?:SPP|CODE)_REJECT,[^\n]*sat=G',trace),'common offset misclassified'
                assert rec['Q'].get('0',0)==0,'common clock step broke filter continuity'
            if path.name=='DOPPLER_ZERO':assert 'psmvalid_G01_F1=0' in (path/'run.log').read_text(errors='replace')
            assert all(math.isfinite(float(r[k])) for r in rows for k in ('x_m','y_m','z_m'))
            records.append(rec)
            print(day,path.name,'Q',rec['Q'],'rejects',rec['SPP_REJECT'],rec['SPATIAL_REJECT'],rec['INPUT_REJECT'],'dXYZmax',round(max(diff),4))
    (OUT/'summary.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    return records

def check_pc_entry():
    out=OUT/'pc_entry';out.mkdir(exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('RTK_PPP_')}
    env.update(RTK_PPP_OUTPUT=str(out/'k80.pos'),RTK_PPP_TRACE=str(out/'pc.trace'))
    with (out/'run.log').open('wb') as log:
        subprocess.run([str(ROOT/'x64/Release/rtklib_Project2.exe')],cwd=ROOT,env=env,
                       input=b'\n',stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
    pc=(out/'k80.pos.stat').read_text().splitlines()
    replay=(OUT/'0317/CLEAN/state.stat').read_text().splitlines()
    assert pc==replay,'PC console / replay state mismatch'
    print('PASS: PC console / replay identical state lines:',len(pc))

if __name__=='__main__':
    import sys
    if '--pc-entry' in sys.argv:
        check_pc_entry();sys.exit(0)
    if '--analyze-only' not in sys.argv:
        fixtures=['PRIMARY_SPIKE','PERSISTENT_PRIMARY','SECONDARY_SPIKE','TWO_SPIKES',
                  'CODE_NAN','CODE_NEGATIVE','PHASE_NAN','PHASE_HUGE','DOPPLER_NAN','DOPPLER_ZERO','CLOCK_OFFSET','SIGNAL_OFFSET','PHASE_SPIKE']
        jobs=[(day,'CLEAN',0,False) for day in ('0317','0321')]
        jobs += [(day,'CLEAN',0,True) for day in ('0317','0321')]
        jobs += [(day,f,40,False) for day in ('0317','0321') for f in fixtures]
        jobs += [('0317','PRIMARY_SPIKE',40,True),('0317','SECONDARY_SPIKE',40,True)]
        jobs += [(day,'DELAYED_START',80,False) for day in ('0317','0321')]
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(lambda x:run(*x),jobs))
    analyze()
    print('PASS: disabled-control filter states unchanged; original-data changes and injected faults checked.')
