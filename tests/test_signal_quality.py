"""Same products/model; isolated outputs, truth never supplied to solver."""
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
from run_height_bias_audit import dataset,options

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'gnss_replay/signal_quality_audit'
EXE=ROOT/'gnss_replay/startup_qc_audit/bin/Release/startup_qc_replay.exe'

def run(day,fixture,off=False,limit=0,phase_auto=False):
    name=fixture+('_OFF' if off else '')+('_PHASE_AUTO' if phase_auto else '')
    out=OUT/day/name;out.mkdir(parents=True,exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith(('STARTUP_TEST_','SIGNAL_TEST_'))}
    if fixture!='CLEAN':env['STARTUP_TEST_FIXTURE']=fixture
    effective=options
    if off:effective=options.replace(' -WGTELCN=0','').replace(' -PPPDIAG=1','')+' -SIGQC=0 -STATMOD=0'
    if phase_auto:effective=options.replace(' -WGTELCN=0','')+' -SIGQCPH=1'
    env['SIGNAL_TEST_OPTIONS']=effective
    cmd=[str(EXE),'--source','rinex']
    for key,path in dataset(day).items():cmd+=['--'+key,path]
    cmd+=['--result-file',str(out/'solution.csv'),'--state-dump',str(out/'state.stat'),
          '--trace',str(out/'ppp.trace'),'--trace-level','3']
    if limit:cmd+=['--max-epochs',str(limit)]
    (out/'inputs.json').write_text(json.dumps({'products':dataset(day),'command':cmd,'pppopt':effective,
      'fixture':fixture,'source_sha256':{f:hashlib.sha256((ROOT/'src'/f).read_bytes()).hexdigest()
      for f in ('ppp.c','rtklib.h','rtkpos.c','pntpos.c')}},indent=2),encoding='utf-8')
    with (out/'run.log').open('wb') as log:
        subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(day,name,'complete',flush=True)

def analyze():
    summary=[]
    for day in ('0317','0321'):
        baseline=list(csv.DictReader((OUT/'before'/day/'CLEAN/solution.csv').open()))
        clean=list(csv.DictReader((OUT/day/'CLEAN/solution.csv').open()))
        for folder in sorted((OUT/day).iterdir()):
            if not (folder/'solution.csv').exists():continue
            rows=list(csv.DictReader((folder/'solution.csv').open()))
            trace=(folder/'ppp.trace').read_text(errors='replace')
            control=baseline if folder.name.startswith('CLEAN') else clean
            if folder.name.endswith('_OFF') and not folder.name.startswith('CLEAN'):
                control=list(csv.DictReader((OUT/day/'CLEAN_OFF/solution.csv').open()))
            diffs=[math.sqrt(sum((float(a[k])-float(b[k]))**2 for k in ('x_m','y_m','z_m')))
                   for a,b in zip(rows,control) if a['Q']!='0' and b['Q']!='0']
            ms=sorted(float(r['processing_ms']) for r in rows)
            record={'day':day,'case':folder.name,'epochs':len(rows),'Q':dict(Counter(r['Q'] for r in rows)),
                    'mean_ms':mean(ms),'p95_ms':ms[math.ceil(.95*len(ms))-1],'max_ms':max(ms),
                    'processing_gt_1000ms':sum(v>1000 for v in ms),
                    'quality_changes':trace.count('$SIG_QUALITY,'),'quarantine':trace.count('$PPP_QUAR,'),
                    'PPP_REJECT':trace.count('$PPP_REJECT,'),'rms_3D_diff_control':math.sqrt(mean(v*v for v in diffs)),
                    'max_3D_diff_control':max(diffs),'final_3D_diff_control':diffs[-1]}
            if folder.name.startswith('CLEAN'):
                ref=(30.759928,103.984162,482.253) if day=='0317' else (30.759991,103.984034,482.202)
                lat,lon=map(math.radians,ref[:2]);s,c=math.sin(lat),math.cos(lat)
                rn=6378137/math.sqrt(1-6.6943799901413165e-3*s*s)
                xyz=((rn+ref[2])*c*math.cos(lon),(rn+ref[2])*c*math.sin(lon),(rn*(1-6.6943799901413165e-3)+ref[2])*s)
                for label,segment in [('first30',rows[:30]),('tail300',rows[-300:]),('all',rows)]:
                    errors=[[float(r[k])-v for k,v in zip(('x_m','y_m','z_m'),xyz)] for r in segment if r['Q']=='6']
                    record[label+'_3D_RMS']=math.sqrt(mean(sum(v*v for v in e) for e in errors))
                    record[label+'_U_RMS']=math.sqrt(mean((c*math.cos(lon)*e[0]+c*math.sin(lon)*e[1]+s*e[2])**2 for e in errors))
            if folder.name=='CLEAN_OFF':
                assert record['max_3D_diff_control']==0,(day,'disabled numerical regression')
                assert (folder/'state.stat').read_bytes()==(OUT/'before'/day/'CLEAN/state.stat').read_bytes(),(day,'disabled state regression')
            if folder.name=='LATE_CODE_PRIMARY':
                assert re.search(r'\$SIG_QUALITY,[^\n]*sat=G01,sig=1C[^\n]*type=CODE[^\n]*var_factor=4',trace)
                assert re.search(r'\$SIG_QUALITY,[^\n]*sat=G01,sig=1C[^\n]*type=CODE[^\n]*var_factor=1',trace)
            if folder.name=='LATE_CODE_SECONDARY' and day=='0317':
                assert re.search(r'\$SIG_QUALITY,[^\n]*sat=G01,sig=5Q[^\n]*type=CODE[^\n]*var_factor=4',trace)
            if folder.name=='LATE_PHASE':
                assert re.search(r'\$PPP_QUAR,[^\n]*sat=G01,sig=1C',trace)
            assert all(math.isfinite(float(r[k])) for r in rows for k in ('x_m','y_m','z_m'))
            assert not record['Q'].get('0',0),(day,folder.name,'unexpected Q0')
            summary.append(record)
            print(json.dumps(record),flush=True)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary

def check_entries():
    out=OUT/'real_replay';out.mkdir(exist_ok=True)
    cmd=[str(ROOT/'gnss_replay/x64/Release/gnss_replay.exe'),'--source','rinex']
    for key,path in dataset('0317').items():cmd+=['--'+key,path]
    cmd+=['--result-file',str(out/'solution.csv'),'--state-dump',str(out/'state.stat')]
    with (out/'run.log').open('wb') as log:
        subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    assert (out/'state.stat').read_text().splitlines()==(OUT/'0317/CLEAN/state.stat').read_text().splitlines()
    print('PASS: production replay / test wrapper full filter-state equality')
    from test_startup_qc import check_pc_entry
    check_pc_entry()

if __name__=='__main__':
    import sys
    if '--check-entries' in sys.argv:
        check_entries();sys.exit(0)
    if '--analyze-only' not in sys.argv:
        jobs=[(d,'CLEAN',off,0) for d in ('0317','0321') for off in (False,True)]
        jobs += [(d,'CLEAN',False,0,True) for d in ('0317','0321')]
        jobs += [(d,f,off,180) for d in ('0317','0321') for f in
                 ('LATE_CODE_PRIMARY','LATE_CODE_SECONDARY','LATE_PHASE') for off in (False,True)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(lambda args:run(*args),jobs))
    analyze()
