"""Controlled statistics revision tests; truth is analysis-only, never solver input."""
import concurrent.futures
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
from statistics import mean
import subprocess
from run_height_bias_audit import dataset, options

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'gnss_replay/statistical_model_audit'
EXE = ROOT / 'gnss_replay/height_bias_audit/bin/Release/height_bias_replay.exe'
CASES = {'LEGACY': 0, 'SMOOTH': 1, 'ATM_BUDGET': 2, 'COMBINED': 3, 'CORRELATED_RESEARCH': 7}

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def run(day, case):
    out = OUT / day / case
    out.mkdir(parents=True, exist_ok=True)
    effective = options.replace(' -WGTELCN=0', '') + ' -STATMOD=' + str(CASES[case])
    assert len(effective.encode()) < 256
    env = {k:v for k,v in os.environ.items() if not k.startswith('HEIGHT_AUDIT_')}
    env.update(HEIGHT_AUDIT_OPTS=effective, HEIGHT_AUDIT_DIAG=str(out/'height_diag.csv'),
               HEIGHT_AUDIT_IONDIAG=str(out/'ion_diag.csv'))
    cmd = [str(EXE), '--source', 'rinex']
    for key,path in dataset(day).items():
        cmd += ['--'+key, path]
    cmd += ['--result-file',str(out/'solution.csv'),'--state-dump',str(out/'state.stat'),
            '--trace',str(out/'ppp.trace'),'--trace-level','3']
    (out/'inputs.json').write_text(json.dumps({'products':dataset(day),'pppopt':effective,
        'command':cmd,'source_sha256':{name:hashlib.sha256((ROOT/'src'/name).read_bytes()).hexdigest()
        for name in ('ppp.c','rtklib.h','rtkpos.c','smartphone_ppp_config.c')}},indent=2),encoding='utf-8')
    with (out/'run.log').open('wb') as log:
        subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(day,case,'complete',flush=True)

def analyze():
    results=[]
    for day in ('0317','0321'):
        ref=(30.759928,103.984162,482.253) if day=='0317' else (30.759991,103.984034,482.202)
        lat,lon=map(math.radians,ref[:2]);s,c=math.sin(lat),math.cos(lat)
        rn=6378137/math.sqrt(1-6.6943799901413165e-3*s*s)
        xyz=((rn+ref[2])*c*math.cos(lon),(rn+ref[2])*c*math.sin(lon),
             (rn*(1-6.6943799901413165e-3)+ref[2])*s)
        legacy=read(OUT/day/'LEGACY/solution.csv')
        for case in CASES:
            folder=OUT/day/case
            rows=read(folder/'solution.csv');diag=read(folder/'height_diag.csv')
            trace=(folder/'ppp.trace').read_text(errors='replace')
            if CASES[case]&2:
                first=re.search(r'\$PPP_PRODUCT_USE,[^\n]+',trace).group(0)
                assert 'ionex_constraints=0,vmf3_zwd_constraints=0' in first,first
                for m in re.finditer(r'\$ATM_STAT,[^\n]+',trace):
                    fields=dict(re.findall(r'(\w+)=([^,\s]+)',m.group(0)))
                    assert float(fields['dt'])>0
                    assert 0<float(fields['information'])<=1
                    assert math.isfinite(float(fields['effective_variance']))
            times=sorted(float(r['processing_ms']) for r in rows)
            result={'day':day,'case':case,'epochs':len(rows),'Q':dict(Counter(r['Q'] for r in rows)),
                'mean_ms':mean(times),'p95_ms':times[math.ceil(.95*len(times))-1],
                'max_ms':max(times),'deadline_miss':sum(v>1000 for v in times),
                'PPP_REJECT':trace.count('$PPP_REJECT,'),'final_ZWD':float(diag[-1]['ZWDest_m']),
                'final_ZTD_sigma':float(diag[-1]['sigma_ZTD_m'])}
            for label,segment in [('first30',rows[:30]),('tail300',rows[-300:]),('all',rows)]:
                errors=[[float(r[k])-v for k,v in zip(('x_m','y_m','z_m'),xyz)] for r in segment if r['Q']=='6']
                result[label+'_3D_RMS']=math.sqrt(mean(sum(v*v for v in e) for e in errors))
                result[label+'_U_RMS']=math.sqrt(mean((c*math.cos(lon)*e[0]+c*math.sin(lon)*e[1]+s*e[2])**2 for e in errors))
            result['delta_legacy_3D_RMS']=math.sqrt(mean(sum((float(a[k])-float(b[k]))**2
                for k in ('x_m','y_m','z_m')) for a,b in zip(rows,legacy)))
            assert len(rows)==(1200 if day=='0317' else 600)
            assert not result['Q'].get('0',0),result
            assert all(math.isfinite(float(r[k])) for r in rows for k in ('x_m','y_m','z_m'))
            if case=='LEGACY':
                old=ROOT/'gnss_replay/signal_quality_audit'/day/'CLEAN'
                assert (folder/'state.stat').read_bytes()==(old/'state.stat').read_bytes(), 'Legacy state differs'
                oldrows=read(old/'solution.csv')
                assert all(a[k]==b[k] for a,b in zip(rows,oldrows)
                    for k in a if k!='processing_ms'), 'Legacy observations/solutions differ'
            results.append(result)
            print(json.dumps(result),flush=True)
    (OUT/'summary.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    return results

if __name__=='__main__':
    import sys
    if '--analyze-only' not in sys.argv:
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(lambda pair:run(*pair),[(d,c) for d in ('0317','0321') for c in CASES]))
    analyze()
