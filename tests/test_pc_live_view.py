"""Read-only display callback and cancellation regression, optional real UI run."""
import argparse
import csv
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'gnss_replay/observation_source_audit/20261002_194141_594504'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--probe', type=Path, default=ROOT/'build/live_observer/Release/live_observer_probe.exe')
    parser.add_argument('--exe', type=Path, default=ROOT/'x64/Release/rtklib_Project2.exe')
    parser.add_argument('--ui', action='store_true')
    args = parser.parse_args()
    inputs = json.loads((BASE/'summary.json').read_text(encoding='utf-8'))['inputs']
    out = ROOT/'gnss_replay/live_view_audit'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out.mkdir(parents=True)
    env = {k:v for k,v in os.environ.items() if not k.startswith('RTK_PPP_')}
    results = []

    def run(name, limit, realtime=False, stop_after=None):
        folder = out/name
        folder.mkdir()
        cmd = [str(args.probe), '--source', 'txt-realtime' if realtime else 'txt-fast',
               '--txt', inputs['txt'], '--signal-mode','ppp-safe','--adr-unc-max','1.0',
               '--trace-level','2','--max-epochs',str(limit),
               '--result-file',str(folder/'solution.csv'), '--state-dump',str(folder/'solution.csv.stat'),
               '--obs-dump',str(folder/'solution.csv.obs.csv'), '--trace',str(folder/'solution.trace')]
        for key,value in inputs['products'].items():
            cmd.extend(['--'+key,value])
        e = dict(env)
        if stop_after is not None:
            e['GNSS_OBSERVER_STOP_AFTER'] = str(stop_after)
        started = time.perf_counter()
        p = subprocess.run(cmd,env=e,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
        elapsed = time.perf_counter()-started
        (folder/'console.log').write_bytes(p.stdout+p.stderr)
        (folder/'stderr.log').write_bytes(p.stderr)
        assert p.returncode == (3 if stop_after else 0), (p.stdout+p.stderr)[-2000:]
        rows = list(csv.DictReader((folder/'solution.csv').open(newline='')))
        callbacks = [s.split(',')[1:] for s in p.stdout.decode(errors='replace').splitlines() if s.startswith('UI_EPOCH,')]
        assert len(rows) == len(callbacks) == (stop_after or limit)
        fields = ['week','tow','Q','ns','nobs','x_m','y_m','z_m','lat_deg','lon_deg','h_m']
        assert callbacks == [[row[key] for key in fields] for row in rows]
        results.append({'test':name,'epochs':len(rows),'callbacks_equal_saved_coordinates':True,
                        'exit':p.returncode,'elapsed_s':round(elapsed,3)})
        return folder,rows

    observed,rows = run('observed_fast639',639)
    old = list(csv.DictReader((BASE/'txt_full_dataset/solution.csv').open(newline='')))
    strip = lambda r:{k:v for k,v in r.items() if k!='processing_ms'}
    assert list(map(strip,rows)) == list(map(strip,old))
    for suffix in ('.stat','.obs.csv'):
        assert (observed/('solution.csv'+suffix)).read_bytes() == (BASE/'txt_full_dataset'/('solution.csv'+suffix)).read_bytes()
    results.append({'test':'observer_vs_previous639','coordinates_states_observations_identical':True})
    realtime,rows = run('observed_realtime20',20,realtime=True)
    old = list(csv.DictReader((BASE/'txt_realtime20/solution.csv').open(newline='')))
    assert list(map(strip,rows)) == list(map(strip,old))
    assert (realtime/'solution.csv.stat').read_bytes() == (BASE/'txt_realtime20/solution.csv.stat').read_bytes()
    assert results[-1]['elapsed_s'] >= 19
    run('cancel_after3',20,realtime=True,stop_after=3)

    if args.ui:
        folder=out/'native_ui20'
        folder.mkdir()
        mapping={'nav':'NAV','sp3':'SP3','clk':'CLK','bia':'BIA','ionex':'IONEX',
                 'vmf0':'VMF_FIRST','vmf6':'VMF_SECOND','orog':'OROGRAPHY','atx':'ATX'}
        e=dict(env)
        e.update({'RTK_PPP_'+mapping[k]:v for k,v in inputs['products'].items()})
        e.update(RTK_PPP_TXT=inputs['txt'],RTK_PPP_OUTPUT=str(folder/'solution.csv'),RTK_PPP_TRACE_LEVEL='2')
        log=(folder/'console.log').open('wb')
        p=subprocess.Popen([str(args.exe),'--source','txt-realtime','--live-view','--max-epochs','20'],
                           env=e,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
        print('UI_TEST_PID='+str(p.pid),flush=True)
        print('UI_TEST_FOLDER='+str(folder),flush=True)
        # User/UI test closes only this test window; no forced process termination.
        rc=p.wait(timeout=180)
        log.close()
        rows=list(csv.DictReader((folder/'solution.csv').open(newline='')))
        assert rc==0 and len(rows)==20, (rc,len(rows))
        old=list(csv.DictReader((realtime/'solution.csv').open(newline='')))
        assert list(map(strip,rows))==list(map(strip,old))
        assert (folder/'solution.csv.stat').read_bytes()==(realtime/'solution.csv.stat').read_bytes()
        results.append({'test':'native_ui_vs_headless20','epochs':20,'coordinates_and_states_identical':True})
    (out/'summary.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(json.dumps(results,indent=2))
    print('RESULT_DIRECTORY='+str(out))


if __name__=='__main__':
    main()
