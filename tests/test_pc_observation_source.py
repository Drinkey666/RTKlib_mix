"""Mutually exclusive PC inputs; isolated outputs, no reference files overwritten."""
import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
PROTECTED = ('ppp.c', 'rtklib.h', 'rtkpos.c', 'smartphone_ppp_config.c',
             'gnss_adapter.c', 'gnss_signal_policy.c')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--replay-exe', type=Path, required=True)
    args = parser.parse_args()
    out = ROOT / 'gnss_replay/observation_source_audit' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out.mkdir(parents=True)
    hashes = {name: hashlib.sha256((ROOT / 'src' / name).read_bytes()).hexdigest() for name in PROTECTED}
    # These are the already-existing standalone replay example products, not PPP overrides.
    text = (ROOT / 'src/gnss_replay.c').read_text(encoding='utf-8')
    products = dict(re.findall(r'a->(nav|sp3|clk|bia|ionex|vmf0|vmf6|orog|atx) = "([^"]+)";', text))
    products = {k: v.replace('\\\\', '\\') for k, v in products.items()}
    txt = Path('E:/RTKLIB_Data/OBS/RawData_20260922_141724.txt')
    assert txt.is_file()
    for value in products.values():
        assert Path(value).is_file(), value
    mapping = {'nav':'NAV', 'sp3':'SP3', 'clk':'CLK', 'bia':'BIA', 'ionex':'IONEX',
               'vmf0':'VMF_FIRST', 'vmf6':'VMF_SECOND', 'orog':'OROGRAPHY', 'atx':'ATX'}
    env = {k: v for k, v in os.environ.items() if not k.startswith('RTK_PPP_')}
    env.update({'RTK_PPP_' + mapping[k]: v for k, v in products.items()})
    env.update(RTK_PPP_TXT=str(txt), RTK_PPP_TRACE_LEVEL='2')
    results = []

    def run(name, flags, changes=None, standalone=False):
        folder = out / name
        folder.mkdir()
        current_env = dict(env)
        current_env.update(RTK_PPP_OUTPUT=str(folder / 'solution.csv'), RTK_PPP_TRACE=str(folder / 'solution.trace'))
        for key, value in (changes or {}).items():
            if value is None:
                current_env.pop(key, None)
            else:
                current_env[key] = value
        if standalone:
            cmd = [str(args.replay_exe), '--source', 'txt-fast', '--txt', str(txt),
                   '--signal-mode', 'ppp-safe', '--adr-unc-max', '1.0',
                   '--trace-level', '2', '--trace', str(folder/'solution.trace'),
                   '--result-file', str(folder/'solution.csv'), '--state-dump', str(folder/'solution.csv.stat'),
                   '--obs-dump', str(folder/'solution.csv.obs.csv'), *flags]
            for key, value in products.items():
                cmd.extend(['--'+key, value])
        else:
            cmd = [str(args.exe), *flags]
        started = time.perf_counter()
        p = subprocess.run(cmd, env=current_env, input=b'\n', stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=180)
        elapsed = time.perf_counter() - started
        (folder/'console.log').write_bytes(p.stdout)
        log = p.stdout.decode('utf-8', errors='replace')
        item = {'test':name, 'exit':p.returncode, 'elapsed_s':round(elapsed, 3)}
        perf = re.search(r'PERF,([^\r\n]+)', log)
        if perf:
            item['performance'] = dict(x.split('=', 1) for x in perf.group(1).split(','))
        path = folder/'solution.csv'
        if path.is_file():
            rows = list(csv.DictReader(path.open(newline='')))
            item['epochs'] = len(rows)
            item['Q5'] = sum(row['Q']=='5' for row in rows)
            item['Q6'] = sum(row['Q']=='6' for row in rows)
        results.append(item)
        return p, folder, log

    protected_dump = out/'protected.csv.obs.csv'
    protected_dump.write_text('path-validation fixture', encoding='utf-8')
    failures = [
        ('both_inputs', ['--check-files'], {'RTK_PPP_OBS':str(txt)}),
        ('rinex_with_txt', ['--source', 'rinex', '--check-files'], {}),
        ('txt_with_obs', ['--source', 'txt-realtime', '--check-files'], {'RTK_PPP_TXT':None, 'RTK_PPP_OBS':str(txt)}),
        ('missing_txt', ['--source', 'txt-realtime', '--check-files'], {'RTK_PPP_TXT':None}),
        ('missing_output', ['--source', 'txt-fast', '--check-files'], {'RTK_PPP_OUTPUT':None}),
        ('invalid_epoch_limit', ['--source', 'txt-fast', '--max-epochs', '-1'], {}),
        ('protect_input', ['--check-files'], {'RTK_PPP_OUTPUT':str(txt)}),
        ('protect_dump', ['--check-files'], {'RTK_PPP_OUTPUT':str(out/'protected.csv'), 'RTK_PPP_NAV':str(protected_dump)}),
    ]
    for name, flags, changes in failures:
        p, folder, _ = run(name, flags, changes)
        assert p.returncode == 2, (name, p.stdout[-1000:])
        assert not list(folder.glob('*.trace')) and not list(folder.glob('*.stat'))
    p, _, log = run('txt_file_check', ['--check-files'])
    assert p.returncode == 0 and 'SOURCE=ANDROID_TXT' in log

    p, fast, log = run('txt_fast20', ['--source', 'txt-fast', '--max-epochs', '20'])
    assert p.returncode == 0, p.stdout[-2000:]
    assert 'rinex_obs=NONE' in log and 'COMPARE,' not in log
    p, realtime, log = run('txt_realtime20', ['--source', 'txt-realtime', '--max-epochs', '20'])
    assert p.returncode == 0, p.stdout[-2000:]
    assert results[-1]['elapsed_s'] >= 19, 'Not a 1x realtime schedule'
    p, baseline, _ = run('standalone20', ['--max-epochs', '20'], standalone=True)
    assert p.returncode == 0, p.stdout[-2000:]

    def compare(a, b):
        left = list(csv.DictReader((a/'solution.csv').open(newline='')))
        right = list(csv.DictReader((b/'solution.csv').open(newline='')))
        assert len(left) == len(right) == 20
        for x, y in zip(left, right):
            assert {k:v for k,v in x.items() if k != 'processing_ms'} == {k:v for k,v in y.items() if k != 'processing_ms'}
        for suffix in ('.stat', '.obs.csv'):
            assert (a/('solution.csv'+suffix)).read_bytes() == (b/('solution.csv'+suffix)).read_bytes(), suffix
        results.append({'test':a.name+'_vs_'+b.name, 'all_coordinates_and_states_equal':True,
                        'obs_P_L_D_SNR_LLI_equal':True, 'epochs':20, 'ENU_RMS_m':[0,0,0]})
    compare(fast, realtime)
    compare(fast, baseline)

    # Verify unchanged PPP environment options really reach the embedded runner.
    p, _, log = run('txt_options_forwarded', ['--source','txt-fast','--max-epochs','2'],
                    {'RTK_PPP_SYSTEMS':'GC'})
    assert p.returncode == 0 and 'navsys=33' in log
    p, _, log = run('txt_full_dataset', ['--source', 'txt-fast'])
    assert p.returncode == 0, p.stdout[-2000:]
    assert results[-1]['epochs'] == 639
    assert hashes == {name:hashlib.sha256((ROOT/'src'/name).read_bytes()).hexdigest() for name in PROTECTED}
    summary = {'tests':results, 'protected_source_sha256':hashes,
               'inputs':{'txt':str(txt), 'products':products},
               'note':'Full dataset fast; realtime equality measured for first 20 epochs only.'}
    (out/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print('RESULT_DIRECTORY='+str(out))


if __name__ == '__main__':
    main()
