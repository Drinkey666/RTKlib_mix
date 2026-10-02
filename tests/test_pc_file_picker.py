"""Isolated PC entry regression. Never overwrites original PPP results."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--baseline-exe', type=Path)
    args = parser.parse_args()
    out = ROOT / 'gnss_replay/file_picker_audit' / ('run_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True)
    results = []
    env = {k: v for k, v in os.environ.items() if not k.startswith('RTK_PPP_')}
    env.update(RTK_PPP_OUTPUT=str(out / 'unused.pos'), RTK_PPP_TRACE=str(out / 'unused.trace'))

    def run(name, flags, changes=None, exe=None, expected=0):
        e = dict(env)
        e.update(changes or {})
        p = subprocess.run([str(exe or args.exe), *flags], input=b'\n', env=e,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        (out / (name + '.log')).write_bytes(p.stdout)
        assert (p.returncode == 0) == (expected == 0), (name, p.returncode, p.stdout[-1500:])
        results.append({'test': name, 'exit': p.returncode})
        return p

    run('help', ['--help'])
    run('unknown_argument', ['--unknown'], expected=2)
    run('check_files', ['--check-files'])
    run('missing_obs', ['--no-gui'], {'RTK_PPP_OBS': str(out/'missing.obs')}, expected=2)
    run('missing_vmf', ['--check-files'], {'RTK_PPP_VMF_FIRST': str(out/'missing.H00')}, expected=2)
    run('directory_as_output', ['--check-files'], {'RTK_PPP_OUTPUT': str(out)}, expected=2)
    run('missing_output_parent', ['--check-files'], {'RTK_PPP_OUTPUT': str(out/'missing_dir'/'result.pos')}, expected=2)
    inputs = json.loads((ROOT/'gnss_replay/statistical_model_audit/0317/COMBINED/inputs.json').read_text(encoding='utf-8'))['products']
    obs = inputs['rinex-obs']
    before = sha(obs)
    run('protect_obs', ['--check-files'], {'RTK_PPP_OUTPUT': obs}, expected=2)
    assert sha(obs) == before
    assert not (out/'unused.pos').exists() and not (out/'unused.trace').exists()

    for day in ('0317', '0321'):
        folder = out / ('with spaces ' + day)
        folder.mkdir()
        products = json.loads((ROOT/f'gnss_replay/statistical_model_audit/{day}/COMBINED/inputs.json').read_text(encoding='utf-8'))['products']
        mapping = {'rinex-obs':'OBS', 'nav':'NAV', 'sp3':'SP3', 'clk':'CLK', 'bia':'BIA',
                   'ionex':'IONEX', 'vmf0':'VMF_FIRST', 'vmf6':'VMF_SECOND', 'orog':'OROGRAPHY', 'atx':'ATX'}
        changes = {'RTK_PPP_'+mapping[k]:v for k,v in products.items()}
        changes.update(RTK_PPP_OUTPUT=str(folder/'selected.pos'), RTK_PPP_TRACE=str(folder/'selected.trace'))
        run(day, ['--no-gui'], changes)
        stat = folder/'selected.pos.stat'
        reference = ROOT/f'gnss_replay/statistical_model_audit/{day}/COMBINED/state.stat'
        # postpos text output is CRLF; replay's binary state dump uses LF.
        equal = stat.read_bytes().splitlines() == reference.read_bytes().splitlines()
        assert equal, f'{day}: complete PPP states differ from verified replay baseline'
        rows = [s for s in (folder/'selected.pos').read_text().splitlines() if s.strip() and not s.startswith('%')]
        q = {}
        for row in rows:
            key = row.split()[5]
            q[key] = q.get(key, 0) + 1
        results.append({'test': day+'_full_state_equal', 'equal': equal, 'rows':len(rows), 'Q':q, 'state_sha256':sha(stat)})
        if day == '0317' and args.baseline_exe:
            changes.update(RTK_PPP_OUTPUT=str(folder/'legacy.pos'), RTK_PPP_TRACE=str(folder/'legacy.trace'))
            run('0317_legacy', [], changes, exe=args.baseline_exe)
            legacy = [s for s in (folder/'legacy.pos').read_text().splitlines() if s.strip() and not s.startswith('%')]
            assert rows == legacy, 'File-picker entry altered normal positional output'
            assert stat.read_bytes() == (folder/'legacy.pos.stat').read_bytes()
            results.append({'test':'0317_legacy_coordinates_and_states_equal', 'equal':True})
    (out/'summary.json').write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(results, indent=2, ensure_ascii=False))
    print('RESULT_DIRECTORY=' + str(out))

if __name__ == '__main__':
    main()
