"""Re-run original fault regressions into NEW folders, preserving prior results."""
import concurrent.futures
import sys
from pathlib import Path
import test_startup_qc as startup
import test_signal_quality as signal

ROOT=Path(__file__).resolve().parents[1]
startup.OUT=ROOT/'gnss_replay/statistical_model_audit/startup_regression'
signal.OUT=ROOT/'gnss_replay/statistical_model_audit/signal_regression'

if __name__=='__main__':
    if '--check-entries' in sys.argv:
        signal.check_entries()
        sys.exit(0)
    fixtures=['PRIMARY_SPIKE','PERSISTENT_PRIMARY','SECONDARY_SPIKE','TWO_SPIKES',
        'CODE_NAN','CODE_NEGATIVE','PHASE_NAN','PHASE_HUGE','DOPPLER_NAN','DOPPLER_ZERO',
        'CLOCK_OFFSET','SIGNAL_OFFSET','PHASE_SPIKE']
    jobs=[(d,'CLEAN',0,off) for d in ('0317','0321') for off in (False,True)]
    jobs += [(d,f,40,False) for d in ('0317','0321') for f in fixtures]
    jobs += [('0317','PRIMARY_SPIKE',40,True),('0317','SECONDARY_SPIKE',40,True)]
    jobs += [(d,'DELAYED_START',80,False) for d in ('0317','0321')]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda args:startup.run(*args),jobs))
    startup.analyze()
    jobs=[(d,'CLEAN',off,0) for d in ('0317','0321') for off in (False,True)]
    jobs += [(d,'CLEAN',False,0,True) for d in ('0317','0321')]
    jobs += [(d,f,off,180) for d in ('0317','0321') for f in
        ('LATE_CODE_PRIMARY','LATE_CODE_SECONDARY','LATE_PHASE') for off in (False,True)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda args:signal.run(*args),jobs))
    signal.analyze()
    print('PASS: 34 startup and 18 sustained-quality real-product regressions')
