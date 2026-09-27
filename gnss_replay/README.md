# Android GNSS TXT → RTKLIB PPP PC replay

This is a separate executable. It does not call the RINEX OBS parser on the
primary path and does not modify `ppp.c` or `rtklib.h`.

Build (Visual Studio 2022, x64):

```powershell
& 'E:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe' 'E:\GNSS\rtklib_app\rtklib_Project2\gnss_replay\gnss_replay.vcxproj' /p:Configuration=Release /p:Platform=x64
```

Fast replay of the sample TXT (the 2026-09-22 product paths are defaults):

```powershell
& 'E:\GNSS\rtklib_app\rtklib_Project2\gnss_replay\x64\Release\gnss_replay.exe' `
  --source txt --txt 'E:\RTKLIB_Data\OBS\RawData_20260922_141724.txt' --mode fast `
  --obs-dump 'E:\RTKLIB_Data\PPP_Result\replay_obs.csv'
```

For a controlled input comparison, use the same executable and product/PPP
configuration, but feed matching RINEX OBS records to `rtkpos()`:

```powershell
& 'E:\GNSS\rtklib_app\rtklib_Project2\gnss_replay\x64\Release\gnss_replay.exe' `
  --source rinex --rinex-obs 'E:\RTKLIB_Data\OBS\GEOL00CHN_R_20262650617_00U_01S_MO.rnx' `
  --mode fast --trace 'E:\RTKLIB_Data\PPP_Result\replay_rinex.trace'
```

`--source txt` is the default. In TXT mode RINEX OBS is only a comparison
reference. RINEX-source mode is a diagnostic control, not part of the Android
raw-data entry path. It retains GPS/Galileo/BeiDou records, one `rtkpos()`
call per epoch, and the same persistent filter state.

Signal selection for controlled PPP equivalence tests:

| `--signal-mode` | Meaning | Needs `--rinex-obs`? |
| --- | --- | --- |
| `ppp-safe` (default) | Independent, fixed PPP signal policy in `gnss_signal_policy.c` | No |
| `full` | Keep all signals currently selected from Android TXT, for research | No |
| `rinex-compatible` | Allow only the system/slot codes actually emitted by this capture's `rinex.c` reader; no hard-coded BeiDou 2I/1P choice | Yes |
| `common-only` | Before slot selection, retain only the same epoch/satellite/code present in RINEX; diagnostic upper-bound comparison | Yes |

`rinex-compatible` and `common-only` use the RINEX file as a reference and
cannot be deployed as-is in a stand-alone Android real-time pipeline. For this
capture, the learned reader codes include BDS `2I` in slot 1 and `7I` in
slot 2. `ppp-safe` permits GPS 1C/5Q, Galileo 1C/5Q and BeiDou 2I/7I/5P.
Galileo E5b remains correctly mapped in the adapter but is excluded from this
PPP policy. BDS 1P cannot displace 2I; BDS 5Q remains available in `full`
but is excluded by the safe policy. A startup `SIGNAL_PRODUCT` summary shows
unique observed satellites and loaded OSB coverage for each signal. OSB
absence alone is not a universal exclusion: current PPP has a legacy DCB
fallback for some codes. `--exclude-raw-signal C:5Q` (or another `G/E/C:code`) is a temporary
ablation test, applied before slot selection, not a production exclusion rule.

Use `--result-file FILE` for per-epoch ECEF/geodetic coordinates, status,
satellite count, formal XYZ standard deviations and time; `--state-dump FILE`
for RTKLIB's `$POS/$CLK/$RCB/$TROP/$ION/$SAT` state stream; and
`--use-dump FILE` for each signal's input flags, elevation, OSB availability,
receiver-code-bias slot label, phase-used flag and cached postfit residual.
The slot label mirrors the current PPP model for diagnostics only; it does
not configure the filter. `--trace-level 3` additionally
emits prefit/postfit per-satellite model residuals but has significant logging
overhead, so compare processing time using the default trace level 2. The
earlier four-run evidence is in `GNSS_SIGNAL_PPP_EQUIV.md`; the independent
safe-mode validation is in `GNSS_PPP_SAFE_REPORT.md`.

`--adr-unc-max 1.0` is the default setting. Numeric settings such
as `0.5`, `1.0`, `2.0`, and `5.0` accept carrier phase only when Android's
ADR uncertainty is finite, nonnegative, and at or below the specified metres.
The uncertainty does not enter PPP weighting. ADR reset/slip records retain
their phase and set `LLI_SLIP`; half-cycle unresolved sets `LLI_HALFC`.
The output includes ADR distribution, phase-missing reasons, strict
P/phase/Doppler and LLI comparisons, and a PPP rejection count from the trace.
See `ADR_PHASE_DIAG.md` for the five-setting experiment on the sample data.
`off` remains available for research. This uncertainty is not used as a PPP
weight.

For 1× playback, replace `fast` with `realtime`. `--max-epochs 10` makes a
short smoke test. `--nav`, `--sp3`, `--clk`, `--bia`, `--ionex`, `--vmf0`,
`--vmf6`, `--orog`, and `--atx` override individual products. SP3 and CLK
must cover the observation time; the program fails early otherwise. The
default VMF3 H06/H12 files bracket this sample's 06:17 GPST observations.

The TXT reader streams `Raw,` records and groups them by `TimeNanos`. The
adapter knows nothing about CSV/TXT; it receives one `android_clock_t` and an
array of `android_meas_t`. All measurements of the same satellite merge into
one `obsd_t`. `rtkinit()` is called once, `rtkpos()` once per usable epoch,
and `rtkfree()` at shutdown. Clock discontinuities and silent `FullBiasNanos`
steps are tracked across epochs. `utcTimeMillis` is not used for GNSS time.

Supported signal slots (`NFREQ=3`):

| Constellation | Slot 1 | Slot 2 | Slot 3 |
| --- | --- | --- | --- |
| GPS | L1 C/A `1C` | — | L5 `5I/5Q/5X` |
| Galileo | E1 `1C/1B/1X` | E5b `7I/7Q/7X` | E5a `5I/5Q/5X` |
| BeiDou | B1I `2I/2Q/2X` or B1C `1D/1P/1X` | B2b `7I/7Q/7X` | B2a `5D/5P/5Q/5X` |

GLONASS and QZSS are deliberately unsupported in v1. Signals sharing one
RTKLIB slot are selected first by the safe policy's signal identity and only
then by measurement quality within that signal; only one is sent to PPP.
The broader table above is adapter recognition, not the safe-policy allowlist.
`P` is metres, `L` cycles, `D` Hz, and
`D = -pseudorange_rate / wavelength`. `SNR` uses RTKLIB's 0.001 dB-Hz unit.
ADR reset/slip, half-cycle unresolved and clock-discontinuity flags map to
RTKLIB LLI. Invalid Android TOW and ADR values are not fabricated.

CSV stdout reports week, TOW, satellite count, `rtkpos()` return, solution
status, WGS-84 latitude/longitude/ellipsoidal height, and processing time.
Processing time includes adapter conversion and `rtkpos()`, not intentional
real-time waiting or optional comparison/dump I/O. Summary stderr includes
mean, P95, max, and processing times over 1000 ms. In TXT mode the optional
RINEX OBS file is only a comparison reference and never feeds the PPP solver.
RINEX conversion may omit additional Android signals or preprocess code and
phase, so the report separates common-signal value differences from
unrepresented/different code signals. Matching common P/L/D values to
0.0005 (RINEX output rounding) validates the adapter, not PPP accuracy.
