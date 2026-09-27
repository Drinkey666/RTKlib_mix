# Android ADR → RTKLIB carrier-phase diagnosis

## Scope and controls

Sample: `E:\RTKLIB_Data\OBS\RawData_20260922_141724.txt` and matching
`E:\RTKLIB_Data\OBS\GEOL00CHN_R_20262650617_00U_01S_MO.rnx`.

All six runs use the same Release `gnss_replay.exe`, 639 matching epochs,
the same NAV/SP3/CLK/BIA/IONEX/VMF3/ATX products, `configure_ppp()` and
one persistent `rtk_t`. Only the observation source or `--adr-unc-max` varies.
No change was made to `ppp.c`, `rtklib.h`, PPP parameters, products, or signal
priority. The RINEX-source control still gives Q5=11, Q6=628, and 17,650
carrier-phase records. Its trace contains 205 `$PPP_REJECT` lines.

`--adr-unc-max off` is the default diagnostic setting. It emits L whenever
the *selected* raw measurement has `ADR_VALID`, finite nonzero ADR, even if
`ADR_RESET` or `ADR_SLIP` is set; either flag instead sets `LLI_SLIP`.
Unresolved half-cycle sets `LLI_HALFC`. Android ADR uncertainty does not enter
PPP weighting. Numeric settings apply only a phase-acceptance threshold.
The earlier extra LLI inference from missing previous phase, signal changes,
or half-cycle transitions was removed after it produced 453 unnecessary slip
flags against matching RINEX carrier phases. Explicit Android reset/slip and
clock discontinuity still flag slips.

## Five sequential runs

The table is from sequential runs, not concurrent timing. `PPP_REJECT` counts
trace lines (possibly multiple iterations), not distinct satellites.

| ADR uncertainty limit (m) | Code | Phase | Q5 | Q6 | PPP_REJECT | LLI_SLIP on emitted phase | Mean ms | P95 ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 33,278 | 542 | 632 | 7 | 162 | 4 | 3.537 | 4.431 |
| 1.0 | 33,278 | 21,599 | 11 | 628 | 235 | 29 | 4.972 | 6.959 |
| 2.0 | 33,278 | 21,599 | 11 | 628 | 235 | 29 | 4.972 | 7.025 |
| 5.0 | 33,278 | 21,599 | 11 | 628 | 235 | 29 | 4.851 | 7.047 |
| off | 33,278 | 21,599 | 11 | 628 | 235 | 29 | 4.913 | 7.041 |

The old adapter emitted 538 phases at the hard 0.5 m limit. After the
requested reset/slip change, the 0.5 m run emits 542: four reset-flagged
measurements now retain their phase. The threshold itself rejects **21,057**
otherwise eligible *selected* phases (21,599 minus 542). Removing it raises
Q6 from **7 to 628** of 639 epochs. All observed valid ADR uncertainties are
below 1.0 m; therefore 1.0, 2.0, 5.0 and off have identical counts here.

## Android ADR uncertainty distribution

The TXT contains 53,627 `Raw` records, including 11 before a usable complete
clock was available; 53,616 reached the adapter. Across all raw records:

| Flag/count | Value |
| --- | ---: |
| ADR_VALID | 35,315 |
| ADR_RESET | 49 |
| ADR_SLIP | 0 |
| ADR_HALF_REPORTED | 53,627 |
| Invalid uncertainty among ADR_VALID | 0 |

The following uncertainty distribution uses the 35,315 `ADR_VALID` records,
including constellations/signals not selected by this version's adapter:

| Statistic | ADR uncertainty (m) |
| --- | ---: |
| Minimum | 0.000225 |
| Median | 0.500625 |
| P90 / P95 | 0.502025 / 0.502025 |
| P99 | 0.504900 |
| Maximum | 0.510000 |

| Cumulative bucket | Count |
| --- | ---: |
| ≤0.5 m | 6,622 |
| ≤1.0 m | 35,315 |
| ≤2.0 m | 35,315 |
| ≤5.0 m | 35,315 |
| >5.0 m | 0 |

Thus 28,693 ADR_VALID raw records lie in (0.5, 1.0] m. Common values such
as 0.500625 and 0.500900 m explain why a strict 0.5 m comparison is so
destructive on this device. The raw-record count is not the phase-loss count:
unsupported/competing signals are removed before phase output.

## Strict TXT/RINEX observation comparison at `off`

Matching requires the same epoch, satellite, and signal code, independent of
frequency-slot position. The sign of each difference is TXT minus RINEX;
maximum is absolute. Phase uses `(L_txt - L_rinex) × wavelength` in metres.
Doppler compares `-D_txt × wavelength` with `-D_rinex × wavelength` in m/s.

| Quantity | Matched pairs | Signed mean | RMS | Maximum absolute |
| --- | ---: | ---: | ---: | ---: |
| Carrier phase (m) | 17,644 | ≈0.000000 m | 0.000066 m | 0.000127 m |
| Pseudorange (m) | 29,110 | ≈0.000000 m | 0.000289 m | 0.000500 m |
| Pseudorange rate from Doppler (m/s) | 29,110 | 0.000001 m/s | 0.000068 m/s | 0.000186 m/s |

All 17,644 matched phases also have identical `LLI` values. The maximum
differences are consistent with RINEX output precision. These checks are more
specific evidence for conversion correctness than Q6 count alone.

The paths are **not yet globally identical**: the adapter emits 21,599 phases,
of which 3,955 have no RINEX phase with the same epoch/satellite/code (3,677
`5Q`, 277 `1P`, one `1C`). The slot-based comparison also reports 4,821 code
mismatches, including 4,544 where the RINEX slot is blank and 277 with another
code. These require a separate signal-mapping/selection audit before claiming
general equivalence or positioning accuracy.

## Why six RINEX phases remain absent

Reference-centric accounting closes exactly: **17,650 = 17,644 matched + 6
missing**. Reasons among the six:

| Cause | Count |
| --- | ---: |
| NO_RAW_SIGNAL | 0 |
| UNSUPPORTED_SIGNAL | 0 |
| ADR_NOT_VALID | 0 |
| ADR_ZERO | 0 |
| ADR_NAN | 0 |
| ADR_UNCERTAINTY | 0 |
| SIGNAL_SELECTION_LOST | **6** |
| CODE_MAPPING_DIFFERENT | 0 |
| OTHER | 0 |

Five are BeiDou C23 `2I` and one is C33 `2I`. The TXT output instead selected
`1P` in their shared frequency slot. Signal priority was intentionally left
unchanged; this is the next-stage mapping/selection issue, not an ADR
uncertainty failure. The 0.5 m experiment had 17,379 missing RINEX phases
classified as `ADR_UNCERTAINTY` plus these six selection losses.

## Interpretation and reproduction

This test establishes that the 0.5 m ADR uncertainty hard limit was the
dominant cause of missing phase and Q5 fallback for this capture. It does
**not** establish centimetre positioning accuracy. More input phases also
raise `$PPP_REJECT` trace occurrences from 162 to 235 (RINEX control: 205),
and the additional TXT-only signals still need quality/mapping assessment.
Do not tune PPP weights from this Q6 result alone.

From the Visual Studio `gnss_replay` project, set command arguments to:

```text
--source txt --txt "E:\RTKLIB_Data\OBS\RawData_20260922_141724.txt" --rinex-obs "E:\RTKLIB_Data\OBS\GEOL00CHN_R_20262650617_00U_01S_MO.rnx" --mode fast --adr-unc-max off --trace "E:\RTKLIB_Data\PPP_Result\adr_off.trace"
```

Replace `off` by `0.5`, `1.0`, `2.0`, or `5.0`, and give each run a distinct
trace path. The sample's 2026-09-22 product paths remain the unchanged
defaults; use product path options when replaying another date.
