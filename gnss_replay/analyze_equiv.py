"""Compare replay outputs without modifying the PPP filter or its products."""

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

GPS_EPOCH = datetime(1980, 1, 6)


def key(row):
    return int(row["week"]), round(float(row["tow"]), 3)


def state_key(epoch):
    # rtkoutstat prints only three decimals and may round a .9995 TOW into
    # the next whole second. The replay epochs are ~1 s apart, so .01 s
    # matching is unambiguous for this capture.
    return epoch[0], round(epoch[1], 2)


def metric(values):
    if not values:
        return {"n": 0}
    return {"n": len(values), "mean": sum(values) / len(values),
            "rms": math.sqrt(sum(v * v for v in values) / len(values)),
            "max_abs": max(abs(v) for v in values)}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_states(path):
    states = defaultdict(lambda: {"rcb": {}, "ion": {}})
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            parts = line.strip().split(",")
            if len(parts) < 5 or parts[0] not in ("$CLK", "$TROP", "$RCB", "$ION"):
                continue
            epoch = int(parts[1]), round(float(parts[2]), 2)
            state = states[epoch]
            if parts[0] == "$CLK" and len(parts) >= 10:
                state["clock"] = [float(v) for v in parts[5:10]]
            elif parts[0] == "$TROP" and len(parts) >= 6 and parts[4] == "1":
                state["ztd"] = float(parts[5])
            elif parts[0] == "$RCB" and len(parts) >= 6:
                state["rcb"][parts[4]] = float(parts[5])
            elif parts[0] == "$ION" and len(parts) >= 8:
                state["ion"][parts[4]] = float(parts[7])
    return states


def enu_delta(reference, candidate):
    dx = [float(candidate[f]) - float(reference[f]) for f in ("x_m", "y_m", "z_m")]
    lat = math.radians(float(reference["lat_deg"]))
    lon = math.radians(float(reference["lon_deg"]))
    sl, cl = math.sin(lat), math.cos(lat)
    so, co = math.sin(lon), math.cos(lon)
    return [-so * dx[0] + co * dx[1],
            -sl * co * dx[0] - sl * so * dx[1] + cl * dx[2],
            cl * co * dx[0] + cl * so * dx[1] + sl * dx[2]]


def compare(ref_rows, rows, ref_states, states):
    ref = {key(row): row for row in ref_rows}
    result = {"epochs": len(rows), "matched": 0, "q6_common": 0,
              "status_mismatch": 0, "status_mismatch_epochs": [],
              "all_enu": [[], [], []], "q6_enu": [[], [], []],
              "q6_xyz": [[], [], []], "q6_std_xyz": [[], [], []],
              "q6_ns": [], "q6_ztd": [], "q6_clock": defaultdict(list),
              "q6_rcb": defaultdict(list), "q6_ion": [], "final": None}
    for row in rows:
        epoch = key(row)
        if epoch not in ref:
            continue
        base = ref[epoch]
        result["matched"] += 1
        if row["Q"] != base["Q"]:
            result["status_mismatch"] += 1
            result["status_mismatch_epochs"].append(epoch)
        diff = enu_delta(base, row)
        for j in range(3):
            result["all_enu"][j].append(diff[j])
        if row["Q"] != "6" or base["Q"] != "6":
            continue
        result["q6_common"] += 1
        for j, field in enumerate(("x_m", "y_m", "z_m")):
            result["q6_xyz"][j].append(float(row[field]) - float(base[field]))
            sf = ("std_x_m", "std_y_m", "std_z_m")[j]
            result["q6_std_xyz"][j].append(float(row[sf]) - float(base[sf]))
            result["q6_enu"][j].append(diff[j])
        result["q6_ns"].append(int(row["ns"]) - int(base["ns"]))
        s0, s1 = ref_states.get(state_key(epoch), {}), states.get(state_key(epoch), {})
        if "ztd" in s0 and "ztd" in s1:
            result["q6_ztd"].append(s1["ztd"] - s0["ztd"])
        if "clock" in s0 and "clock" in s1:
            for j, name in ((0, "GPS_clock_ns"), (2, "GAL_ISB_ns"),
                            (3, "BDS_ISB_ns")):
                result["q6_clock"][name].append(s1["clock"][j] - s0["clock"][j])
        for name in s0.get("rcb", {}).keys() & s1.get("rcb", {}).keys():
            result["q6_rcb"][name].append(s1["rcb"][name] - s0["rcb"][name])
        for sat in s0.get("ion", {}).keys() & s1.get("ion", {}).keys():
            result["q6_ion"].append(s1["ion"][sat] - s0["ion"][sat])
        result["final"] = {"epoch": epoch, "enu_m": diff,
                           "xyz_m": [result["q6_xyz"][j][-1] for j in range(3)],
                           "height_m": float(row["h_m"]) - float(base["h_m"])}
    for name in ("all_enu", "q6_enu", "q6_xyz", "q6_std_xyz"):
        result[name] = dict(zip(("E", "N", "U") if "enu" in name else
                                ("X", "Y", "Z"), map(metric, result[name])))
    for name in ("q6_ns", "q6_ztd", "q6_ion"):
        result[name] = metric(result[name])
    for name in ("q6_clock", "q6_rcb"):
        result[name] = {label: metric(v) for label, v in result[name].items()}
    return result


def trace_rejects(path):
    pattern = re.compile(r"\$PPP_REJECT,([^,]+),iter=\d+,sat=([^,]+),F(\d+),(PHASE|CODE)")
    rejects = defaultdict(lambda: [0, 0])
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = pattern.search(line)
            if not match:
                continue
            stamp = datetime.strptime(match[1], "%Y/%m/%d %H:%M:%S.%f")
            total = (stamp - GPS_EPOCH).total_seconds()
            k = (round(total, 2), match[2], int(match[3]))
            rejects[k][0 if match[4] == "PHASE" else 1] += 1
    return rejects


def accepted_residual_groups(path):
    pattern = re.compile(
        r"\$PPP_DIAG_SIG,time=[^,]+,sys=(.),sig=([^,]+),"
        r"nP=(\d+),meanP=[^,]+,rmsP=([^,]+),"
        r"nL=(\d+),meanL=[^,]+,rmsL=([^,]+)")
    groups = defaultdict(lambda: [0, 0.0, 0, 0.0])
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = pattern.search(line)
            if not match:
                continue
            np, rp, nl, rl = (int(match[3]), float(match[4]),
                              int(match[5]), float(match[6]))
            group = groups[f"{match[1]}/{match[2]}"]
            group[0] += np
            group[1] += np * rp * rp
            group[2] += nl
            group[3] += nl * rl * rl
    return {name: {"accepted_code": g[0],
                   "accepted_code_postfit_rms_m": math.sqrt(g[1] / g[0]) if g[0] else None,
                   "accepted_phase": g[2],
                   "accepted_phase_postfit_rms_m": math.sqrt(g[3] / g[2]) if g[2] else None}
            for name, g in sorted(groups.items())}


def model_residuals(path):
    pattern = re.compile(
        r"^3 (\d{4}/\d\d/\d\d \d\d:\d\d:\d\d\.\d\d) "
        r"sat=\s*(\d+) ([LP])(\d) res=\s*([-+\d.]+)")
    records = {}
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = pattern.search(line)
            if not match:
                continue
            k = (match[1], int(match[2]), int(match[4]), match[3])
            value = float(match[5])
            if k not in records:
                records[k] = [value, value, 1]
            else:
                records[k][1] = value
                records[k][2] += 1
    return records


def trace_stamp(row):
    seconds = round(int(row["week"]) * 604800 + float(row["tow"]), 2)
    stamp = GPS_EPOCH + timedelta(seconds=seconds)
    return stamp.strftime("%Y/%m/%d %H:%M:%S.%f")[:22]


def extra_signals(full_rows, rinex_rows, rejects, model):
    ref_phase = {(key(r), r["sat"], r["signal"]) for r in rinex_rows
                 if r["has_L"] == "1"}
    groups = defaultdict(lambda: {"input": 0, "cn0": [], "elevation": [],
                                  "adr_unc": [], "osb_valid": 0,
                                  "used_phase": 0, "post_phase": [],
                                  "prefit_phase": [], "postfit_phase_trace": [],
                                  "prefit_code": [], "postfit_code_trace": [],
                                  "post_code_cache": [], "reject_phase": 0,
                                  "reject_code": 0})
    for row in full_rows:
        if row["has_L"] != "1" or (key(row), row["sat"], row["signal"]) in ref_phase:
            continue
        label = (row["sat"][0], row["signal"], int(row["slot"]))
        group = groups[label]
        group["input"] += 1
        group["cn0"].append(float(row["cn0_dbhz"]))
        group["elevation"].append(float(row["elevation_deg"]))
        if float(row["adr_unc_m"]) >= 0:
            group["adr_unc"].append(float(row["adr_unc_m"]))
        group["osb_valid"] += int(row["osb_valid"])
        if row["phase_vsat"] == "1":
            group["used_phase"] += 1
            group["post_phase"].append(float(row["postfit_phase_cache_m"]))
            mk = (trace_stamp(row), int(row["sat_no"]), int(row["slot"]), "L")
            if mk in model:
                group["prefit_phase"].append(model[mk][0])
                group["postfit_phase_trace"].append(model[mk][1])
        if row["has_P"] == "1":
            group["post_code_cache"].append(float(row["postfit_code_cache_m"]))
            mk = (trace_stamp(row), int(row["sat_no"]), int(row["slot"]), "P")
            if mk in model:
                group["prefit_code"].append(model[mk][0])
                group["postfit_code_trace"].append(model[mk][1])
        rk = (round(int(row["week"]) * 604800 + float(row["tow"]), 2),
              row["sat"], int(row["slot"]))
        group["reject_phase"] += rejects[rk][0]
        group["reject_code"] += rejects[rk][1]
    return {"/".join(map(str, label)): {
        "input": group["input"], "mean_cn0_dbhz": metric(group["cn0"])["mean"],
        "mean_elevation_deg": metric(group["elevation"])["mean"],
        "adr_unc_m": metric(group["adr_unc"]),
        "osb_valid": group["osb_valid"],
        "osb_missing": group["input"] - group["osb_valid"],
        "phase_vsat": group["used_phase"],
        "phase_postfit_cache_rms_m": metric(group["post_phase"]).get("rms"),
        "phase_prefit_trace_rms_m": metric(group["prefit_phase"]).get("rms"),
        "phase_prefit_trace_n": len(group["prefit_phase"]),
        "phase_postfit_trace_rms_m": metric(group["postfit_phase_trace"]).get("rms"),
        "code_cache_rms_m": metric(group["post_code_cache"]).get("rms"),
        "code_prefit_trace_rms_m": metric(group["prefit_code"]).get("rms"),
        "code_prefit_trace_n": len(group["prefit_code"]),
        "code_postfit_trace_rms_m": metric(group["postfit_code_trace"]).get("rms"),
        "reject_phase_lines": group["reject_phase"],
        "reject_code_lines": group["reject_code"]}
        for label, group in sorted(groups.items())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    d = args.directory
    ref_rows = read_csv(d / "equiv_rinex.csv")
    ref_states = read_states(d / "equiv_rinex.stat")
    output = {"comparisons": {}}
    for name in ("full", "compatible", "common", "no_c5q", "no_c1p"):
        if not (d / f"equiv_{name}.csv").exists():
            continue
        output["comparisons"][name] = compare(
            ref_rows, read_csv(d / f"equiv_{name}.csv"), ref_states,
            read_states(d / f"equiv_{name}.stat"))
    output["extra_full"] = extra_signals(
        read_csv(d / "equiv_full3_use.csv"),
        read_csv(d / "equiv_rinex_use.csv"),
        trace_rejects(d / "equiv_full.trace"),
        model_residuals(d / "equiv_full3.trace"))
    output["accepted_full"] = accepted_residual_groups(d / "equiv_full.trace")
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
