"""Compare three gnss_replay PPP runs without changing the solver.

Run after txt_fast_solution.csv, txt_realtime_solution.csv and
rinex_solution.csv (plus matching .stat and .stderr.log files) exist.
"""

import argparse
import csv
import hashlib
import math
import re
from pathlib import Path

from analyze_equiv import read_states

NAMES = ("txt_fast", "txt_realtime", "rinex")
PAIRS = (("txt_fast", "txt_realtime"),
         ("rinex", "txt_realtime"), ("rinex", "txt_fast"))
RCB = ("GPS_C5Q", "GAL_C5Q", "BDS_C7I", "BDS_C5P")


def stamp(row):
    return int(row["week"]), round(float(row["tow"]) * 1000)


def read_solutions(path):
    rows = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            epoch = stamp(row)
            if epoch in rows:
                raise ValueError(f"Duplicate epoch {epoch} in {path}")
            rows[epoch] = row
    if not rows:
        raise ValueError(f"No solution epochs in {path}")
    return rows


def metric(values):
    if not values:
        return None
    return (sum(values) / len(values),
            math.sqrt(sum(x * x for x in values) / len(values)),
            max(abs(x) for x in values))


def fmt(values, digits=6):
    m = metric(values)
    return "n/a" if m is None else "/".join(f"{v:.{digits}f}" for v in m)


def state_at(states, epoch):
    # rtkoutstat prints three decimals; analyze_equiv stores two decimals.
    return states.get((epoch[0], round(epoch[1] / 1000, 2)), {})


def enu_rotation(reference):
    lat = math.radians(float(reference["lat_deg"]))
    lon = math.radians(float(reference["lon_deg"]))
    s, c = math.sin, math.cos
    return ((-s(lon), c(lon), 0.0),
            (-s(lat) * c(lon), -s(lat) * s(lon), c(lat)),
            (c(lat) * c(lon), c(lat) * s(lon), s(lat)))


def perf(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"PERF,epochs=(\d+),mean_ms=([\d.]+),p95_ms=([\d.]+),"
                      r"max_ms=([\d.]+),deadline_miss_gt_1000ms=(\d+)", text)
    if not match:
        raise ValueError(f"Missing PERF line in {path}")
    return tuple(float(x) for x in match.groups())


def write_enriched(root, name, rows, states):
    """Keep the original solver CSV intact and publish joined state columns."""
    extra = ("ztd_m", "gps_clock_ns", "gal_isb_ns", "bds_isb_ns",
             *("rcb_" + signal + "_m" for signal in RCB))
    fields = (*next(iter(rows.values())).keys(), *extra)
    path = root / f"{name}_solution_with_states.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for epoch in sorted(rows):
            row = dict(rows[epoch])
            state = state_at(states, epoch)
            if "ztd" in state:
                row["ztd_m"] = state["ztd"]
            if "clock" in state:
                for index, field in ((0, "gps_clock_ns"), (2, "gal_isb_ns"),
                                     (3, "bds_isb_ns")):
                    row[field] = state["clock"][index]
            for signal in RCB:
                if signal in state.get("rcb", {}):
                    row["rcb_" + signal + "_m"] = state["rcb"][signal]
            writer.writerow(row)


def compare(root, base_name, other_name, solutions, states, anchor):
    base, other = solutions[base_name], solutions[other_name]
    common = sorted(base.keys() & other.keys())
    q6 = [k for k in common if base[k]["Q"] == other[k]["Q"] == "6"]
    if not common:
        raise ValueError(f"No common epochs: {base_name}, {other_name}")
    rotation = enu_rotation(anchor)
    values = {name: [] for name in ("dX", "dY", "dZ", "dE", "dN", "dU",
                                     "d3D", "dNs", "dZTD", "dGPSclock_ns",
                                     "dGAL_ISB_ns", "dBDS_ISB_ns",
                                     *("dRCB_" + x for x in RCB))}
    all_solver_equal = True
    state_differences = 0
    status_mismatch = 0
    detail = root / f"{other_name}_minus_{base_name}.csv"
    fieldnames = ("week", "tow", "base_Q", "other_Q", "base_ns", "other_ns",
                  "dX", "dY", "dZ", "dE", "dN", "dU", "d3D", "dZTD",
                  "dGPSclock_ns", "dGAL_ISB_ns", "dBDS_ISB_ns",
                  *("dRCB_" + x for x in RCB))
    with detail.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for epoch in common:
            a, b = base[epoch], other[epoch]
            all_solver_equal &= all(a[f] == b[f] for f in a if f != "processing_ms")
            status_mismatch += a["Q"] != b["Q"]
            xyz = [float(b[f]) - float(a[f]) for f in ("x_m", "y_m", "z_m")]
            enu = [sum(r[j] * xyz[j] for j in range(3)) for r in rotation]
            row = {"week": epoch[0], "tow": f"{epoch[1] / 1000:.3f}",
                   "base_Q": a["Q"], "other_Q": b["Q"],
                   "base_ns": a["ns"], "other_ns": b["ns"]}
            row.update(dict(zip(("dX", "dY", "dZ", "dE", "dN", "dU"),
                                xyz + enu)))
            row["d3D"] = math.sqrt(sum(x * x for x in xyz))
            sa, sb = state_at(states[base_name], epoch), state_at(states[other_name], epoch)
            if "ztd" in sa and "ztd" in sb:
                row["dZTD"] = sb["ztd"] - sa["ztd"]
            if "clock" in sa and "clock" in sb:
                for idx, label in ((0, "dGPSclock_ns"), (2, "dGAL_ISB_ns"),
                                   (3, "dBDS_ISB_ns")):
                    row[label] = sb["clock"][idx] - sa["clock"][idx]
            for name in RCB:
                if name in sa.get("rcb", {}) and name in sb.get("rcb", {}):
                    row["dRCB_" + name] = sb["rcb"][name] - sa["rcb"][name]
            if epoch in q6:
                for label in values:
                    if label == "dNs":
                        values[label].append(int(b["ns"]) - int(a["ns"]))
                    elif label in row:
                        values[label].append(row[label])
            if sa != sb:
                state_differences += 1
            writer.writerow(row)
    final = q6[-1] if q6 else None
    a, b = (base[final], other[final]) if final else (None, None)
    xyz = [float(b[f]) - float(a[f]) for f in ("x_m", "y_m", "z_m")] if final else None
    final_enu = [sum(r[j] * xyz[j] for j in range(3)) for r in rotation] if final else None
    return {"base": base_name, "other": other_name, "common": len(common),
            "common_q6": len(q6), "base_only": len(base.keys() - other.keys()),
            "other_only": len(other.keys() - base.keys()),
            "status_mismatch": status_mismatch,
            "all_solver_fields_equal": all_solver_equal,
            "full_state_file_equal": (hashlib.sha256((root / f"{base_name}.stat").read_bytes()).digest() ==
                                      hashlib.sha256((root / f"{other_name}.stat").read_bytes()).digest()),
            "state_differences": state_differences,
            "values": values, "final_epoch": final,
            "final_enu": final_enu, "final_xyz": xyz,
            "final_height": float(b["h_m"]) - float(a["h_m"]) if final else None,
            "detail": detail.name}


def write_report(root, solutions, comparisons, perfs):
    a, b, c = comparisons
    lines = ["# 模拟实时 TXT 与 RINEX 后处理定位对比", "",
             "## 方法与边界", "",
             "同一个 Release `gnss_replay.exe` 分别运行 TXT_FAST、TXT_REALTIME、RINEX。"
             "三组共用程序内同一个 `configure_ppp()` 及同一组 NAV/SP3/CLK/BIA/IONEX/VMF3/ATX；"
             "每一组仅在启动时初始化一次 `rtk_t`，逐历元调用一次 `rtkpos()`。"
             "TXT_REALTIME 使用单调时钟按原始 GPST 差值调度，不在每历元固定 `sleep(1)`。",
             "坐标差均为“后一组－前一组”；ENU 使用首个 RINEX Q6 坐标的"
             "纬经度构造一套固定参考轴，不使用独立真值。"
             "状态取 `rtkoutstat()` 的逐历元 `$CLK/$TROP/$RCB`；"
             "缺失状态不按零填充。耗时列不参与滤波等价性判断。", "",
             "解算 CSV 中 `Q` 即 `solution_status`，`x_m/y_m/z_m` 即 ECEF，"
             "`lat_deg/lon_deg/h_m` 即纬经高，`std_x_m/std_y_m/std_z_m` "
             "即 `sd_x/sd_y/sd_z`。ZTD、钟差、ISB、RCB 保存在同名 `.stat`，"
             "并汇入下方逐历元差值 CSV。附加状态列的解算文件为 "
             "[TXT_REALTIME](txt_realtime_solution_with_states.csv)、"
             "[TXT_FAST](txt_fast_solution_with_states.csv)、"
             "[RINEX](rinex_solution_with_states.csv)；原始解算 CSV 保持不变。", "",
             "## 历元与性能", "",
             "| 路径 | total epochs | Q5 | Q6 | mean/P95/max (ms) | >1000 ms deadline miss |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for name in NAMES:
        rows = solutions[name]
        p = perfs[name]
        lines.append(f"| {name.upper()} | {len(rows)} | "
                     f"{sum(r['Q'] == '5' for r in rows.values())} | "
                     f"{sum(r['Q'] == '6' for r in rows.values())} | "
                     f"{p[1]:.3f}/{p[2]:.3f}/{p[3]:.3f} | {int(p[4])} |")
    lines += ["", "| 比较（后－前） | matched | 共同 Q6 | 仅前/仅后 | "
              "Q 不同 | E RMS | N RMS | U RMS | 3D RMS/max |",
              "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for result in comparisons:
        v = result["values"]
        rms = lambda name: metric(v[name])[1] if v[name] else float("nan")
        max3 = metric(v["d3D"])[2] if v["d3D"] else float("nan")
        lines.append(f"| {result['other'].upper()} − {result['base'].upper()} | "
                     f"{result['common']} | {result['common_q6']} | "
                     f"{result['base_only']}/{result['other_only']} | "
                     f"{result['status_mismatch']} | {rms('dE'):.6f} | "
                     f"{rms('dN'):.6f} | {rms('dU'):.6f} | "
                     f"{rms('d3D'):.6f}/{max3:.6f} |")
    lines += ["", "共同 Q6 历元详细统计，clock/ISB 单位为 ns，dNs 为卫星数，其他为 m；"
              "每格为 mean/RMS/max absolute：", ""]
    for result in comparisons:
        lines += [f"### {result['other'].upper()} − {result['base'].upper()}", "",
                  "| 指标 | mean/RMS/max absolute |", "| --- | ---: |"]
        for label in ("dX", "dY", "dZ", "dE", "dN", "dU", "d3D",
                      "dNs", "dZTD", "dGPSclock_ns", "dGAL_ISB_ns",
                      "dBDS_ISB_ns", *("dRCB_" + x for x in RCB)):
            lines.append(f"| {label} | {fmt(result['values'][label])} |")
        if result["final_epoch"]:
            e, n, u = result["final_enu"]
            x, y, z = result["final_xyz"]
            lines += ["", f"最后共同 Q6 历元 "
                      f"{result['final_epoch'][0]}:{result['final_epoch'][1]/1000:.3f}："
                      f"dE/dN/dU = {e:+.6f}/{n:+.6f}/{u:+.6f} m；"
                      f"dX/dY/dZ = {x:+.6f}/{y:+.6f}/{z:+.6f} m；"
                      f"高程差 {result['final_height']:+.6f} m。"]
        else:
            lines += ["", "没有共同 Q6 历元；最终坐标差不适用。"]
        lines += [f"逐历元明细：[ {result['detail']} ]({result['detail']})。", ""]
    no_loss = all(x["base_only"] == x["other_only"] == 0 for x in comparisons)
    fast_equal = (a["base_only"] == a["other_only"] == 0 and
                  a["all_solver_fields_equal"] and a["full_state_file_equal"])
    lines += ["## 判定", "",
              f"- 实时回放是否丢历元：{'否' if no_loss else '是，见上表仅前/仅后计数'}。",
              f"- TXT_FAST 与 TXT_REALTIME 解算字段及状态是否一致："
              f"{'是' if fast_equal else '否，需检查调度、历元顺序与状态文件'}。",
              f"- 两次 TXT 的完整状态文件逐字节一致：{'是' if a['full_state_file_equal'] else '否'}；"
              f"所选 ZTD/CLK/RCB/ION 状态有差异的历元数：{a['state_differences']}。",
              "- TXT_REALTIME 与 RINEX 的差异只证明两条入口的相对一致性，"
              "不证明相对于独立真值的绝对厘米级精度。",
              "- 输入产品及可执行文件的路径和 SHA-256 见 [inputs.json](inputs.json)；"
              "此实验不修改现有 PPP 参数。", ""]
    (root / "REALTIME_VS_RINEX_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    root = args.directory
    solutions = {name: read_solutions(root / f"{name}_solution.csv") for name in NAMES}
    reference = solutions["rinex"]
    anchor = next((reference[k] for k in sorted(reference) if reference[k]["Q"] == "6"),
                  reference[min(reference)])
    states = {name: read_states(root / f"{name}.stat") for name in NAMES}
    perfs = {name: perf(root / f"{name}.stderr.log") for name in NAMES}
    for name in NAMES:
        write_enriched(root, name, solutions[name], states[name])
    comparisons = [compare(root, first, second, solutions, states, anchor)
                   for first, second in PAIRS]
    write_report(root, solutions, comparisons, perfs)
    print(root / "REALTIME_VS_RINEX_REPORT.md")


if __name__ == "__main__":
    main()
