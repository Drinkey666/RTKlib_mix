"""Read-only inventory and exact-signal BDS B1I/B1C product gate.

Run from the repository root. All times reported here are GPST as stored in the
RINEX, SP3, CLK and Bias-SINEX files. An OSB interval is [start, end).
"""
import bisect
import csv
import datetime as dt
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gnss_replay" / "bds_b1_product_search"
SEARCH_ROOTS = [Path(r"E:\GNSS"), Path(r"E:\RTKLIB_Data"),
                Path(r"E:\collect data"), Path(r"E:\BaiduNetdiskDownload"),
                Path(r"C:\Users\Drinkey\Documents\Codex")]
EXTS = {".bia": "BIAS", ".bsx": "BIAS", ".bias": "BIAS",
        ".atx": "ATX", ".clk": "CLK", ".sp3": "SP3"}
GPS0 = dt.datetime(1980, 1, 6)
CODES = ("C2I", "L2I", "C1D", "C1P", "C1X", "L1D", "L1P", "L1X")
SESSIONS = ("raw_20260922", "rinex_20260321")


def iso(t):
    return t.isoformat(sep=" ", timespec="milliseconds") if t else ""


def bias_time(s):
    y, doy, sec = map(int, s.split(":"))
    return dt.datetime(y, 1, 1) + dt.timedelta(days=doy - 1, seconds=sec)


def epoch(fields):
    return dt.datetime(*map(int, fields[:5]), 0) + dt.timedelta(seconds=float(fields[5]))


def family(path):
    name = path.name.upper()
    match = re.match(r"([A-Z]{3}0MGX(?:RTS|NRT|FIN))_", name)
    return match.group(1) if match else name.split("_")[0]


def center(path):
    product_id = family(path)
    return "WHU" if product_id.startswith(("WUM", "WHU")) else product_id[:3]


def sampling(path):
    match = re.search(r"_[0-9]{2}[SMHD]_([0-9]{2}[SMHD])_", path.name.upper())
    return match.group(1) if match else ""


def inventory_paths():
    paths = []
    for root in SEARCH_ROOTS:
        if root.exists():
            paths.extend(p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in EXTS)
    return sorted(set(paths), key=lambda p: str(p).lower())


def read_bias(path):
    records = defaultdict(list)
    starts, ends, systems = [], [], set()
    with path.open(encoding="ascii", errors="replace") as source:
        for line in source:
            f = line.split()
            if len(f) < 8 or f[0] != "OSB" or not re.fullmatch(r"[GRECJ][0-9]{2}", f[2]):
                continue
            try:
                start, end = bias_time(f[4]), bias_time(f[5])
            except (ValueError, IndexError):
                continue
            if end <= start:
                continue
            systems.add(f[2][0]); starts.append(start); ends.append(end)
            if f[2][0] == "C" and f[3] in CODES:
                records[(f[2], f[3])].append((start, end))
    for intervals in records.values():
        intervals.sort()
    return dict(records=records, start=min(starts) if starts else None,
                end=max(ends) if ends else None, systems=systems)


def read_sp3(path):
    records = defaultdict(list)
    current = None
    all_times, systems = [], set()
    with path.open(encoding="ascii", errors="replace") as source:
        for line in source:
            if line.startswith("*  "):
                try:
                    current = epoch(line[1:].split())
                    all_times.append(current)
                except (ValueError, IndexError):
                    current = None
            elif current and line.startswith("P") and re.match(r"[GRECJ][0-9]{2}", line[1:4]):
                sat = line[1:4]
                systems.add(sat[0])
                try:
                    xyz = [float(line[i:i + 14]) for i in (4, 18, 32)]
                    valid = all(abs(x) < 99999 and x != 0 for x in xyz)
                except ValueError:
                    valid = False
                if valid:
                    records[sat].append(current)
    return dict(records=records, start=min(all_times) if all_times else None,
                end=max(all_times) if all_times else None, systems=systems)


def read_clk(path):
    records = defaultdict(list)
    all_times, systems = [], set()
    with path.open(encoding="ascii", errors="replace") as source:
        for line in source:
            if not line.startswith("AS "):
                continue
            f = line.split()
            if len(f) < 10 or not re.fullmatch(r"[GRECJ][0-9]{2}", f[1]):
                continue
            try:
                t = epoch(f[2:8])
                float(f[9])
            except (ValueError, IndexError):
                continue
            records[f[1]].append(t); all_times.append(t); systems.add(f[1][0])
    return dict(records=records, start=min(all_times) if all_times else None,
                end=max(all_times) if all_times else None, systems=systems)


def atx_date(line):
    f = line[:60].split()
    if len(f) < 6:
        return None
    try:
        return epoch(f)
    except ValueError:
        return None


def read_atx(path):
    blocks = defaultdict(list)
    systems = set()
    sat = None; freqs = set(); start = None; end = None
    def save():
        if sat:
            blocks[sat].append((start, end, frozenset(freqs)))
    with path.open(encoding="ascii", errors="replace") as source:
        for line in source:
            label = line[60:].strip()
            if label == "START OF ANTENNA":
                sat = None; freqs = set(); start = None; end = None
            elif label == "TYPE / SERIAL NO":
                candidate = line[20:40].strip()
                if re.fullmatch(r"[GRECJ][0-9]{2,3}", candidate):
                    systems.add(candidate[0])
                sat = candidate if re.fullmatch(r"C[0-9]{2}", candidate) else None
            elif label == "VALID FROM" and sat:
                start = atx_date(line)
            elif label == "VALID UNTIL" and sat:
                end = atx_date(line)
            elif label == "START OF FREQUENCY" and sat:
                freqs.add(line[:10].strip())
            elif label == "END OF ANTENNA":
                save(); sat = None
    starts = [a for values in blocks.values() for a, _, _ in values if a]
    ends = [b for values in blocks.values() for _, b, _ in values if b]
    return dict(records=blocks, start=min(starts) if starts else None,
                end=max(ends) if ends else None, systems=systems,
                open_ended=any(b is None for values in blocks.values() for _, b, _ in values))


def osb_at(product, sat, code, t):
    return next(((a, b) for a, b in product["records"].get((sat, code), ())
                 if a <= t < b), None)


def atx_at(product, sat, freq, t):
    return any((a is None or a <= t) and (b is None or t < b) and freq in fs
               for a, b, fs in product["records"].get(sat, ()))


def bracketing(product, sat, t, max_gap):
    times = product["records"].get(sat, ())
    index = bisect.bisect_left(times, t)
    if index < len(times) and times[index] == t:
        return True
    return 0 < index < len(times) and (times[index] - times[index - 1]).total_seconds() <= max_gap


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = inventory_paths()
    parsed = {}
    inventory = []
    for p in paths:
        kind = EXTS[p.suffix.lower()]
        try:
            data = {"BIAS": read_bias, "SP3": read_sp3, "CLK": read_clk,
                    "ATX": read_atx}[kind](p)
            error = ""
        except (OSError, UnicodeError, ValueError) as exc:
            data = dict(records={}, start=None, end=None, systems=set())
            error = str(exc)
        parsed[p] = data
        inventory.append(dict(path=str(p), filename=p.name, product_type=kind,
                              analysis_center=center(p), product_family=family(p),
                              start_time=iso(data["start"]),
                              end_time="OPEN" if kind == "ATX" and data.get("open_ended") else iso(data["end"]),
                              sampling=sampling(p), GNSS_systems="".join(sorted(data["systems"])),
                              bytes=p.stat().st_size, parse_error=error))
    write_csv(OUT / "product_inventory.csv", inventory,
              ["path", "filename", "product_type", "analysis_center", "product_family",
               "start_time", "end_time", "sampling", "GNSS_systems", "bytes", "parse_error"])

    bias_catalog = []
    for path in paths:
        if EXTS[path.suffix.lower()] != "BIAS":
            continue
        records = parsed[path]["records"]
        sats = sorted({sat for sat, _ in records})
        for sat in sats:
            for code in CODES:
                intervals = records.get((sat, code), ())
                gaps = [(intervals[i][0] - intervals[i - 1][1]).total_seconds()
                        for i in range(1, len(intervals))]
                bias_catalog.append(dict(path=str(path), PRN=sat, signal=code,
                                         intervals=len(intervals),
                                         first_valid=iso(intervals[0][0]) if intervals else "",
                                         last_valid=iso(intervals[-1][1]) if intervals else "",
                                         largest_internal_gap_s=max([0, *gaps])))
    write_csv(OUT / "bias_signal_catalog.csv", bias_catalog,
              ["path", "PRN", "signal", "intervals", "first_valid", "last_valid",
               "largest_internal_gap_s"])

    with (ROOT / "gnss_replay/all_signal_usage/b1i_b1c_coexist.csv").open(newline="", encoding="utf-8") as source:
        coexist = list(csv.DictReader(source))
    groups = {s: [] for s in SESSIONS}
    for row in coexist:
        if row["session"] in groups:
            row["epoch"] = GPS0 + dt.timedelta(weeks=int(row["week"]), seconds=float(row["tow"]))
            groups[row["session"]].append(row)

    bias_candidates = {}
    orbit_candidates = {}
    clk_candidates = {}
    atx_candidates = [p for p in paths if EXTS[p.suffix.lower()] == "ATX"]
    for session, rows in groups.items():
        first, last = min(r["epoch"] for r in rows), max(r["epoch"] for r in rows)
        def overlapping(p):
            d = parsed[p]
            return d["start"] is not None and d["end"] is not None and d["start"] <= last and d["end"] >= first
        bias_candidates[session] = [p for p in paths if EXTS[p.suffix.lower()] == "BIAS" and overlapping(p)]
        orbit_candidates[session] = [p for p in paths if EXTS[p.suffix.lower()] == "SP3" and overlapping(p)]
        clk_candidates[session] = [p for p in paths if EXTS[p.suffix.lower()] == "CLK" and overlapping(p)]

    # Score each single bias product against each session. No time extrapolation or cross-file stitching.
    bias_scores = []
    for session, rows in groups.items():
        for p in bias_candidates[session]:
            d = parsed[p]; counts = Counter()
            for row in rows:
                sat, t = row["sat"], row["epoch"]
                present = {c: bool(osb_at(d, sat, c, t)) for c in CODES}
                counts.update({c: int(present[c]) for c in CODES})
                counts["pair"] += int(all(present[c] for c in ("C2I", "L2I", "C1P", "L1P")))
            bias_scores.append(dict(session=session, bias_path=str(p), family=family(p),
                                    **{f"{c}_epochs": counts[c] for c in CODES},
                                    pair_complete_epochs=counts["pair"]))
    write_csv(OUT / "bias_candidate_coverage.csv", bias_scores,
              ["session", "bias_path", "family"] + [f"{c}_epochs" for c in CODES] + ["pair_complete_epochs"])

    # Choose a single best product per type. Exact family alignment outranks coverage ties.
    choices = {}
    for session, rows in groups.items():
        biases = sorted((r for r in bias_scores if r["session"] == session),
                        key=lambda r: (r["pair_complete_epochs"],
                                       r["bias_path"].startswith(r"E:\RTKLIB_Data"),
                                       "FIN" in r["family"]), reverse=True)
        selected_bias = Path(biases[0]["bias_path"]) if biases else None
        selected_family = family(selected_bias) if selected_bias else ""
        def best_product(candidates, kind):
            scores = []
            for p in candidates:
                d = parsed[p]
                max_gap = 900 if kind == "SP3" else (360 if sampling(p) == "05M" else 120)
                count = sum(bracketing(d, r["sat"], r["epoch"], max_gap) for r in rows)
                scores.append((count, family(p) == selected_family, p))
            return max(scores, key=lambda item: (item[0], item[1], str(item[2])))[2] if scores else None
        orbit = best_product(orbit_candidates[session], "SP3")
        clk = best_product(clk_candidates[session], "CLK")
        atx_scores = [(sum(atx_at(parsed[p], r["sat"], f, r["epoch"]) for r in rows for f in ("C01", "C02")), p)
                      for p in atx_candidates]
        atx = max(atx_scores, key=lambda item: (item[0],
                   str(item[1]).lower() == r"e:\rtklib_data\tables\igs20.atx",
                   str(item[1])))[1] if atx_scores else None
        choices[session] = dict(bias=selected_bias, sp3=orbit, clk=clk, atx=atx)

    matrix = []; pair_detail = []; summaries = {}
    for session, rows in groups.items():
        choice = choices[session]
        by_sat = defaultdict(list)
        counts = Counter()
        for row in rows:
            sat, t = row["sat"], row["epoch"]
            per_sig = {}
            for sig, code, phase, freq in (("B1I_2I", "C2I", "L2I", "C02"),
                                           ("B1C_1P", "C1P", "L1P", "C01")):
                b = parsed[choice["bias"]] if choice["bias"] else {"records": {}}
                ci, li = osb_at(b, sat, code, t), osb_at(b, sat, phase, t)
                ao = atx_at(parsed[choice["atx"]], sat, freq, t) if choice["atx"] else False
                so = bracketing(parsed[choice["sp3"]], sat, t, 900) if choice["sp3"] else False
                cg = 360 if choice["clk"] and sampling(choice["clk"]) == "05M" else 120
                co = bracketing(parsed[choice["clk"]], sat, t, cg) if choice["clk"] else False
                per_sig[sig] = dict(code=bool(ci), phase=bool(li), atx=ao, sp3=so, clk=co,
                                    osb_start=max(ci[0], li[0]) if ci and li else None,
                                    osb_end=min(ci[1], li[1]) if ci and li else None)
                by_sat[(sat, sig)].append(per_sig[sig])
            a, b = per_sig["B1I_2I"], per_sig["B1C_1P"]
            bias_ok = all(q["code"] and q["phase"] for q in (a, b))
            atx_ok = a["atx"] and b["atx"]
            orbit_ok = a["sp3"] and b["sp3"]
            clock_ok = a["clk"] and b["clk"]
            full = bias_ok and atx_ok and orbit_ok and clock_ok
            counts["raw"] += 1; counts["bias"] += bias_ok; counts["atx"] += atx_ok
            counts["orbit"] += orbit_ok; counts["clock"] += clock_ok; counts["full"] += full
            pair_detail.append(dict(session=session, week=row["week"], tow=row["tow"], PRN=sat,
                                    B1I_CODE_OK=int(a["code"]), B1I_PHASE_OK=int(a["phase"]),
                                    B1C_CODE_OK=int(b["code"]), B1C_PHASE_OK=int(b["phase"]),
                                    B1I_ATX_OK=int(a["atx"]), B1C_ATX_OK=int(b["atx"]),
                                    orbit_ok=int(orbit_ok), clock_ok=int(clock_ok),
                                    BIAS_complete=int(bias_ok), ATX_complete=int(atx_ok),
                                    FULL_MODEL_COMPLETE=int(full)))
        for (sat, sig), items in sorted(by_sat.items()):
            code = "C2I" if sig == "B1I_2I" else "C1P"
            phase = "L2I" if sig == "B1I_2I" else "L1P"
            all_ok = all(q["code"] and q["phase"] and q["atx"] and q["sp3"] and q["clk"] for q in items)
            intervals = [(q["osb_start"], q["osb_end"]) for q in items if q["osb_start"]]
            matrix.append(dict(session=session, PRN=sat, signal=sig, observed_pairs=len(items),
                               code_osb=f"{code}:{sum(q['code'] for q in items)}/{len(items)}",
                               phase_osb=f"{phase}:{sum(q['phase'] for q in items)}/{len(items)}",
                               osb_start=iso(min((a for a, _ in intervals), default=None)),
                               osb_end=iso(max((b for _, b in intervals), default=None)),
                               ATX_frequency=("C02" if sig == "B1I_2I" else "C01"),
                               atx_ok=sum(q["atx"] for q in items), orbit_ok=sum(q["sp3"] for q in items),
                               clock_ok=sum(q["clk"] for q in items),
                               SP3=str(choice["sp3"] or ""), CLK=str(choice["clk"] or ""),
                               BIAS=str(choice["bias"] or ""), ATX_file=str(choice["atx"] or ""),
                               product_family=family(choice["bias"]) if choice["bias"] else "",
                               combination_status=("SAME_FAMILY" if choice["bias"] and choice["sp3"] and choice["clk"]
                                                   and len({family(choice[k]) for k in ("bias", "sp3", "clk")}) == 1
                                                   else "CANDIDATE_UNVERIFIED"),
                               complete=int(all_ok), complete_pairs=sum(q["code"] and q["phase"] and q["atx"] and q["sp3"] and q["clk"] for q in items)))
        summaries[session] = dict(counts, first=iso(min(r["epoch"] for r in rows)),
                                  last=iso(max(r["epoch"] for r in rows)),
                                  chosen={k: str(v or "") for k, v in choice.items()},
                                  product_family_status=("SAME_FAMILY" if all(choice.values()) and
                                    len({family(choice[k]) for k in ("bias", "sp3", "clk")}) == 1
                                    else "CANDIDATE_UNVERIFIED"))
    write_csv(OUT / "BDS_B1_PRODUCT_MATRIX.csv", matrix,
              ["session", "PRN", "signal", "observed_pairs", "code_osb", "phase_osb",
               "osb_start", "osb_end", "ATX_frequency", "atx_ok", "SP3", "CLK", "BIAS",
               "ATX_file", "orbit_ok", "clock_ok", "product_family", "combination_status",
               "complete", "complete_pairs"])
    write_csv(OUT / "pair_completeness.csv", pair_detail,
              ["session", "week", "tow", "PRN", "B1I_CODE_OK", "B1I_PHASE_OK",
               "B1C_CODE_OK", "B1C_PHASE_OK", "B1I_ATX_OK", "B1C_ATX_OK",
               "orbit_ok", "clock_ok", "BIAS_complete", "ATX_complete", "FULL_MODEL_COMPLETE"])

    antenna_rows = []
    for session, rows in groups.items():
        bysat = defaultdict(list)
        for r in rows: bysat[r["sat"]].append(r["epoch"])
        for path in atx_candidates:
            d = parsed[path]
            for sat, times in sorted(bysat.items()):
                antenna_rows.append(dict(session=session, PRN=sat,
                    **{f"{f}_present": int(all(atx_at(d, sat, f, t) for t in times))
                       for f in ("C01", "C02", "C05", "C07")}, ATX_file=str(path)))
    write_csv(OUT / "atx_frequency_matrix.csv", antenna_rows,
              ["session", "PRN", "C01_present", "C02_present", "C05_present", "C07_present", "ATX_file"])

    # Exact signal availability by each individual candidate bias product, for every observed PRN.
    exact_rows = []
    for session, rows in groups.items():
        bysat = defaultdict(list)
        for r in rows: bysat[r["sat"]].append(r["epoch"])
        for path in bias_candidates[session]:
            d = parsed[path]
            for sat, times in sorted(bysat.items()):
                for code in CODES:
                    matches = [osb_at(d, sat, code, t) for t in times]
                    valid = [x for x in matches if x]
                    exact_rows.append(dict(session=session, bias_path=str(path), PRN=sat,
                        signal=code, covered_epochs=len(valid), observed_epochs=len(times),
                        full_session=int(len(valid) == len(times)),
                        first_valid=iso(min((x[0] for x in valid), default=None)),
                        last_valid=iso(max((x[1] for x in valid), default=None))))
    write_csv(OUT / "bias_exact_signal_coverage.csv", exact_rows,
              ["session", "bias_path", "PRN", "signal", "covered_epochs", "observed_epochs",
               "full_session", "first_valid", "last_valid"])

    combo_rows = []
    for session, rows in groups.items():
        def coverage(path, kind):
            d = parsed[path]
            if kind == "BIAS":
                return [all(osb_at(d, r["sat"], code, r["epoch"])
                            for code in ("C2I", "L2I", "C1P", "L1P")) for r in rows]
            if kind == "ATX":
                return [atx_at(d, r["sat"], "C01", r["epoch"]) and
                        atx_at(d, r["sat"], "C02", r["epoch"]) for r in rows]
            gap = 900 if kind == "SP3" else (360 if sampling(path) == "05M" else 120)
            return [bracketing(d, r["sat"], r["epoch"], gap) for r in rows]
        flags = {p: coverage(p, EXTS[p.suffix.lower()]) for p in
                 bias_candidates[session] + orbit_candidates[session] +
                 clk_candidates[session] + atx_candidates}
        for b in bias_candidates[session]:
            for s in orbit_candidates[session]:
                for c in clk_candidates[session]:
                    for a in atx_candidates:
                        matches = sum(all(x) for x in zip(flags[b], flags[s], flags[c], flags[a]))
                        coherent = len({family(b), family(s), family(c)}) == 1
                        combo_rows.append(dict(session=session, BIAS=str(b), SP3=str(s), CLK=str(c),
                                               ATX=str(a), BIAS_center=center(b),
                                               SP3_center=center(s), CLK_center=center(c),
                                               BIAS_family=family(b), SP3_family=family(s),
                                               CLK_family=family(c),
                                               compatibility="SAME_FAMILY" if coherent else "CANDIDATE_UNVERIFIED",
                                               FULL_MODEL_COMPLETE_pairs=matches))
    write_csv(OUT / "product_combinations.csv", combo_rows,
              ["session", "BIAS", "SP3", "CLK", "ATX", "BIAS_center", "SP3_center",
               "CLK_center", "BIAS_family", "SP3_family", "CLK_family",
               "compatibility", "FULL_MODEL_COMPLETE_pairs"])
    write_csv(OUT / "pair_counts.csv", [dict(session=s, RAW_coexist=v.get("raw", 0),
              BIAS_complete=v.get("bias", 0), ATX_complete=v.get("atx", 0),
              orbit_complete=v.get("orbit", 0), clock_complete=v.get("clock", 0),
              FULL_MODEL_COMPLETE=v.get("full", 0), combination_status=v["product_family_status"])
              for s, v in summaries.items()],
              ["session", "RAW_coexist", "BIAS_complete", "ATX_complete",
               "orbit_complete", "clock_complete", "FULL_MODEL_COMPLETE", "combination_status"])
    (OUT / "summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    print(json.dumps(dict(inventory_files=len(paths), sessions=summaries), indent=2))


if __name__ == "__main__":
    main()
