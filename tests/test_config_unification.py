"""Regression evidence for the shared PPP configuration (no independent truth).
Inputs are the default K80 PC and replay outputs in gnss_replay/config_unification.
This test checks all printed state lines, not Q counts alone.
"""
import csv
import math
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gnss_replay/config_unification"

def text(path):
    return path.read_text(encoding="utf-8-sig")

def pos_lines(path):
    return [line for line in text(path).splitlines() if line and not line.startswith("%")]

def stat_lines(path):
    return text(path).splitlines()

def test_configuration():
    pc = text(OUT / "k80_pc.log").splitlines()[0].rsplit(",entry=", 1)[0]
    replay = text(OUT / "k80_replay.log").splitlines()[0].rsplit(",entry=", 1)[0]
    assert pc == replay
    assert "NFREQ=4" in pc and "nf=4" in pc
    assert "-DOPPWARM=0" in pc and "-GALE5BPHASE=0" in pc
    assert "-BDSCODEVAR=1" in pc
    assert "NSATQZS=0" in pc
    return re.search(r"signature=([a-f0-9]+)", pc).group(1)

def test_solution():
    current = pos_lines(OUT / "k80_pc.pos")
    prior = pos_lines(ROOT / "gnss_replay/coarse_experiment/k80_final_default.pos")
    assert current == prior, "PC numerical regression"
    pc_state = stat_lines(OUT / "k80_pc.pos.stat")
    assert pc_state == stat_lines(
        ROOT / "gnss_replay/coarse_experiment/k80_final_default.pos.stat")
    replay_state = stat_lines(OUT / "k80_replay.stat")
    assert pc_state == replay_state, "PC/replay printed filter state mismatch"
    # The PC .pos has fewer decimals. Compare after matching its output precision.
    replay = list(csv.DictReader((OUT / "k80_replay.csv").open()))
    state_epochs = [line.split(",")[1:3] for line in pc_state if line.startswith("$POS,")]
    assert len(current) == len(replay) == 1200
    assert len(state_epochs) == len(replay)
    for line, row, state_epoch in zip(current, replay, state_epochs):
        cols = line.split()
        assert int(cols[5]) == int(row["Q"]) and int(cols[6]) == int(row["ns"])
        assert f"{float(row['lat_deg']):.9f}" == cols[2]
        assert f"{float(row['lon_deg']):.9f}" == cols[3]
        assert abs(float(row["h_m"]) - float(cols[4])) <= 0.000051
        assert int(row["week"]) == int(state_epoch[0])
        # CSV labels the raw input epoch; .stat labels sol.time, which SPP can
        # correct for receiver clock offset. Match the same 1 Hz event, rather
        # than demanding identical raw and clock-corrected timestamps.
        assert abs(float(row["tow"]) - float(state_epoch[1])) < 0.5
    counts = dict(Counter(row["Q"] for row in replay))
    assert counts == {"6": 1197, "5": 3}
    return counts, len(pc_state)

if __name__ == "__main__":
    print("profile signature:", test_configuration())
    print("solution Q counts / identical state records:", test_solution())
    print("PASS: PC unchanged; PC and replay match at exported precision.")
