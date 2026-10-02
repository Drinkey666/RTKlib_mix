"""Analyze controlled height experiments without feeding truth to the solver."""
import csv
import json
import math
from pathlib import Path
from collections import Counter
from statistics import mean
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gnss_replay/height_bias_audit"
# Provisional date-specific references; confirm antenna/ellipsoid realization.
REFERENCES = {"0317": (30.759928,103.984162,482.2530),
              "0321": (30.759991,103.984034,482.2020)}

def ecef(llh):
    lat,lon,h=llh;lat=math.radians(lat);lon=math.radians(lon)
    a=6378137.;e2=6.6943799901413165e-3;N=a/math.sqrt(1-e2*math.sin(lat)**2)
    return ((N+h)*math.cos(lat)*math.cos(lon),(N+h)*math.cos(lat)*math.sin(lon),(N*(1-e2)+h)*math.sin(lat))
def enu(row,ref):
    x,y,z=ecef(ref);lat,lon=map(math.radians,ref[:2])
    dx,dy,dz=float(row["x_m"])-x,float(row["y_m"])-y,float(row["z_m"])-z
    return (-math.sin(lon)*dx+math.cos(lon)*dy,
        -math.sin(lat)*math.cos(lon)*dx-math.sin(lat)*math.sin(lon)*dy+math.cos(lat)*dz,
        math.cos(lat)*math.cos(lon)*dx+math.cos(lat)*math.sin(lon)*dy+math.sin(lat)*dz)
def read(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def rms(xs):return math.sqrt(mean(x*x for x in xs))
def correlation(a,b):
    ma,mb=mean(a),mean(b);cov=sum((x-ma)*(y-mb) for x,y in zip(a,b))
    den=math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
    return cov/den if den else None

results=[]
for day,ref in REFERENCES.items():
    baseline=read(OUT/day/"BASE/solution.csv")
    base_bytime={(r["week"],r["tow"]):r for r in baseline}
    for case_dir in sorted((OUT/day).iterdir()):
        if not (case_dir/"solution.csv").exists():continue
        rows=read(case_dir/"solution.csv");diag=read(case_dir/"height_diag.csv")
        diag_bytime={(r["week"],r["tow"]):r for r in diag}
        q6=[r for r in rows if r["Q"]=="6"]
        tail=[r for r in q6 if float(r["tow"])>=float(rows[-1]["tow"])-299.5]
        errors=[enu(r,ref) for r in q6];tail_errors=[enu(r,ref) for r in tail]
        final=diag[-1]
        diffs=[float(r["h_m"])-float(base_bytime[(r["week"],r["tow"])]["h_m"]) for r in q6]
        rec={"day":day,"case":case_dir.name,"epochs":len(rows),"Q":dict(Counter(r["Q"] for r in rows)),
             "U_mean":mean(e[2] for e in errors),"U_RMS":rms([e[2] for e in errors]),
             "tail_U_mean":mean(e[2] for e in tail_errors),"tail_U_RMS":rms([e[2] for e in tail_errors]),
             "tail_E_RMS":rms([e[0] for e in tail_errors]),"tail_N_RMS":rms([e[1] for e in tail_errors]),
             "tail_3D_RMS":math.sqrt(mean(sum(x*x for x in e) for e in tail_errors)),
             "final_U":enu(rows[-1],ref)[2],"final_ZTD":float(final["ZTD_m"]),
             "final_ZHD":float(final["ZHD_m"]),"final_ZWD":float(final["ZWDest_m"]),
             "final_ZWD_prior":float(final["ZWDprior_m"]),"final_corr_U_ZTD":float(final["corr_U_ZTD"]),
             "H_diff_base_RMS":rms(diffs),"vmf_valid_epochs":sum(math.isfinite(float(r["ZHD_m"])) for r in diag),
             "minimum_H":min(float(r["h_m"]) for r in rows),"maximum_H":max(float(r["h_m"]) for r in rows),
             "temporal_H_ZTD_corr":correlation([float(r["h_m"]) for r in q6],[float(diag_bytime[(r["week"],r["tow"])]["ZTD_m"]) for r in q6])}
        # Accepted, final-postfit grouped residuals, not stale ssat caches.
        sigs={}
        text=(case_dir/"ppp.trace").read_text(errors="replace")
        for line in text.splitlines():
            if "$PPP_DIAG_SIG," not in line:continue
            fields=dict(re.findall(r"([A-Za-z0-9_]+)=([^,\n]+)",line))
            key=fields["sys"]+":"+fields["sig"]
            agg=sigs.setdefault(key,{"nP":0,"sumP":0.,"sumP2":0.,"nL":0,"sumL2":0.})
            np,nl=int(fields["nP"]),int(fields["nL"])
            agg["nP"]+=np;agg["sumP"]+=np*float(fields["meanP"]);agg["sumP2"]+=np*float(fields["rmsP"])**2
            agg["nL"]+=nl;agg["sumL2"]+=nl*float(fields["rmsL"])**2
        rec["signal_residuals"]={key:{"nP":v["nP"],"meanP":v["sumP"]/v["nP"] if v["nP"] else None,
             "rmsP":math.sqrt(v["sumP2"]/v["nP"]) if v["nP"] else None,
             "nL":v["nL"],"rmsL":math.sqrt(v["sumL2"]/v["nL"]) if v["nL"] else None}
             for key,v in sigs.items()}
        if (case_dir/"ion_diag.csv").exists():
            ions=read(case_dir/"ion_diag.csv")
            end=max(float(r["tow"]) for r in ions)
            ions=[r for r in ions if float(r["tow"])>=end-299.5]
            rec["tail_ion_diag"]={"sat_epochs":len(ions),
                "state_minus_prior_mean_m":mean(float(r["delta_m"]) for r in ions),
                "state_minus_prior_rms_m":rms([float(r["delta_m"]) for r in ions]),
                "mean_prior_m":mean(float(r["Iprior_m"]) for r in ions),
                "mean_state_m":mean(float(r["Iest_m"]) for r in ions)}
        results.append(rec)
(OUT/"summary.json").write_text(json.dumps(results,indent=2),encoding="utf-8")
for r in results:
    print(r["day"],r["case"],"U(mean/RMS/tailmean/tailRMS):",
          *(f"{r[k]:+.3f}" for k in ("U_mean","U_RMS","tail_U_mean","tail_U_RMS")),
          "tail3D:",f'{r["tail_3D_RMS"]:.3f}',"ZWD:",f'{r["final_ZWD"]:.3f}')
