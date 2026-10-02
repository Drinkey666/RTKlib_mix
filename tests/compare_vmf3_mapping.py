"""Compare the deployed mapping formula with TU Wien's VMF3 coefficients.
Official .m files are read as coefficient DATA, never executed. This is a
standalone audit, not a production model substitution.
"""
import csv
import json
import math
import re
from pathlib import Path
from datetime import datetime,timedelta

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"gnss_replay/height_bias_audit"
source=(OUT/"reference/vmf3_ht.m").read_text(encoding="latin1")
coeff={}
for name in ("bh","bw","ch","cw"):
    for mode in ("anm","bnm"):
        blob=re.search(r"^"+mode+"_"+name+r"\s*=\s*\[(.*?)\];",source,re.M|re.S).group(1)
        values=[[float(x.replace("D","e").replace("d","e")) for x in row.split()] for row in blob.split(";")]
        assert len(values)==91 and all(len(row)==5 for row in values)
        coeff[mode,name]=values

def bc(lat,lon,doy):
    lat,lon=map(math.radians,(lat,lon));x=math.cos(lat)*math.cos(lon);y=math.cos(lat)*math.sin(lon);z=math.sin(lat)
    V=[[0.]*13 for _ in range(13)];W=[[0.]*13 for _ in range(13)]
    V[0][0]=1.;V[1][0]=z
    for n in range(2,13):V[n][0]=((2*n-1)*z*V[n-1][0]-(n-1)*V[n-2][0])/n
    for m in range(1,13):
        V[m][m]=(2*m-1)*(x*V[m-1][m-1]-y*W[m-1][m-1])
        W[m][m]=(2*m-1)*(x*W[m-1][m-1]+y*V[m-1][m-1])
        if m<12:V[m+1][m]=(2*m+1)*z*V[m][m];W[m+1][m]=(2*m+1)*z*W[m][m]
        for n in range(m+2,13):
            V[n][m]=((2*n-1)*z*V[n-1][m]-(n+m-1)*V[n-2][m])/(n-m)
            W[n][m]=((2*n-1)*z*W[n-1][m]-(n+m-1)*W[n-2][m])/(n-m)
    seasonal=(1,math.cos(doy/365.25*2*math.pi),math.sin(doy/365.25*2*math.pi),
              math.cos(doy/365.25*4*math.pi),math.sin(doy/365.25*4*math.pi))
    out={}
    for name in ("bh","bw","ch","cw"):
        out[name]=0.;i=0
        for n in range(13):
            for m in range(n+1):
                out[name]+=sum((coeff["anm",name][i][k]*V[n][m]+coeff["bnm",name][i][k]*W[n][m])*seasonal[k] for k in range(5))
                i+=1
    return out

def mapf(el,a,b,c):
    s=math.sin(math.radians(el))
    return (1+a/(1+b/(1+c)))/(s+a/(s+b/(s+c)))
def grid(path):
    return [list(map(float,l.split())) for l in path.read_text().splitlines() if l and not l.startswith("!")]
def rows(path):
    with path.open() as f:return list(csv.DictReader(f))

result=[]
for day in ("0317","0321"):
    diag=rows(OUT/day/"BASE/height_diag.csv")
    inputs=json.loads((OUT/day/"BASE/inputs.json").read_text())["products"]
    g0,g1=grid(Path(inputs["vmf0"])),grid(Path(inputs["vmf6"]))
    for row in (diag[0],diag[len(diag)//2],diag[-1]):
        lat,lon,h=(float(row[k]) for k in ("lat_deg","lon_deg","h_m"))
        utc=datetime(1980,1,6)+timedelta(weeks=int(row["week"]),seconds=float(row["tow"])-18)
        doy=utc.timetuple().tm_yday+(utc.hour*3600+utc.minute*60+utc.second+utc.microsecond/1e6)/86400
        frac=(utc.hour*3600+utc.minute*60+utc.second+utc.microsecond/1e6)/21600
        ii=int(math.floor((87.5-lat)/5));jj=int(math.floor((lon-2.5)/5))
        u=(87.5-lat)/5-ii;v=(lon-2.5)/5-jj
        corners=((ii,jj,(1-u)*(1-v)),(ii,jj+1,(1-u)*v),(ii+1,jj,u*(1-v)),(ii+1,jj+1,u*v))
        data=[]
        for i,j,weight in corners:
            k=i*72+j
            ah=g0[k][2]+frac*(g1[k][2]-g0[k][2])
            aw=g0[k][3]+frac*(g1[k][3]-g0[k][3])
            data.append((weight,ah,aw,bc(g0[k][0],g0[k][1],doy)))
        for el in (15,20,30,45,60,90):
            correction=(1/math.sin(math.radians(el))-mapf(el,2.53e-5,5.49e-3,1.14e-3))*h/1000
            hydro=sum(w*mapf(el,ah,c["bh"],c["ch"]) for w,ah,aw,c in data)+correction
            wet=sum(w*mapf(el,aw,c["bw"],c["cw"]) for w,ah,aw,c in data)
            deployed_h=mapf(el,float(row["ah"]),0.0029,0.0620)+correction
            deployed_w=mapf(el,float(row["aw"]),0.00146,0.04391)
            result.append({"day":day,"tow":float(row["tow"]),"el":el,
                "dmfh":deployed_h-hydro,"dmfw":deployed_w-wet,
                "dslant_prior_m":(deployed_h-hydro)*float(row["ZHD_m"])+(deployed_w-wet)*float(row["ZWDprior_m"]),
                "dslant_estimated_m":(deployed_h-hydro)*float(row["ZHD_m"])+(deployed_w-wet)*float(row["ZWDest_m"])})
(OUT/"vmf_mapping_comparison.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
for day in ("0317","0321"):
    data=[r for r in result if r["day"]==day]
    print(day,"max |slant product difference|",max(abs(r["dslant_prior_m"]) for r in data),
        "max |slant using estimated ZWD|",max(abs(r["dslant_estimated_m"]) for r in data))
    print([r for r in data if r["el"]==15])
