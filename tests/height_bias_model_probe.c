/* Read-only access to internal model quantities for the separate audit exe.
 * Include, do not edit, the canonical PPP implementation. */
#include "../src/ppp.c"
int height_probe_trop_index(const prcopt_t *opt) { return IT(opt); }
int height_probe_clock_index(int sys,const prcopt_t *opt) { return IC(ppp_clk_index(sys),opt); }
int height_probe_ion_index(int sat,const prcopt_t *opt) { return II(sat,opt); }
int height_probe_ion_prior(gtime_t time,const nav_t *nav,const double *pos,
    const double *azel,const prcopt_t *opt,double *ion,double *variance)
{ return ionex_vertical_prior(time,nav,pos,azel,opt,ion,variance); }
int height_probe_atmos(gtime_t time,const double *pos,double out[8])
{
    double ah,aw,zhd,zwd,mfh,mfw;
    double azel[2]={0.0,15.0*D2R};
    if(!vmf3_grid_interp(time,pos,&ah,&aw,&zhd,&zwd))return 0;
    if(!vmf3_trop(time,pos,azel,&mfh,&mfw,&zhd,&zwd))return 0;
    out[0]=zhd;out[1]=zwd;out[2]=ah;out[3]=aw;
    out[4]=mfh;out[5]=mfw;out[6]=(double)vmf3_grid[0].time.time;out[7]=(double)vmf3_grid[1].time.time;
    return 1;
}
