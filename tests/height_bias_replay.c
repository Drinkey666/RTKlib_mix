/* Isolated diagnostic runner. Production entry points/defaults are untouched. */
#include "../src/smartphone_ppp_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int height_probe_trop_index(const prcopt_t *opt);
int height_probe_clock_index(int sys, const prcopt_t *opt);
int height_probe_atmos(gtime_t time, const double *pos, double out[8]);
int height_probe_ion_index(int sat,const prcopt_t *opt);
int height_probe_ion_prior(gtime_t time,const nav_t *nav,const double *pos,
    const double *azel,const prcopt_t *opt,double *ion,double *variance);
static FILE *height_diag;
static FILE *ion_diag;

static void audit_configure(prcopt_t *opt, solopt_t *sol)
{
    const char *value;
    char error[256];
    smartphone_ppp_configure(opt, sol);
    if ((value=getenv("HEIGHT_AUDIT_OPTS")) && *value) {
        if (strlen(value)>=sizeof(opt->pppopt)) { fprintf(stderr,"Audit options overflow\n"); exit(2); }
        strcpy(opt->pppopt,value);
    }
    if ((value=getenv("HEIGHT_AUDIT_NF")) && *value) opt->nf=atoi(value);
    if ((value=getenv("HEIGHT_AUDIT_NAVSYS")) && *value) opt->navsys=atoi(value);
    if ((value=getenv("HEIGHT_AUDIT_PCV")) && *value) opt->posopt[0]=atoi(value);
    if ((value=getenv("HEIGHT_AUDIT_IONO")) && *value) opt->ionoopt=atoi(value);
    if (!smartphone_ppp_validate(opt,error,sizeof(error))) {
        fprintf(stderr,"Audit configuration: %s\n",error); exit(2);
    }
}
static double height_cov(const rtk_t *rtk, const double u[3], int j)
{
    int i; double result=0.0;
    for(i=0;i<3;i++) result+=u[i]*rtk->P[i+j*rtk->nx];
    return result;
}
static int audit_rtkpos(rtk_t *rtk,const obsd_t *obs,int nobs,const nav_t *nav)
{
    int ok=rtkpos(rtk,obs,nobs,nav), i,j,w,it;
    double pos[3],u[3],varu=0.0,atmos[8],pzt,cuz;
    if(!height_diag||nobs<=0||!rtk->x||norm(rtk->x,3)<1e6) return ok;
    ecef2pos(rtk->x,pos);
    u[0]=cos(pos[0])*cos(pos[1]);u[1]=cos(pos[0])*sin(pos[1]);u[2]=sin(pos[0]);
    for(i=0;i<3;i++)for(j=0;j<3;j++)varu+=u[i]*u[j]*rtk->P[i+j*rtk->nx];
    it=height_probe_trop_index(&rtk->opt);
    pzt=rtk->P[it+it*rtk->nx];cuz=height_cov(rtk,u,it);
    for(i=0;i<8;i++)atmos[i]=NAN;
    height_probe_atmos(obs[0].time,pos,atmos);
    {
        double tow=time2gpst(obs[0].time,&w);
        fprintf(height_diag,"%d,%.9f,%d,%d,%.9f,%.9f,%.6f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f",
          w,tow,rtk->sol.stat,rtk->sol.ns,pos[0]*R2D,pos[1]*R2D,pos[2],
          sqrt(fmax(varu,0.0)),rtk->x[it],sqrt(fmax(pzt,0.0)),
          varu>0&&pzt>0?cuz/sqrt(varu*pzt):NAN,
          atmos[0],atmos[1],rtk->x[it]-atmos[0],atmos[2],atmos[3],atmos[4],atmos[5]);
        for(i=0;i<3;i++) {
            const int sys[3]={SYS_GPS,SYS_GAL,SYS_CMP};
            j=height_probe_clock_index(sys[i],&rtk->opt);
            fprintf(height_diag,",%.9f,%.9f",rtk->x[j],height_cov(rtk,u,j));
        }
        fprintf(height_diag,"\n");
        if(ion_diag && rtk->sol.stat==SOLQ_PPP && rtk->opt.ionoopt==IONOOPT_EST) {
            for(i=0;i<nobs;i++) {
                int sat=obs[i].sat,used=0,f,ii;
                double prior,pvar,ivar;char sid[16];
                if(sat<1||sat>MAXSAT)continue;
                for(f=0;f<rtk->opt.nf;f++)used|=rtk->ssat[sat-1].vsat[f];
                if(!used)continue;
                ii=height_probe_ion_index(sat,&rtk->opt);
                if(!height_probe_ion_prior(obs[i].time,nav,pos,rtk->ssat[sat-1].azel,
                    &rtk->opt,&prior,&pvar))continue;
                ivar=rtk->P[ii+ii*rtk->nx];satno2id(sat,sid);
                fprintf(ion_diag,"%d,%.9f,%s,%.3f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n",
                    w,tow,sid,rtk->ssat[sat-1].azel[1]*R2D,rtk->x[ii],prior,
                    rtk->x[ii]-prior,sqrt(fmax(pvar,0.0)),sqrt(fmax(ivar,0.0)),
                    varu>0&&ivar>0?height_cov(rtk,u,ii)/sqrt(varu*ivar):NAN);
            }
        }
    }
    return ok;
}
#define smartphone_ppp_configure audit_configure
#define rtkpos audit_rtkpos
#define main replay_diagnostic_main
#include "../src/gnss_replay.c"
#undef main
#undef rtkpos
#undef smartphone_ppp_configure

int main(int argc,char **argv)
{
    const char *path=getenv("HEIGHT_AUDIT_DIAG");
    int result;
    if(path&&*path) {
        height_diag=fopen(path,"wb");
        if(!height_diag){fprintf(stderr,"Cannot create audit output\n");return 2;}
        fprintf(height_diag,"week,tow,Q,ns,lat_deg,lon_deg,h_m,sigma_U_m,ZTD_m,sigma_ZTD_m,corr_U_ZTD,ZHD_m,ZWDprior_m,ZWDest_m,ah,aw,mfh15,mfw15,GPS_clock_m,cov_U_GPSclock,GAL_clock_m,cov_U_GALclock,BDS_clock_m,cov_U_BDSclock\n");
    }
    path=getenv("HEIGHT_AUDIT_IONDIAG");
    if(path&&*path) {
        ion_diag=fopen(path,"wb");
        if(!ion_diag){if(height_diag)fclose(height_diag);return 2;}
        fprintf(ion_diag,"week,tow,sat,el_deg,Iest_m,Iprior_m,delta_m,sigma_prior_m,sigma_state_m,corr_U_ion\n");
    }
    result=replay_diagnostic_main(argc,argv);
    if(height_diag)fclose(height_diag);
    if(ion_diag)fclose(ion_diag);
    return result;
}
