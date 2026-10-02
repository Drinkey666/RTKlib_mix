/* Fault injection only in this test runner; never alters input files. */
#include "../src/smartphone_ppp_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static int calls;
static const char *fixture;
static void test_configure(prcopt_t *opt,solopt_t *sol)
{
    char *p;
    smartphone_ppp_configure(opt,sol);
    if(getenv("SIGNAL_TEST_OPTIONS")) {
        const char *options=getenv("SIGNAL_TEST_OPTIONS");
        if(strlen(options)>=sizeof(opt->pppopt))exit(9);
        strcpy(opt->pppopt,options);
    }
    if(getenv("STARTUP_TEST_OFF")) {
        /* Make room without changing the omitted option's default (=0). */
        if((p=strstr(opt->pppopt," -WGTELCN=0")))memmove(p,p+11,strlen(p+11)+1);
        /* Disable the subsequent persistent-QC feature too when comparing
         * against the pre-startup historical filter-state baseline. PPPDIAG
         * changes trace diagnostics only, not observations or model. */
        if((p=strstr(opt->pppopt," -PPPDIAG=1")))memmove(p,p+11,strlen(p+11)+1);
        if((p=strstr(opt->pppopt," -DOPPWARM=0")))memmove(p,p+12,strlen(p+12)+1);
        if(strlen(opt->pppopt)+strlen(" -STARTQC=0 -SIGQC=0 -STATMOD=0")>=sizeof(opt->pppopt))exit(9);
        strcat(opt->pppopt," -STARTQC=0 -SIGQC=0 -STATMOD=0");
    }
}
static int test_rtkpos(rtk_t *rtk,const obsd_t *obs,int n,const nav_t *nav)
{
    obsd_t input[MAXOBS],saved[MAXOBS];int i,f,result;
    if(n>MAXOBS)exit(10);
    memcpy(input,obs,n*sizeof(*obs));
    if(fixture && (calls==0 || (!strcmp(fixture,"PERSISTENT_PRIMARY")&&calls<10))) {
        for(i=0;i<n;i++) {
            if(input[i].sat==satno(SYS_GPS,1)) {
                if(!strcmp(fixture,"PRIMARY_SPIKE")||!strcmp(fixture,"PERSISTENT_PRIMARY")||!strcmp(fixture,"TWO_SPIKES"))input[i].P[0]+=500.0;
                if(!strcmp(fixture,"SECONDARY_SPIKE"))input[i].P[2]+=500.0;
                if(!strcmp(fixture,"CODE_NAN"))input[i].P[0]=NAN;
                if(!strcmp(fixture,"CODE_NEGATIVE"))input[i].P[0]=-1.0;
                if(!strcmp(fixture,"PHASE_NAN"))input[i].L[0]=NAN;
                if(!strcmp(fixture,"PHASE_HUGE"))input[i].L[0]=1E100;
                if(!strcmp(fixture,"DOPPLER_NAN"))input[i].D[0]=(float)NAN;
                if(!strcmp(fixture,"DOPPLER_ZERO"))input[i].D[0]=0.0;
            }
            if(!strcmp(fixture,"TWO_SPIKES")&&input[i].sat==satno(SYS_GPS,3))input[i].P[0]-=400.0;
            if(!strcmp(fixture,"CLOCK_OFFSET"))for(f=0;f<NFREQ;f++) {
                double hz=sat2freq(input[i].sat,input[i].code[f],nav);
                if(input[i].P[f])input[i].P[f]+=1000.0;
                if(input[i].L[f]&&hz>0.0)input[i].L[f]+=1000.0*hz/CLIGHT;
            }
            if(!strcmp(fixture,"SIGNAL_OFFSET")&&satsys(input[i].sat,NULL)==SYS_GPS&&input[i].P[2])input[i].P[2]+=300.0;
        }
    }
    if(fixture&&!strcmp(fixture,"PHASE_SPIKE")&&calls==5)
        for(i=0;i<n;i++)if(input[i].sat==satno(SYS_GPS,1))input[i].L[0]+=500.0;
    if(fixture&&calls>=60&&calls<100)
        for(i=0;i<n;i++)if(input[i].sat==satno(SYS_GPS,1)) {
            if(!strcmp(fixture,"LATE_CODE_PRIMARY"))input[i].P[0]+=40.0;
            if(!strcmp(fixture,"LATE_CODE_SECONDARY"))input[i].P[2]+=40.0;
            if(!strcmp(fixture,"LATE_PHASE")) {
                double hz=sat2freq(input[i].sat,input[i].code[0],nav);
                if(input[i].L[0]&&hz>0)input[i].L[0]+=(calls%2?1.0:-1.0)*hz/CLIGHT;
            }
        }
    if(fixture&&!strcmp(fixture,"DELAYED_START")&&calls<40)
        for(i=0;i<n;i++)for(f=0;f<NFREQ;f++)input[i].P[f]=input[i].L[f]=input[i].D[f]=0.0;
    memcpy(saved,input,n*sizeof(*obs));
    result=rtkpos(rtk,input,n,nav);
    if(memcmp(saved,input,n*sizeof(*obs))){fprintf(stderr,"FAIL: mutated const input\n");exit(11);}
    for(i=0;i<3;i++)if(!isfinite(rtk->sol.rr[i])) {fprintf(stderr,"FAIL: nonfinite solution\n");exit(12);}
    if(calls==0) {
        ssat_t *s=rtk->ssat+satno(SYS_GPS,1)-1;
        fprintf(stderr,"STARTUP_TEST,fixture=%s,Q=%d,psmvalid_G01_F1=%d,code_bad_G01_F3=%d\n",
            fixture?fixture:"CLEAN",rtk->sol.stat,s->psmvalid[0],s->ppp_code_bad[2]);
    }
    calls++;return result;
}
#define smartphone_ppp_configure test_configure
#define rtkpos test_rtkpos
#define main replay_main
#include "../src/gnss_replay.c"
#undef main
#undef rtkpos
#undef smartphone_ppp_configure
int main(int argc,char **argv)
{
    fixture=getenv("STARTUP_TEST_FIXTURE");
    return replay_main(argc,argv);
}
