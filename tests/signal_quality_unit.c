/* Compile the actual private helpers, not a separate algorithm copy. */
#undef NDEBUG /* Assertions must also run in the Release test executable. */
#include "../src/ppp.c"
#include <assert.h>
#include <stdio.h>

/* Independent coefficient expansion checks the production recursion,
 * including the shared old Doppler and all past smoothed-code covariances. */
static void test_statistical_model(void)
{
    double v=4.0,rv=.09,c=0.0,past=0.0,w[128]={0},sum[128]={0};
    double a=.9,dt=.7,b=dt*.5,exact,cross,ratecross,raw=4.0;
    gtime_t t=gpst2time(2438,100.0),empty={0};
    prcopt_t opt=prcopt_default;
    int k,j;
    w[0]=1.0;
    for(k=1;k<=50;k++) {
        for(j=0;j<128;j++){sum[j]+=w[j];w[j]*=a;}
        w[k]=1.0-a;
        w[64+k-1]-=a*b;w[64+k]-=a*b;
        smoother_covariance_ppp(a,dt,raw,.09,&v,&rv,&c,&past);
        exact=cross=0.0;
        for(j=0;j<128;j++) {
            double noise=j<64?raw:.09;
            exact+=w[j]*w[j]*noise;cross+=w[j]*sum[j]*noise;
        }
        ratecross=w[64+k]*.09;
        assert(fabs(v-exact)<1E-10);
        assert(fabs(past-cross)<1E-10);
        assert(fabs(c-ratecross)<1E-12);
        assert(v*rv-c*c>=-1E-12); /* PSD of current code/rate covariance */
    }
    assert(fabs(pow(.9,.5)*pow(.9,.5)-.9)<1E-12);
    assert(atmosphere_information_ppp(t,empty,300)==1.0);
    assert(atmosphere_information_ppp(t,t,300)==0.0);
    assert(atmosphere_information_ppp(timeadd(t,-1),t,300)==0.0);
    assert(atmosphere_information_ppp(timeadd(t,900),t,300)==1.0);
    assert(fabs(atmosphere_information_ppp(timeadd(t,1),t,300)-1.0/300)<1E-12);
    assert(fabs(10*atmosphere_information_ppp(timeadd(t,.1),t,300)-1.0/300)<1E-10);
    strcpy(opt.pppopt,"-DOPPSM=nan");assert(doppsm_factor(&opt)==0.0);
    strcpy(opt.pppopt,"-STATMOD=0");assert(statistical_model_ppp(&opt)==0);
    {
        rtk_t *rtk=(rtk_t *)calloc(1,sizeof(rtk_t));
        nav_t *nav=(nav_t *)calloc(1,sizeof(nav_t));
        obsd_t obs={0};double p,expected,lambda=CLIGHT/1575.42E6;
        opt.mode=PMODE_PPP_STATIC;opt.nf=1;
        strcpy(opt.pppopt,"-DOPPSM=0.90 -STATMOD=3");
        rtkinit(rtk,&opt);obs.sat=1;obs.code[0]=CODE_L1C;
        obs.time=t;obs.P[0]=2E7;obs.D[0]=-1000;obs.SNR[0]=(uint16_t)(40/SNR_UNIT);
        assert(uddoppsm_ppp(rtk,&obs,1,nav,NULL)==0);
        obs.time=timeadd(t,.1);obs.P[0]+=2.0;
        expected=(1.0-pow(.9,.1))*obs.P[0]+pow(.9,.1)*(2E7+lambda*1000*.1);
        assert(uddoppsm_ppp(rtk,&obs,1,nav,NULL)==1);
        assert(fabs(rtk->ssat[0].psmP[0]-expected)<1E-8);
        p=rtk->ssat[0].psmP[0];
        assert(uddoppsm_ppp(rtk,&obs,1,nav,NULL)==0);assert(rtk->ssat[0].psmP[0]==p);
        obs.time=timeadd(t,10.0);obs.P[0]+=2;
        assert(uddoppsm_ppp(rtk,&obs,1,nav,NULL)==0);
        assert(rtk->ssat[0].psmP[0]==obs.P[0]&&rtk->ssat[0].psm_rate_cov[0]==0);
        obs.time=timeadd(t,11);obs.LLI[0]=LLI_SLIP;
        uddoppsm_ppp(rtk,&obs,1,nav,NULL);assert(!rtk->ssat[0].psmvalid[0]);
        obs.LLI[0]=0;obs.D[0]=0;
        uddoppsm_ppp(rtk,&obs,1,nav,NULL);assert(!rtk->ssat[0].psmvalid[0]);
        rtkfree(rtk);free(rtk);free(nav);
    }
    printf("PASS: smoother independent covariance expansion (50 steps), PSD, time scaling and prior budgets\n");
}

int main(void)
{
    test_statistical_model();
    prcopt_t opt=prcopt_default;
    ppp_signal_quality_t q={0};
    rtk_t *rtk=(rtk_t *)calloc(1,sizeof(rtk_t));
    ppp_diag_epoch_t evidence={0},final={0};
    obsd_t obs[6]={{0}};
    unsigned char mexc[MAXOBS*2*NFREQ]={0};
    gtime_t start=gpst2time(2438,100.0);
    int i,j;
    opt.mode=PMODE_PPP_STATIC;opt.nf=4;
    strcpy(opt.pppopt,"-PPPQUAR=3,120 -SIGW=G5Q:2:3;E1C:4:1");
    rtkinit(rtk,&opt);
    assert(signal_noise_index_ppp(SYS_GPS,CODE_L5Q,0)==2);
    assert(signal_noise_index_ppp(SYS_GAL,CODE_L1C,3)==0);
    assert(signal_noise_index_ppp(SYS_CMP,CODE_L1P,0)==3);
    assert(signal_noise_index_ppp(SYS_CMP,CODE_L2I,3)==0);
    assert(signal_weight_ppp(SYS_GPS,CODE_L5Q,1,&opt)==2);
    assert(signal_weight_ppp(SYS_GPS,CODE_L5Q,0,&opt)==3);
    assert(signal_weight_ppp(SYS_GAL,CODE_L1C,1,&opt)==4);
    assert(signal_weight_ppp(SYS_CMP,CODE_L5Q,1,&opt)==1);
    assert(varerr(1,SYS_GPS,.8,40,signal_noise_index_ppp(SYS_GPS,CODE_L5Q,0),1,&opt)==
           varerr(1,SYS_GPS,.8,40,signal_noise_index_ppp(SYS_GPS,CODE_L5Q,3),1,&opt));
    strcpy(opt.pppopt,"-SIGW=G5Q:nan:1");assert(signal_weight_ppp(SYS_GPS,CODE_L5Q,1,&opt)==1);
    strcpy(opt.pppopt,"-SIGW=G5Q:0.1:1");assert(signal_weight_ppp(SYS_GPS,CODE_L5Q,1,&opt)==1);
    strcpy(opt.pppopt,"-SIGW=G5Q:2:1garbage");assert(signal_weight_ppp(SYS_GPS,CODE_L5Q,1,&opt)==1);
    for(i=0;i<5;i++) {
        signal_quality_step_ppp(&q,timeadd(start,i),1,0);
        signal_quality_step_ppp(&q,timeadd(start,i),1,0);
    }
    assert(q.weak && q.bad_count==5); /* duplicate timestamp is not new evidence */
    assert(persistent_variance_factor_ppp(&opt,&q,1)==4);
    assert(persistent_variance_factor_ppp(&opt,&q,0)==1);
    strcpy(opt.pppopt,"-SIGQCPH=1");assert(persistent_variance_factor_ppp(&opt,&q,0)==4);
    strcpy(opt.pppopt,"-SIGQC=0");assert(persistent_variance_factor_ppp(&opt,&q,1)==1);
    for(i=0;i<9;i++)signal_quality_step_ppp(&q,timeadd(start,5+i),0,1);
    assert(q.weak);
    signal_quality_step_ppp(&q,timeadd(start,14),0,1);assert(!q.weak);
    memset(&q,0,sizeof(q));
    for(i=0;i<8;i++)signal_quality_step_ppp(&q,timeadd(start,i*.1),1,0);
    assert(!q.weak); /* callback rate cannot replace elapsed time */
    signal_quality_step_ppp(&q,timeadd(start,20),1,0);assert(q.bad_count==1);
    for(i=0;i<6;i++) {
        obs[i].sat=(uint8_t)satno(SYS_GPS,i+1);obs[i].code[0]=CODE_L1C;
        obs[i].L[0]=100.0;
        rtk->ssat[obs[i].sat-1].ppp_track_code[0]=CODE_L1C;
        evidence.obs[i][0].signal=CODE_L1C;
        evidence.obs[i][0].code_ok=evidence.obs[i][0].phase_ok=1;
        evidence.obs[i][0].physical_sigma[0]=evidence.obs[i][0].physical_sigma[1]=1;
        evidence.obs[i][0].code_res=200+(i?0:8);
        evidence.obs[i][0].phase_res=200;
    }
    for(j=0;j<5;j++) {
        for(i=0;i<6;i++)obs[i].time=timeadd(start,j);
        update_signal_quality_ppp(rtk,obs,6,&evidence);
    }
    assert(rtk->ssat[0].ppp_res_quality[0][1].weak);
    assert(!rtk->ssat[0].ppp_res_quality[0][0].weak);
    for(i=1;i<6;i++)assert(!rtk->ssat[i].ppp_res_quality[0][1].weak);
    for(j=0;j<10;j++) {
        evidence.obs[0][0].code_res=200;
        for(i=0;i<6;i++)obs[i].time=timeadd(start,5+j);
        update_signal_quality_ppp(rtk,obs,6,&evidence);
    }
    assert(!rtk->ssat[0].ppp_res_quality[0][1].weak);
    /* Four peers are insufficient: sparse groups abstain. */
    evidence.obs[0][0].code_res=208;
    for(j=0;j<8;j++) {
        for(i=0;i<6;i++)obs[i].time=timeadd(start,15+j);
        update_signal_quality_ppp(rtk,obs,4,&evidence);
    }
    assert(!rtk->ssat[0].ppp_res_quality[0][1].weak);
    /* A switch clears both quality and ambiguity history. */
    rtk->ssat[0].ppp_res_quality[0][1].weak=1;
    obs[0].code[0]=CODE_L1P;
    prepare_signal_quality_ppp(rtk,obs,6);
    assert(!rtk->ssat[0].ppp_res_quality[0][1].weak);
    /* Untested phase must not clear a streak. An actual accepted test can. */
    rtk->ssat[0].ppp_phase_reject_streak[0]=2;
    update_phase_quarantine_ppp(rtk,obs,1,mexc,&final);
    assert(rtk->ssat[0].ppp_phase_reject_streak[0]==2);
    final.obs[0][0].phase_ok=1;
    update_phase_quarantine_ppp(rtk,obs,1,mexc,&final);
    assert(rtk->ssat[0].ppp_phase_reject_streak[0]==0);
    rtkfree(rtk);free(rtk);
    puts("PASS: actual-signal weights, duplicate/rate guards, common-mode, persistence/recovery, sparse peers, signal switch, tested-only quarantine");
    return 0;
}
