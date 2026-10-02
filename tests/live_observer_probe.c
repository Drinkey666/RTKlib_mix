/* Exercise the same read-only callback/cancellation API used by the window. */
#include "../src/gnss_replay.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct { int count,stop_after; } test_context_t;
static void epoch(void *context,const gnss_replay_epoch_t *r)
{
    test_context_t *c=(test_context_t *)context;
    double pos[3]={0},tow;int week;
    c->count++;
    tow=time2gpst(r->time,&week);
    if(r->solution.stat) ecef2pos(r->solution.rr,pos);
    printf("UI_EPOCH,%d,%.9f,%d,%d,%d,%.9f,%.9f,%.9f,%.12f,%.12f,%.6f\n",
        week,tow,r->solution.stat,r->solution.ns,r->nobs,
        r->solution.rr[0],r->solution.rr[1],r->solution.rr[2],
        pos[0]*R2D,pos[1]*R2D,pos[2]);
}
static int cancelled(void *context)
{
    test_context_t *c=(test_context_t *)context;
    return c->stop_after>0 && c->count>=c->stop_after;
}
int main(int argc,char **argv)
{
    test_context_t context={0};
    const char *stop=getenv("GNSS_OBSERVER_STOP_AFTER");
    gnss_replay_observer_t observer={epoch,cancelled,&context};
    int result;
    if(stop) context.stop_after=atoi(stop);
    result=gnss_replay_run_observed(argc,argv,NULL,&observer);
    printf("UI_CALLBACK_COUNT,%d,exit=%d\n",context.count,result);
    return result;
}
