#ifndef GNSS_REPLAY_H
#define GNSS_REPLAY_H
#include "rtklib.h"
/* Reusable runner; NULL options retains the standalone replay defaults. */
int gnss_replay_run(int argc, char **argv, const prcopt_t *options);
typedef struct {
    gtime_t time;
    sol_t solution; /* snapshot only; the UI cannot access/mutate rtk_t */
    int nobs;
    double processing_ms;
} gnss_replay_epoch_t;
typedef struct {
    void (*epoch)(void *context, const gnss_replay_epoch_t *result);
    int (*cancelled)(void *context);
    void *context;
} gnss_replay_observer_t;
/* Returns 3 on cooperative user cancellation, after saving/closing products. */
int gnss_replay_run_observed(int argc, char **argv, const prcopt_t *options,
                            const gnss_replay_observer_t *observer);
#endif
