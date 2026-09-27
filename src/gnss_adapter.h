#ifndef GNSS_ADAPTER_H
#define GNSS_ADAPTER_H

#include "gnss_raw.h"
#include "rtklib.h"

typedef struct {
    int input_meas, unsupported, code_valid, phase_valid, doppler_valid;
    int invalid_time, duplicate, overflow;
    int phase_rejected_unc, phase_lli_slip;
    int lli_android, lli_clock;
} gnss_adapter_stats_t;

typedef struct {
    int last_hcdc, have_hcdc;
    double adr_unc_max;
    int adr_unc_off;
    int rinex_compatible;
    int ppp_safe;
    uint8_t rinex_code[3][NFREQ]; /* GPS, GAL, BDS; learned from rinex.c output */
} gnss_adapter_t;

void gnss_adapter_init(gnss_adapter_t *ctx);
void gnss_adapter_set_adr_unc_max(gnss_adapter_t *ctx, double max, int off);
void gnss_adapter_set_rinex_codes(gnss_adapter_t *ctx,
                                  const uint8_t codes[3][NFREQ]);
void gnss_adapter_set_ppp_safe(gnss_adapter_t *ctx, int enabled);
int gnss_adapter_map_signal(const android_meas_t *m, int *sat,
                            uint8_t *code, int *slot);
int gnss_adapter_clock_time(const android_clock_t *clock, gtime_t *time);
int gnss_adapter_convert(gnss_adapter_t *ctx, const android_clock_t *clock,
                         const android_meas_t *meas, int nmeas,
                         obsd_t *obs, int capacity,
                         gnss_adapter_stats_t *stats);

#endif
