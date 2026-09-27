#include "gnss_adapter.h"
#include "gnss_signal_policy.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NS_WEEK 604800000000000LL
#define NS_HALF_WEEK 302400000000000LL
#define NS_BDS_GPS 14000000000LL
#define STATE_CODE_LOCK 1u
#define STATE_TOW_DECODED 8u
#define STATE_TOW_KNOWN 16384u

typedef struct { int sys, slot; uint8_t code; double freq; } signal_t;

static int near_frequency(double actual, double nominal)
{
    return isfinite(actual) && fabs(actual - nominal) < 10000.0;
}

static int signal_of(const android_meas_t *m, signal_t *s)
{
    double f = m->carrier_frequency_hz;
    char c = m->code_type[0];
    memset(s, 0, sizeof(*s));
    if (m->constellation_type == 1) {
        s->sys = SYS_GPS;
        if (near_frequency(f, 1575420000.0) && c == 'C')
            s->slot = 0, s->code = CODE_L1C, s->freq = 1575420000.0;
        else if (near_frequency(f, 1176450000.0) && (c == 'I' || c == 'Q' || c == 'X'))
            s->slot = 2, s->code = obs2code(c == 'I' ? "5I" : c == 'Q' ? "5Q" : "5X"), s->freq = 1176450000.0;
    }
    else if (m->constellation_type == 6) {
        s->sys = SYS_GAL;
        if (near_frequency(f, 1575420000.0) && (c == 'C' || c == 'B' || c == 'X'))
            s->slot = 0, s->code = obs2code(c == 'C' ? "1C" : c == 'B' ? "1B" : "1X"), s->freq = 1575420000.0;
        else if (near_frequency(f, 1207140000.0) && (c == 'I' || c == 'Q' || c == 'X'))
            s->slot = 1, s->code = obs2code(c == 'I' ? "7I" : c == 'Q' ? "7Q" : "7X"), s->freq = 1207140000.0;
        else if (near_frequency(f, 1176450000.0) && (c == 'I' || c == 'Q' || c == 'X'))
            s->slot = 2, s->code = obs2code(c == 'I' ? "5I" : c == 'Q' ? "5Q" : "5X"), s->freq = 1176450000.0;
    }
    else if (m->constellation_type == 5) {
        s->sys = SYS_CMP;
        if (near_frequency(f, 1561098000.0) && (c == 'I' || c == 'Q' || c == 'X'))
            s->slot = 0, s->code = obs2code(c == 'I' ? "2I" : c == 'Q' ? "2Q" : "2X"), s->freq = 1561098000.0;
        else if (near_frequency(f, 1575420000.0) && (c == 'D' || c == 'P' || c == 'X'))
            s->slot = 0, s->code = obs2code(c == 'D' ? "1D" : c == 'P' ? "1P" : "1X"), s->freq = 1575420000.0;
        else if (near_frequency(f, 1207140000.0) && (c == 'I' || c == 'Q' || c == 'X'))
            s->slot = 1, s->code = obs2code(c == 'I' ? "7I" : c == 'Q' ? "7Q" : "7X"), s->freq = 1207140000.0;
        else if (near_frequency(f, 1176450000.0) && (c == 'D' || c == 'P' || c == 'Q' || c == 'X'))
            s->slot = 2, s->code = obs2code(c == 'D' ? "5D" : c == 'P' ? "5P" : c == 'Q' ? "5Q" : "5X"), s->freq = 1176450000.0;
    }
    return s->sys && s->code;
}

void gnss_adapter_init(gnss_adapter_t *ctx)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->adr_unc_max = 1.0;
}

void gnss_adapter_set_adr_unc_max(gnss_adapter_t *ctx, double max, int off)
{
    ctx->adr_unc_max = max;
    ctx->adr_unc_off = off;
}

void gnss_adapter_set_rinex_codes(gnss_adapter_t *ctx,
                                  const uint8_t codes[3][NFREQ])
{
    memcpy(ctx->rinex_code, codes, sizeof(ctx->rinex_code));
    ctx->rinex_compatible = 1;
}

void gnss_adapter_set_ppp_safe(gnss_adapter_t *ctx, int enabled)
{
    ctx->ppp_safe = enabled != 0;
}

static int permitted(const gnss_adapter_t *ctx, const signal_t *s)
{
    int sys = s->sys == SYS_GPS ? 0 : s->sys == SYS_GAL ? 1 :
              s->sys == SYS_CMP ? 2 : -1;
    if (ctx->ppp_safe)
        return gnss_ppp_signal_allowed(s->sys, s->slot, s->code);
    return !ctx->rinex_compatible ||
           (sys >= 0 && ctx->rinex_code[sys][s->slot] == s->code);
}

int gnss_adapter_map_signal(const android_meas_t *m, int *sat,
                            uint8_t *code, int *slot)
{
    signal_t s;
    int id;
    if (!signal_of(m, &s) || !(id = satno(s.sys, m->svid))) return 0;
    if (sat) *sat = id;
    if (code) *code = s.code;
    if (slot) *slot = s.slot;
    return 1;
}

int gnss_adapter_clock_time(const android_clock_t *clock, gtime_t *time)
{
    int week;
    int64_t gps_ns, tow_ns;
    double tow;
    if (!clock || !time || !clock->has_full_bias) return 0;
    gps_ns = clock->time_nanos - clock->full_bias_nanos -
             clock->clock_correction_nanos;
    if (gps_ns <= 0) return 0;
    week = (int)(gps_ns / NS_WEEK);
    tow_ns = gps_ns % NS_WEEK;
    tow = (double)tow_ns * 1e-9 - clock->bias_nanos * 1e-9;
    if (tow < 0.0) week--, tow += 604800.0;
    if (tow >= 604800.0) week++, tow -= 604800.0;
    *time = gpst2time(week, tow);
    return 1;
}

static int cmp_sat(const void *a, const void *b)
{
    return (int)((const obsd_t *)a)->sat - (int)((const obsd_t *)b)->sat;
}

int gnss_adapter_convert(gnss_adapter_t *ctx, const android_clock_t *clock,
                         const android_meas_t *meas, int nmeas,
                         obsd_t *obs, int capacity,
                         gnss_adapter_stats_t *stats)
{
    int i, j, n = 0, week, hcdc_changed;
    int best[MAXSAT][NFREQ], best_score[MAXSAT][NFREQ];
    int best_priority[MAXSAT][NFREQ];
    int64_t gps_ns, tow_ns;
    double tow;
    gtime_t time;
    if (stats) memset(stats, 0, sizeof(*stats));
    if (!ctx || !clock || !meas || !obs || capacity <= 0 ||
        !clock->has_full_bias || nmeas <= 0) return 0;
    gps_ns = clock->time_nanos - clock->full_bias_nanos - clock->clock_correction_nanos;
    if (gps_ns <= 0) return 0;
    week = (int)(gps_ns / NS_WEEK);
    tow_ns = gps_ns % NS_WEEK;
    tow = (double)tow_ns * 1e-9 - clock->bias_nanos * 1e-9;
    if (tow < 0.0) week--, tow += 604800.0;
    if (tow >= 604800.0) week++, tow -= 604800.0;
    time = gpst2time(week, tow);
    hcdc_changed = clock->clock_discontinuity ||
        (ctx->have_hcdc &&
         clock->hardware_clock_discontinuity_count != ctx->last_hcdc);
    ctx->last_hcdc = clock->hardware_clock_discontinuity_count;
    ctx->have_hcdc = 1;

    for (i = 0; i < MAXSAT; i++) for (j = 0; j < NFREQ; j++) {
        best[i][j] = -1;
        best_score[i][j] = -1;
        best_priority[i][j] = -1;
    }
    /* Select one signal per RTKLIB frequency slot before touching arc state.
     * This makes B1I/B1C and repeated raw records independent of TXT order. */
    for (i = 0; i < nmeas; i++) {
        signal_t s;
        const android_meas_t *m = meas + i;
        int sat, score, priority;
        if (!signal_of(m, &s) || !permitted(ctx, &s) ||
            !(sat = satno(s.sys, m->svid))) continue;
        priority = ctx->ppp_safe ?
            gnss_ppp_signal_priority(s.sys, s.slot, s.code) : 0;
        score = ((m->state & STATE_CODE_LOCK) &&
                 (m->state & (STATE_TOW_DECODED | STATE_TOW_KNOWN)) &&
                 m->received_sv_time_nanos > 0 &&
                 m->received_sv_time_uncertainty_nanos <= 500.0) ? 100 : 0;
        if ((m->adr_state & ANDROID_ADR_VALID) &&
            !(m->adr_state & (ANDROID_ADR_RESET | ANDROID_ADR_SLIP)) &&
            fabs(m->adr_meters) > 1e-4) score += 40;
        if (isfinite(m->pseudorange_rate_mps)) score += 10;
        if (!ctx->ppp_safe && s.sys == SYS_CMP && s.slot == 0 &&
            s.code == CODE_L2I) score += 5;
        if (ctx->ppp_safe && isfinite(m->cn0_dbhz) && m->cn0_dbhz > 0.0)
            score += (int)fmin(m->cn0_dbhz, 60.0);
        if (priority > best_priority[sat - 1][s.slot] ||
            (priority == best_priority[sat - 1][s.slot] &&
             score > best_score[sat - 1][s.slot])) {
            best_priority[sat - 1][s.slot] = priority;
            best_score[sat - 1][s.slot] = score;
            best[sat - 1][s.slot] = i;
        }
    }

    for (i = 0; i < nmeas; i++) {
        const android_meas_t *m = meas + i;
        signal_t sig;
        obsd_t *o;
        int sat, k, valid_adr, half_unresolved;
        int64_t flight_ns;
        double flight_s, lambda, q;
        if (stats) stats->input_meas++;
        if (!signal_of(m, &sig) || !permitted(ctx, &sig) || m->svid <= 0 ||
            !(sat = satno(sig.sys, m->svid))) {
            if (stats) stats->unsupported++;
            continue;
        }
        if (best[sat - 1][sig.slot] != i) {
            if (stats) stats->duplicate++;
            continue;
        }
        for (j = 0; j < n && obs[j].sat != sat; j++);
        if (j == n) {
            if (n >= capacity) { if (stats) stats->overflow++; continue; }
            memset(obs + n, 0, sizeof(obsd_t));
            obs[n].sat = (uint8_t)sat;
            obs[n].rcv = 1;
            obs[n].time = time;
            n++;
        }
        o = obs + j;
        k = sig.slot;
        o->code[k] = sig.code;
        lambda = CLIGHT / sig.freq;
        if (isfinite(m->cn0_dbhz) && m->cn0_dbhz > 0.0) {
            q = floor(m->cn0_dbhz / SNR_UNIT + 0.5);
            o->SNR[k] = (uint16_t)(q > 65535.0 ? 65535 : q);
        }
        if (m->received_sv_time_nanos > 0 &&
            m->received_sv_time_uncertainty_nanos >= 0.0 &&
            m->received_sv_time_uncertainty_nanos <= 500.0 &&
            (m->state & STATE_CODE_LOCK) &&
            (m->state & (STATE_TOW_DECODED | STATE_TOW_KNOWN))) {
            flight_ns = tow_ns - m->received_sv_time_nanos -
                (sig.sys == SYS_CMP ? NS_BDS_GPS : 0);
            if (flight_ns > NS_HALF_WEEK) flight_ns -= NS_WEEK;
            if (flight_ns < -NS_HALF_WEEK) flight_ns += NS_WEEK;
            flight_s = ((double)flight_ns + m->time_offset_nanos -
                        clock->bias_nanos) * 1e-9;
            if (flight_s >= 0.001 && flight_s <= 0.200) {
                o->P[k] = flight_s * CLIGHT;
                if (stats) stats->code_valid++;
            }
            else if (stats) stats->invalid_time++;
        }
        if (isfinite(m->pseudorange_rate_mps) &&
            isfinite(m->pseudorange_rate_uncertainty_mps) &&
            m->pseudorange_rate_uncertainty_mps >= 0.0 &&
            m->pseudorange_rate_uncertainty_mps <= 10.0) {
            o->D[k] = (float)(-m->pseudorange_rate_mps / lambda);
            if (stats) stats->doppler_valid++;
        }
        valid_adr = (m->adr_state & ANDROID_ADR_VALID) != 0;
        half_unresolved = (m->adr_state & ANDROID_ADR_HALF_REPORTED) &&
                          !(m->adr_state & ANDROID_ADR_HALF_RESOLVED);
        if (valid_adr && isfinite(m->adr_meters) &&
            fabs(m->adr_meters) > 1e-4 &&
            (ctx->adr_unc_off ||
             (isfinite(m->adr_uncertainty_meters) &&
              m->adr_uncertainty_meters >= 0.0 &&
              m->adr_uncertainty_meters <= ctx->adr_unc_max))) {
            o->L[k] = m->adr_meters / lambda;
            if (m->adr_state & (ANDROID_ADR_RESET | ANDROID_ADR_SLIP)) {
                o->LLI[k] |= LLI_SLIP;
                if (stats) stats->lli_android++;
            }
            if (hcdc_changed) {
                o->LLI[k] |= LLI_SLIP;
                if (stats) stats->lli_clock++;
            }
            if (half_unresolved) o->LLI[k] |= LLI_HALFC;
            if (stats) stats->phase_valid++;
        }
        else if (valid_adr && isfinite(m->adr_meters) &&
                 fabs(m->adr_meters) > 1e-4 && stats)
            stats->phase_rejected_unc++;
        if (o->L[k] && (o->LLI[k] & LLI_SLIP) && stats)
            stats->phase_lli_slip++;
    }
    qsort(obs, n, sizeof(obsd_t), cmp_sat);
    return n;
}
