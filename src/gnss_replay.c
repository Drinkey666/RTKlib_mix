/* PC Android-Raw replay entry. The PPP engine and product readers are reused,
 * but observation records never pass through the RINEX OBS reader. */
#include "gnss_adapter.h"
#include "gnss_signal_policy.h"
#include <windows.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RAW_PER_EPOCH 256
#define CSV_COLS 80
#define LINE_SIZE 8192
#define NS_CLOCK_JUMP 50000LL

typedef struct {
    const char *txt, *nav, *sp3, *clk, *bia, *ionex;
    const char *vmf0, *vmf6, *orog, *atx, *rinex, *dump, *trace;
    const char *result_file, *state_dump, *use_dump;
    const char *exclude_raw_signal; /* diagnostic ablation, e.g. C:5Q */
    int realtime, max_epochs, rinex_source, adr_unc_off, trace_level;
    int signal_mode; /* 0=full, 1=rinex-compatible, 2=common-only, 3=ppp-safe */
    double adr_unc_max;
} replay_args_t;

typedef struct {
    int have_previous;
    int64_t time_nanos, gps_nanos, correction_nanos;
    int hcdc;
} clock_history_t;

typedef struct {
    double *ms;
    int n, cap, deadline_miss, processing_gt_1000ms;
    int status_count[8];
    int raw_meas, code_meas, phase_meas, doppler_meas, unsupported, duplicates;
    int phase_rejected_unc, phase_lli_slip;
    int lli_android, lli_clock;
    int adr_raw, adr_valid, adr_reset, adr_slip, adr_half, adr_unc_invalid;
    int adr_unc_bins[5];
    double *adr_unc;
    int adr_unc_n, adr_unc_cap;
    double sum, max;
} perf_t;

typedef struct {
    int epochs, matched_epochs, missing_epochs, matched_signals, missing_signals;
    int code_mismatch, unrepresented_signals, different_code, time_mismatch;
    double max_p, max_l, max_d, max_snr, max_time;
    double sum_p, sum_l, sum_d, sum_snr;
    int np, nl, nd, ns;
    int rinex_phase, phase_matched, phase_missing[9];
    int strict_p_n, strict_d_n, lli_match, lli_mismatch;
    double phase_sum, phase_sq, phase_max;
    double prr_sum, prr_sq, prr_max;
    double p_sum, p_sq, p_max;
} compare_t;

enum {
    MISS_NO_RAW_SIGNAL, MISS_UNSUPPORTED_SIGNAL, MISS_ADR_NOT_VALID,
    MISS_ADR_ZERO, MISS_ADR_NAN, MISS_ADR_UNCERTAINTY,
    MISS_SIGNAL_SELECTION_LOST, MISS_CODE_MAPPING_DIFFERENT, MISS_OTHER
};

static double monotonic_ms(void)
{
    LARGE_INTEGER counter, freq;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart * 1000.0 / (double)freq.QuadPart;
}

static void wait_until(double deadline_ms)
{
    for (;;) {
        double remaining = deadline_ms - monotonic_ms();
        if (remaining <= 0.0) return;
        if (remaining > 3.0) Sleep((DWORD)(remaining - 2.0));
        else SwitchToThread();
    }
}

static void usage(const char *exe)
{
    fprintf(stderr,
        "Usage: %s --source txt|txt-fast|txt-realtime|rinex --txt RawData.txt [options]\n"
        "  --source txt-fast/txt-realtime selects fast/1x monotonic-clock TXT replay\n"
        "  --mode fast|realtime remains available with --source txt\n"
        "  --nav FILE --sp3 FILE --clk FILE --bia FILE --ionex FILE\n"
        "  --vmf0 FILE --vmf6 FILE --orog FILE --atx FILE\n"
        "  --rinex-obs FILE   TXT comparison, or PPP input with --source rinex\n"
        "  --signal-mode ppp-safe|full|rinex-compatible|common-only (default ppp-safe)\n"
        "  --obs-dump FILE    output accepted observation CSV\n"
        "  --result-file FILE --state-dump FILE --use-dump FILE\n"
        "  --exclude-raw-signal C:5Q (diagnostic ablation only)\n"
        "  --adr-unc-max 0.5|1.0|2.0|5.0|off (default 1.0)\n"
        "  --trace FILE --trace-level 2|3 --max-epochs N\n", exe);
}

static int parse_args(int argc, char **argv, replay_args_t *a)
{
    int i;
    memset(a, 0, sizeof(*a));
    a->signal_mode = 3;
    a->adr_unc_max = 1.0;
    a->trace_level = 2;
    a->nav = "E:\\RTKLIB_Data\\NAV\\BRDC00WRD_S_20262650000_01D_MN.rnx";
    a->sp3 = "E:\\RTKLIB_Data\\SP3\\WUM0MGXNRT_20262641000_02D_05M_ORB.SP3";
    a->clk = "E:\\RTKLIB_Data\\CLK\\WUM0MGXNRT_20262641000_02D_05M_CLK.CLK";
    a->bia = "E:\\RTKLIB_Data\\BIA\\WUM0MGXRTS_20262650000_01D_05M_OSB.BIA";
    a->ionex = "E:\\RTKLIB_Data\\IONEX\\COD0OPSP0D_20262650000_01D_01H_GIM.INX";
    /* This capture is after 06:00 UTC (06:00:18 GPST); bracket it with
     * the 06 and 12 UTC grids, not the 00 and 06 grids. */
    a->vmf0 = "E:\\RTKLIB_Data\\tro\\VMF3_20260922.H06";
    a->vmf6 = "E:\\RTKLIB_Data\\tro\\VMF3_20260922.H12";
    a->orog = "E:\\RTKLIB_Data\\tro\\orography_ell_5x5";
    a->atx = "E:\\RTKLIB_Data\\Tables\\igs20.atx";
    a->trace = "gnss_replay.trace";
    for (i = 1; i < argc; i++) {
        const char *v;
        if (i + 1 >= argc) return 0;
        v = argv[++i];
        if (!strcmp(argv[i - 1], "--txt")) a->txt = v;
        else if (!strcmp(argv[i - 1], "--source")) {
            if (!strcmp(v, "txt")) a->rinex_source = 0;
            else if (!strcmp(v, "txt-fast"))
                a->rinex_source = 0, a->realtime = 0;
            else if (!strcmp(v, "txt-realtime"))
                a->rinex_source = 0, a->realtime = 1;
            else if (!strcmp(v, "rinex")) a->rinex_source = 1;
            else return 0;
        }
        else if (!strcmp(argv[i - 1], "--mode")) {
            if (!strcmp(v, "fast")) a->realtime = 0;
            else if (!strcmp(v, "realtime")) a->realtime = 1;
            else return 0;
        }
        else if (!strcmp(argv[i - 1], "--nav")) a->nav = v;
        else if (!strcmp(argv[i - 1], "--sp3")) a->sp3 = v;
        else if (!strcmp(argv[i - 1], "--clk")) a->clk = v;
        else if (!strcmp(argv[i - 1], "--bia")) a->bia = v;
        else if (!strcmp(argv[i - 1], "--ionex")) a->ionex = v;
        else if (!strcmp(argv[i - 1], "--vmf0")) a->vmf0 = v;
        else if (!strcmp(argv[i - 1], "--vmf6")) a->vmf6 = v;
        else if (!strcmp(argv[i - 1], "--orog")) a->orog = v;
        else if (!strcmp(argv[i - 1], "--atx")) a->atx = v;
        else if (!strcmp(argv[i - 1], "--rinex-obs")) a->rinex = v;
        else if (!strcmp(argv[i - 1], "--signal-mode")) {
            if (!strcmp(v, "full")) a->signal_mode = 0;
            else if (!strcmp(v, "rinex-compatible")) a->signal_mode = 1;
            else if (!strcmp(v, "common-only")) a->signal_mode = 2;
            else if (!strcmp(v, "ppp-safe")) a->signal_mode = 3;
            else return 0;
        }
        else if (!strcmp(argv[i - 1], "--obs-dump")) a->dump = v;
        else if (!strcmp(argv[i - 1], "--result-file")) a->result_file = v;
        else if (!strcmp(argv[i - 1], "--state-dump")) a->state_dump = v;
        else if (!strcmp(argv[i - 1], "--use-dump")) a->use_dump = v;
        else if (!strcmp(argv[i - 1], "--exclude-raw-signal")) {
            if (strlen(v) != 4 || v[1] != ':' ||
                (v[0] != 'G' && v[0] != 'E' && v[0] != 'C') ||
                !obs2code(v + 2)) return 0;
            a->exclude_raw_signal = v;
        }
        else if (!strcmp(argv[i - 1], "--adr-unc-max")) {
            char *end;
            if (!strcmp(v, "off")) a->adr_unc_off = 1;
            else {
                a->adr_unc_max = strtod(v, &end);
                if (!*v || *end || !isfinite(a->adr_unc_max) ||
                    a->adr_unc_max < 0.0) return 0;
                a->adr_unc_off = 0;
            }
        }
        else if (!strcmp(argv[i - 1], "--trace")) a->trace = v;
        else if (!strcmp(argv[i - 1], "--trace-level")) {
            a->trace_level = atoi(v);
            if (a->trace_level < 2 || a->trace_level > 3) return 0;
        }
        else if (!strcmp(argv[i - 1], "--max-epochs")) a->max_epochs = atoi(v);
        else return 0;
    }
    return (a->rinex_source ? a->rinex != NULL : a->txt != NULL) &&
           ((a->signal_mode != 1 && a->signal_mode != 2) || a->rinex != NULL);
}

static int exists(const char *path)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "Missing input: %s\n", path); return 0; }
    fclose(fp);
    return 1;
}

static int split_csv(char *line, char **fields, int cap)
{
    int n = 0;
    char *p = line;
    while (n < cap) {
        fields[n++] = p;
        p = strchr(p, ',');
        if (!p) break;
        *p++ = 0;
    }
    return n;
}

static int64_t as_i64(const char *s) { return s && *s ? _strtoi64(s, NULL, 10) : 0; }
static double as_double(const char *s) { return s && *s ? strtod(s, NULL) : 0.0; }
static int as_int(const char *s) { return s && *s ? (int)strtol(s, NULL, 10) : 0; }

static int parse_raw(char *line, android_clock_t *clock, android_meas_t *m,
                     const android_clock_t *last)
{
    char *f[CSV_COLS];
    int n = split_csv(line, f, CSV_COLS);
    if (n < 37 || strcmp(f[0], "Raw")) return 0;
    memset(clock, 0, sizeof(*clock));
    memset(m, 0, sizeof(*m));
    clock->time_nanos = as_i64(f[2]);
    clock->full_bias_nanos = *f[5] ? as_i64(f[5]) : last->full_bias_nanos;
    clock->has_full_bias = *f[5] != 0 || last->has_full_bias;
    clock->bias_nanos = *f[6] ? as_double(f[6]) : last->bias_nanos;
    clock->hardware_clock_discontinuity_count = *f[10] ? as_int(f[10]) :
        last->hardware_clock_discontinuity_count;
    m->svid = as_int(f[11]);
    m->time_offset_nanos = as_double(f[12]);
    m->state = (uint32_t)as_int(f[13]);
    m->received_sv_time_nanos = as_i64(f[14]);
    m->received_sv_time_uncertainty_nanos = as_double(f[15]);
    m->cn0_dbhz = as_double(f[16]);
    m->pseudorange_rate_mps = as_double(f[17]);
    m->pseudorange_rate_uncertainty_mps = as_double(f[18]);
    m->adr_state = (uint32_t)as_int(f[19]);
    m->adr_meters = as_double(f[20]);
    m->adr_uncertainty_meters = as_double(f[21]);
    m->carrier_frequency_hz = as_double(f[22]);
    m->constellation_type = as_int(f[28]);
    if (*f[35]) {
        strncpy(m->code_type, f[35], sizeof(m->code_type) - 1);
        m->code_type[0] = (char)toupper((unsigned char)m->code_type[0]);
    }
    return 1;
}

static void normalize_clock(clock_history_t *h, android_clock_t *c)
{
    int64_t raw_gps, elapsed, step;
    if (!c->has_full_bias) return;
    raw_gps = c->time_nanos - c->full_bias_nanos;
    if (h->have_previous) {
        elapsed = c->time_nanos - h->time_nanos;
        step = (raw_gps - h->gps_nanos) - elapsed;
        if (c->hardware_clock_discontinuity_count == h->hcdc &&
            elapsed > 0 && elapsed <= 5000000000LL &&
            (step >= NS_CLOCK_JUMP || step <= -NS_CLOCK_JUMP)) {
            h->correction_nanos += step;
            c->clock_discontinuity = 1;
            fprintf(stderr, "CLOCK_STEP,raw_ns=%lld,correction_ns=%lld\n",
                    (long long)step, (long long)h->correction_nanos);
        }
    }
    h->time_nanos = c->time_nanos;
    h->gps_nanos = raw_gps;
    h->hcdc = c->hardware_clock_discontinuity_count;
    h->have_previous = 1;
    c->clock_correction_nanos = h->correction_nanos;
}

static void configure_ppp(prcopt_t *opt)
{
    *opt = prcopt_default;
    opt->mode = PMODE_PPP_STATIC;
    opt->dynamics = 0;
    opt->navsys = SYS_GPS | SYS_GAL | SYS_CMP;
    opt->nf = 3;
    opt->sateph = EPHOPT_PREC;
    opt->ionoopt = IONOOPT_EST;
    opt->tropopt = TROPOPT_EST;
    opt->err[1] = 0.008; opt->err[2] = 0.012;
    opt->eratio[0] = 100.0; opt->eratio[1] = 70.0; opt->eratio[2] = 60.0;
    opt->modear = ARMODE_OFF;
    opt->thresslip = 0.20;
    opt->elmin = 15.0 * D2R;
    opt->prn[0] = 1e-4; opt->prn[1] = 1e-3; opt->prn[2] = 1e-4;
    opt->maxout = 30;
    opt->tidecorr = 1;
    opt->posopt[0] = 1; opt->posopt[1] = 0;
    opt->posopt[2] = 1; opt->posopt[3] = 1;
    opt->snrmask.ena[0] = 0;
    strcpy(opt->pppopt,
        "-GAP_RESION=120 -IONCONS=1.5 -IONCONSINT=1 -VMF3SIG=0.15 -VMF3ZWDSIG=0.30 -VMF3ZWDINT=300 -DOPPSM=0.90 -DOPPWARM=10 -PREPROC=1 -DOPPSLIP=0.50 -CODEJUMP=30 -MWTHRES=5 -BDSCODEVAR=1 -BDSCODEWARM=120 -WGTELCN=0 -PPPDIAG=1 -PPPQUAR=3,120");
}

static void free_products(nav_t *nav)
{
    int i;
    free(nav->eph); free(nav->geph); free(nav->seph);
    free(nav->peph); free(nav->pclk); free(nav->osbs);
    for (i = 0; i < nav->nt; i++) {
        free(nav->tec[i].data); free(nav->tec[i].rms);
    }
    free(nav->tec); free(nav->erp.data);
}

/* Inspect actual TXT satellite/code availability after loading the BIA.
 * This is diagnostic only: OSB absence does not itself gate a signal because
 * the existing PPP engine may use its legacy DCB fallback. */
static int scan_signal_products(const char *path, const nav_t *nav)
{
    static const struct { int sys, slot; const char *name; } target[] = {
        {SYS_GPS,0,"1C"}, {SYS_GPS,2,"5Q"},
        {SYS_GAL,0,"1C"}, {SYS_GAL,1,"7Q"}, {SYS_GAL,2,"5Q"},
        {SYS_CMP,0,"2I"}, {SYS_CMP,0,"1P"}, {SYS_CMP,1,"7I"},
        {SYS_CMP,2,"5P"}, {SYS_CMP,2,"5Q"}
    };
    uint8_t seen[MAXSAT][MAXCODE + 1] = {{0}};
    android_clock_t clock, last = {0};
    android_meas_t meas;
    char line[LINE_SIZE];
    FILE *fp = fopen(path, "rb");
    int i, sat, slot;
    uint8_t code;
    if (!fp) return 0;
    while (fgets(line, sizeof(line), fp)) {
        if (!parse_raw(line, &clock, &meas, &last)) continue;
        last = clock;
        if (gnss_adapter_map_signal(&meas, &sat, &code, &slot) &&
            sat >= 1 && sat <= MAXSAT && code <= MAXCODE)
            seen[sat - 1][code] = 1;
    }
    fclose(fp);
    for (i = 0; i < (int)(sizeof(target) / sizeof(target[0])); i++) {
        int s, observed = 0, valid = 0, missing = 0;
        int wanted = obs2code(target[i].name);
        const char *sys = target[i].sys == SYS_GPS ? "GPS" :
                          target[i].sys == SYS_GAL ? "GAL" : "BDS";
        for (s = 1; s <= MAXSAT; s++) {
            if (satsys(s, NULL) != target[i].sys || !seen[s - 1][wanted])
                continue;
            observed++;
            if (nav->osb_valid[s - 1][wanted]) valid++;
            else missing++;
        }
        fprintf(stderr,
            "SIGNAL_PRODUCT,%s,%s,slot=%d,observed_satellites=%d,osb_valid_satellites=%d,osb_missing_satellites=%d,ppp_safe=%d,reason=%s\n",
            sys, target[i].name, target[i].slot + 1, observed, valid,
            missing, gnss_ppp_signal_allowed(target[i].sys, target[i].slot,
                                             (uint8_t)wanted),
            gnss_ppp_signal_reason(target[i].sys, target[i].slot,
                                   (uint8_t)wanted));
    }
    return 1;
}

static int load_products(const replay_args_t *a, nav_t *nav, prcopt_t *opt)
{
    obs_t unused_obs = {0};
    sta_t sta = {0};
    pcvs_t pcvs = {0};
    int i;
    if (!exists(a->nav) || !exists(a->sp3) || !exists(a->clk) ||
        !exists(a->bia) || !exists(a->ionex) || !exists(a->vmf0) ||
        !exists(a->vmf6) || !exists(a->orog) || !exists(a->atx)) return 0;
    if (!readrnx(a->nav, 1, "", &unused_obs, nav, &sta)) return 0;
    freeobs(&unused_obs);
    uniqnav(nav);
    readsp3(a->sp3, nav, 0);
    if (!readrnxc(a->clk, nav)) return 0;
    if (!readdcb(a->bia, nav, &sta)) return 0;
    readtec(a->ionex, nav, 1);
    if (!readpcv(a->atx, &pcvs)) return 0;
    for (i = 0; i < MAXSAT; i++) {
        pcv_t *pcv;
        if (!(satsys(i + 1, NULL) & opt->navsys)) continue;
        pcv = searchpcv(i + 1, "", nav->peph && nav->ne ? nav->peph[0].time : gpst2time(0, 0), &pcvs);
        if (pcv) nav->pcvs[i] = *pcv;
    }
    free(pcvs.pcv);
    if (!pppvmf3load(a->vmf0, a->vmf6, a->orog)) return 0;
    fprintf(stderr, "PRODUCTS,nav=%d,sp3=%d,clk=%d,bia=%d,ionex=%d\n",
            nav->n, nav->ne, nav->nc, nav->nosb, nav->nt);
    return nav->n > 0 && nav->ne > 0 && nav->nc > 0 && nav->nosb > 0 && nav->nt > 0;
}

static int rinex_sys_index(int sys)
{
    return sys == SYS_GPS ? 0 : sys == SYS_GAL ? 1 :
           sys == SYS_CMP ? 2 : -1;
}

/* Diagnostic mirror of ppp.c::ppp_ifb_index(), not an estimator setting.
 * Keep this label in sync with that function when the PPP model changes. */
static const char *receiver_code_bias_slot(int sys, int f, uint8_t code,
                                           int nf)
{
    const char *sig = code2obs(code);
    if (f == 0) return "reference_code";
    if (nf < 3 || !sig) return "unmodelled";
    if (f == 2 && !strcmp(sig, "5Q")) {
        if (sys == SYS_GPS || sys == SYS_QZS) return "GPS_C5Q";
        if (sys == SYS_GAL) return "GAL_C5Q";
    }
    if (sys == SYS_CMP && f == 1 && !strcmp(sig, "7I"))
        return "BDS_C7I";
    if (sys == SYS_CMP && f == 2 && !strcmp(sig, "5P"))
        return "BDS_C5P";
    return "unmodelled";
}

/* rinex.c has already applied its header-based getcodepri()/set_index()
 * selection. Learn the actually selected code in each system/slot instead of
 * imposing a BDS 2I/1P preference in the adapter. */
static int learn_rinex_codes(const obs_t *rnx, uint8_t codes[3][NFREQ])
{
    int i, f;
    memset(codes, 0, sizeof(uint8_t) * 3 * NFREQ);
    for (i = 0; i < rnx->n; i++) {
        int sys = satsys(rnx->data[i].sat, NULL);
        int si = rinex_sys_index(sys);
        if (si < 0) continue;
        for (f = 0; f < NFREQ; f++) {
            uint8_t code = rnx->data[i].code[f];
            if (!code) continue;
            if (codes[si][f] && codes[si][f] != code) {
                fprintf(stderr, "RINEX_SELECTION_CONFLICT,sys=%d,slot=%d,%s,%s\n",
                        sys, f + 1, code2obs(codes[si][f]), code2obs(code));
                return 0;
            }
            codes[si][f] = code;
        }
    }
    for (i = 0; i < 3; i++) for (f = 0; f < NFREQ; f++) {
        int sys = i == 0 ? SYS_GPS : i == 1 ? SYS_GAL : SYS_CMP;
        if (codes[i][f]) fprintf(stderr,
            "RINEX_SELECTION,sys=%d,slot=%d,code=%s,reader_priority=%d\n",
            sys, f + 1, code2obs(codes[i][f]),
            getcodepri(sys, codes[i][f], ""));
    }
    return 1;
}

static int filter_common_only(obsd_t *obs, int n, const obs_t *rnx,
                              int cursor)
{
    int i, j, f, out = 0, start = cursor;
    if (!n || !rnx || !rnx->data) return 0;
    while (start < rnx->n &&
           timediff(rnx->data[start].time, obs[0].time) < -0.01) start++;
    if (start == rnx->n ||
        fabs(timediff(rnx->data[start].time, obs[0].time)) > 0.01) return 0;
    for (i = 0; i < n; i++) {
        int any = 0;
        for (f = 0; f < NFREQ; f++) {
            int found = 0, k;
            if (!obs[i].code[f]) continue;
            for (j = start; j < rnx->n &&
                 fabs(timediff(rnx->data[j].time, obs[0].time)) <= 0.01; j++) {
                if (rnx->data[j].sat != obs[i].sat) continue;
                for (k = 0; k < NFREQ; k++) {
                    if (rnx->data[j].code[k] == obs[i].code[f]) found = 1;
                }
            }
            if (found) any = 1;
            else {
                obs[i].code[f] = CODE_NONE;
                obs[i].P[f] = obs[i].L[f] = obs[i].D[f] = 0.0;
                obs[i].SNR[f] = obs[i].LLI[f] = 0;
            }
        }
        if (any) obs[out++] = obs[i];
    }
    return out;
}

static int select_common_raw(const android_clock_t *clock,
                             const android_meas_t *raw, int nraw,
                             const obs_t *rnx, int cursor,
                             android_meas_t *selected)
{
    gtime_t time;
    int i, j, f, start = cursor, count = 0;
    if (!rnx || !gnss_adapter_clock_time(clock, &time)) return 0;
    while (start < rnx->n && timediff(rnx->data[start].time, time) < -0.01)
        start++;
    if (start == rnx->n || fabs(timediff(rnx->data[start].time, time)) > 0.01)
        return 0;
    for (i = 0; i < nraw; i++) {
        int sat, found = 0;
        uint8_t code;
        if (!gnss_adapter_map_signal(raw + i, &sat, &code, NULL)) continue;
        for (j = start; j < rnx->n &&
             fabs(timediff(rnx->data[j].time, time)) <= 0.01; j++) {
            if (rnx->data[j].sat != sat) continue;
            for (f = 0; f < NFREQ; f++)
                if (rnx->data[j].code[f] == code) found = 1;
        }
        if (found && count < RAW_PER_EPOCH) selected[count++] = raw[i];
    }
    return count;
}

static int cmp_double(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

static void add_perf(perf_t *p, double ms)
{
    if (p->n == p->cap) {
        int cap = p->cap ? p->cap * 2 : 1024;
        double *next = (double *)realloc(p->ms, sizeof(double) * cap);
        if (!next) return;
        p->ms = next; p->cap = cap;
    }
    p->ms[p->n++] = ms;
    p->sum += ms;
    if (ms > p->max) p->max = ms;
    if (ms > 1000.0) p->processing_gt_1000ms++;
}

static void add_adr_stats(perf_t *p, const android_meas_t *raw, int nraw)
{
    int i, bin;
    for (i = 0; i < nraw; i++) {
        double u = raw[i].adr_uncertainty_meters;
        uint32_t state = raw[i].adr_state;
        p->adr_raw++;
        if (state & ANDROID_ADR_RESET) p->adr_reset++;
        if (state & ANDROID_ADR_SLIP) p->adr_slip++;
        if (state & ANDROID_ADR_HALF_REPORTED) p->adr_half++;
        if (!(state & ANDROID_ADR_VALID)) continue;
        p->adr_valid++;
        if (!isfinite(u) || u < 0.0) { p->adr_unc_invalid++; continue; }
        bin = u <= 0.5 ? 0 : u <= 1.0 ? 1 : u <= 2.0 ? 2 :
              u <= 5.0 ? 3 : 4;
        p->adr_unc_bins[bin]++;
        if (p->adr_unc_n == p->adr_unc_cap) {
            int cap = p->adr_unc_cap ? p->adr_unc_cap * 2 : 1024;
            double *next = (double *)realloc(p->adr_unc, sizeof(double) * cap);
            if (!next) continue;
            p->adr_unc = next;
            p->adr_unc_cap = cap;
        }
        p->adr_unc[p->adr_unc_n++] = u;
    }
}

static void print_adr_stats(perf_t *p)
{
    int n = p->adr_unc_n;
    fprintf(stderr,
        "ADR_FLAGS,raw=%d,valid=%d,reset=%d,slip=%d,half_reported=%d,unc_invalid=%d\n",
        p->adr_raw, p->adr_valid, p->adr_reset, p->adr_slip, p->adr_half,
        p->adr_unc_invalid);
    fprintf(stderr,
        "ADR_UNC_BINS,le_0.5=%d,le_1=%d,le_2=%d,le_5=%d,gt_5=%d\n",
        p->adr_unc_bins[0],
        p->adr_unc_bins[0] + p->adr_unc_bins[1],
        p->adr_unc_bins[0] + p->adr_unc_bins[1] + p->adr_unc_bins[2],
        p->adr_unc_bins[0] + p->adr_unc_bins[1] +
            p->adr_unc_bins[2] + p->adr_unc_bins[3], p->adr_unc_bins[4]);
    if (!n) return;
    qsort(p->adr_unc, n, sizeof(double), cmp_double);
    fprintf(stderr,
        "ADR_UNC_DISTRIBUTION,n=%d,min=%.6f,median=%.6f,p90=%.6f,p95=%.6f,p99=%.6f,max=%.6f\n",
        n, p->adr_unc[0],
        n & 1 ? p->adr_unc[n / 2] :
            (p->adr_unc[n / 2 - 1] + p->adr_unc[n / 2]) / 2.0,
        p->adr_unc[(int)ceil(0.90 * n) - 1],
        p->adr_unc[(int)ceil(0.95 * n) - 1],
        p->adr_unc[(int)ceil(0.99 * n) - 1], p->adr_unc[n - 1]);
}

static int count_trace_tag(const char *path, const char *tag)
{
    FILE *fp;
    char line[2048];
    int count = 0;
    if (!path || !(fp = fopen(path, "rb"))) return -1;
    while (fgets(line, sizeof(line), fp)) if (strstr(line, tag)) count++;
    fclose(fp);
    return count;
}

static void delta(double a, double b, double *max_abs, double *sum, int *n)
{
    double v = fabs(a - b);
    if (v > *max_abs) *max_abs = v;
    *sum += v; (*n)++;
}

static void signed_delta(double diff, double *sum, double *sum_sq,
                         double *max_abs)
{
    *sum += diff;
    *sum_sq += diff * diff;
    if (fabs(diff) > *max_abs) *max_abs = fabs(diff);
}

static int raw_sat(const android_meas_t *m)
{
    int sys = m->constellation_type == 1 ? SYS_GPS :
              m->constellation_type == 5 ? SYS_CMP :
              m->constellation_type == 6 ? SYS_GAL : 0;
    return sys ? satno(sys, m->svid) : 0;
}

static int missing_phase_reason(int sat, uint8_t code, int output_has_code,
                                const android_meas_t *raw, int nraw,
                                const gnss_adapter_t *adapter,
                                const nav_t *nav)
{
    int i, exact = 0, same_band = 0, unsupported = 0, different = 0;
    int valid = 0, finite = 0, nonzero = 0, unc_pass = 0;
    double freq = sat2freq(sat, code, nav);
    for (i = 0; i < nraw; i++) {
        int mapped_sat, slot;
        uint8_t mapped_code;
        int mapped;
        const android_meas_t *m = raw + i;
        if (raw_sat(m) != sat) continue;
        mapped = gnss_adapter_map_signal(m, &mapped_sat, &mapped_code, &slot);
        if (mapped && mapped_sat == sat && mapped_code == code) {
            exact = 1;
            if (!(m->adr_state & ANDROID_ADR_VALID)) continue;
            valid = 1;
            if (!isfinite(m->adr_meters)) continue;
            finite = 1;
            if (fabs(m->adr_meters) <= 1e-4) continue;
            nonzero = 1;
            if (adapter->adr_unc_off ||
                (isfinite(m->adr_uncertainty_meters) &&
                 m->adr_uncertainty_meters >= 0.0 &&
                 m->adr_uncertainty_meters <= adapter->adr_unc_max))
                unc_pass = 1;
        }
        if (freq > 0.0 && isfinite(m->carrier_frequency_hz) &&
            fabs(m->carrier_frequency_hz - freq) < 10000.0) {
            same_band = 1;
            if (!mapped) unsupported = 1;
            else if (mapped_code != code) different = 1;
        }
    }
    if (exact) {
        if (!output_has_code) return MISS_SIGNAL_SELECTION_LOST;
        if (!valid) return MISS_ADR_NOT_VALID;
        if (!finite) return MISS_ADR_NAN;
        if (!nonzero) return MISS_ADR_ZERO;
        if (!unc_pass) return MISS_ADR_UNCERTAINTY;
        return MISS_OTHER;
    }
    if (!same_band) return MISS_NO_RAW_SIGNAL;
    if (unsupported) return MISS_UNSUPPORTED_SIGNAL;
    if (different) return MISS_CODE_MAPPING_DIFFERENT;
    return MISS_OTHER;
}

static void compare_epoch(const obsd_t *obs, int n, const obs_t *rnx,
                          int *cursor, compare_t *c,
                          const android_meas_t *raw, int nraw,
                          const gnss_adapter_t *adapter, const nav_t *nav)
{
    int i, j, k, start;
    double dt;
    if (!rnx || !rnx->data || !rnx->n) return;
    c->epochs++;
    while (*cursor < rnx->n && timediff(rnx->data[*cursor].time, obs[0].time) < -0.01)
        (*cursor)++;
    if (*cursor >= rnx->n ||
        fabs(timediff(rnx->data[*cursor].time, obs[0].time)) > 0.01) {
        c->missing_epochs++; return;
    }
    c->matched_epochs++;
    start = *cursor;
    dt = fabs(timediff(rnx->data[start].time, obs[0].time));
    if (dt > c->max_time) c->max_time = dt;
    if (dt > 0.001) c->time_mismatch++;
    for (i = 0; i < n; i++) {
        const obsd_t *b = NULL;
        for (j = start; j < rnx->n && fabs(timediff(rnx->data[j].time, obs[0].time)) <= 0.01; j++) {
            if (rnx->data[j].sat == obs[i].sat) { b = rnx->data + j; break; }
        }
        if (!b) { c->missing_signals++; continue; }
        for (k = 0; k < NFREQ; k++) {
            if (!obs[i].code[k]) continue;
            if (obs[i].code[k] != b->code[k]) {
                c->code_mismatch++;
                if (!b->code[k]) c->unrepresented_signals++;
                else c->different_code++;
                if (c->code_mismatch <= 12) {
                    char sat[8];
                    satno2id(obs[i].sat, sat);
                    fprintf(stderr, "COMPARE_CODE,%s,F%d,TXT=%s,RINEX=%s\n",
                            sat, k + 1, code2obs(obs[i].code[k]),
                            code2obs(b->code[k]));
                }
                continue;
            }
            c->matched_signals++;
            if (obs[i].P[k] && b->P[k]) delta(obs[i].P[k], b->P[k], &c->max_p, &c->sum_p, &c->np);
            if (obs[i].L[k] && b->L[k]) delta(obs[i].L[k], b->L[k], &c->max_l, &c->sum_l, &c->nl);
            if (obs[i].D[k] && b->D[k]) delta(obs[i].D[k], b->D[k], &c->max_d, &c->sum_d, &c->nd);
            if (obs[i].SNR[k] && b->SNR[k]) delta(obs[i].SNR[k] * SNR_UNIT,
                b->SNR[k] * SNR_UNIT, &c->max_snr, &c->sum_snr, &c->ns);
        }
    }
    /* Reference-centric pass: one reason for every RINEX carrier phase that
     * the TXT path did not reproduce, matched by epoch/satellite/code. */
    for (j = start; j < rnx->n &&
         fabs(timediff(rnx->data[j].time, obs[0].time)) <= 0.01; j++) {
        const obsd_t *b = rnx->data + j;
        if (!(satsys(b->sat, NULL) & (SYS_GPS | SYS_GAL | SYS_CMP))) continue;
        for (k = 0; k < NFREQ; k++) {
            const obsd_t *a = NULL;
            int af = -1, has_code = 0;
            double freq, lambda;
            if (!b->code[k]) continue;
            for (i = 0; i < n; i++) {
                int f;
                if (obs[i].sat != b->sat) continue;
                for (f = 0; f < NFREQ; f++) if (obs[i].code[f] == b->code[k]) {
                    a = obs + i; af = f; has_code = 1; break;
                }
                if (a) break;
            }
            freq = sat2freq(b->sat, b->code[k], nav);
            if (a && a->P[af] && b->P[k]) {
                signed_delta(a->P[af] - b->P[k], &c->p_sum,
                             &c->p_sq, &c->p_max);
                c->strict_p_n++;
            }
            if (freq > 0.0 && a && a->D[af] && b->D[k]) {
                lambda = CLIGHT / freq;
                signed_delta((-a->D[af] + b->D[k]) * lambda,
                             &c->prr_sum, &c->prr_sq, &c->prr_max);
                c->strict_d_n++;
            }
            if (!b->L[k]) continue;
            c->rinex_phase++;
            if (!a || !a->L[af]) {
                int reason = missing_phase_reason(b->sat, b->code[k],
                                                  has_code, raw, nraw,
                                                  adapter, nav);
                c->phase_missing[reason]++;
                continue;
            }
            if (freq <= 0.0) { c->phase_missing[MISS_OTHER]++; continue; }
            lambda = CLIGHT / freq;
            signed_delta((a->L[af] - b->L[k]) * lambda,
                         &c->phase_sum, &c->phase_sq, &c->phase_max);
            c->phase_matched++;
            if (a->LLI[af] == b->LLI[k]) c->lli_match++;
            else c->lli_mismatch++;
        }
    }
}

static void dump_obs(FILE *fp, const obsd_t *obs, int n)
{
    int i, f, week;
    for (i = 0; i < n; i++) for (f = 0; f < NFREQ; f++) {
        char sat[8];
        double tow;
        if (!obs[i].code[f]) continue;
        satno2id(obs[i].sat, sat);
        tow = time2gpst(obs[i].time, &week);
        fprintf(fp, "%d,%.9f,%s,%s,%d,%.4f,%.6f,%.4f,%.3f,%u,%u\n",
                week, tow, sat, code2obs(obs[i].code[f]),
                f + 1, obs[i].P[f], obs[i].L[f], obs[i].D[f],
                obs[i].SNR[f] * SNR_UNIT, obs[i].LLI[f], obs[i].code[f]);
    }
}

static double selected_adr_unc(const android_meas_t *raw, int nraw,
                               int sat, uint8_t code, double phase_m)
{
    int i;
    double best = 1E99, unc = -1.0;
    for (i = 0; i < nraw; i++) {
        int rsat;
        uint8_t rcode;
        double error;
        if (!gnss_adapter_map_signal(raw + i, &rsat, &rcode, NULL) ||
            rsat != sat || rcode != code) continue;
        error = fabs(raw[i].adr_meters - phase_m);
        if (error < best) best = error, unc = raw[i].adr_uncertainty_meters;
    }
    return unc;
}

static int solve_obs_epoch(const obsd_t *obs, int nobs, double adapter_ms,
                           rtk_t *rtk, const nav_t *nav, FILE *dump,
                           FILE *result_file, FILE *state_dump, FILE *use_dump,
                           const obs_t *rnx, int *rnx_cursor,
                           compare_t *comparison, perf_t *perf,
                           const android_meas_t *raw, int nraw,
                           const gnss_adapter_t *adapter,
                           int realtime, gtime_t *first_time,
                           double *first_ms)
{
    double t0, elapsed, pos[3] = {0}, target_ms = 0.0;
    int week, status, i, f;
    double tow;
    if (!nobs) return 0;
    for (i = 0; i < nobs; i++) for (f = 0; f < NFREQ; f++) {
        if (obs[i].P[f]) perf->code_meas++;
        if (obs[i].L[f]) perf->phase_meas++;
        if (obs[i].D[f]) perf->doppler_meas++;
    }
    if (timediff(obs[0].time, nav->peph[0].time) < 0.0 ||
        timediff(obs[0].time, nav->peph[nav->ne - 1].time) > 0.0 ||
        timediff(obs[0].time, nav->pclk[0].time) < 0.0 ||
        timediff(obs[0].time, nav->pclk[nav->nc - 1].time) > 0.0) {
        fprintf(stderr, "PRODUCT_TIME_ERROR,observation outside SP3/CLK coverage\n");
        return -1;
    }
    if (!first_time->time) {
        *first_time = obs[0].time;
        *first_ms = monotonic_ms();
    }
    if (realtime) {
        target_ms = *first_ms + timediff(obs[0].time, *first_time) * 1000.0;
        wait_until(target_ms);
    }
    if (dump) dump_obs(dump, obs, nobs);
    if (rnx) compare_epoch(obs, nobs, rnx, rnx_cursor, comparison,
                           raw, nraw, adapter, nav);
    t0 = monotonic_ms();
    status = rtkpos(rtk, obs, nobs, nav); /* exactly one filter call per epoch */
    elapsed = adapter_ms + monotonic_ms() - t0;
    add_perf(perf, elapsed);
    /* A 1 Hz deadline is missed when the epoch finishes over 1 s after its
     * GNSS-time target, even if the filter itself took less than 1000 ms. */
    if (realtime && monotonic_ms() - target_ms > 1000.0)
        perf->deadline_miss++;
    if (rtk->sol.stat >= 0 && rtk->sol.stat < 8)
        perf->status_count[rtk->sol.stat]++;
    if (rtk->sol.stat) ecef2pos(rtk->sol.rr, pos);
    tow = time2gpst(obs[0].time, &week);
    printf("%d,%.3f,%d,%d,%d,%.9f,%.9f,%.4f,%.3f\n",
           week, tow, nobs, status, rtk->sol.stat,
           pos[0] * R2D, pos[1] * R2D, pos[2], elapsed);
    if (result_file) fprintf(result_file,
        "%d,%.9f,%d,%d,%d,%.9f,%.9f,%.9f,%.12f,%.12f,%.6f,%.6f,%.6f,%.6f,%.3f\n",
        week, tow, rtk->sol.stat, rtk->sol.ns, nobs,
        rtk->sol.rr[0], rtk->sol.rr[1], rtk->sol.rr[2],
        pos[0] * R2D, pos[1] * R2D, pos[2],
        sqrt(fmax(rtk->sol.qr[0], 0.0)),
        sqrt(fmax(rtk->sol.qr[1], 0.0)),
        sqrt(fmax(rtk->sol.qr[2], 0.0)), elapsed);
    if (state_dump) {
        char *buff = (char *)malloc(131072);
        if (buff) {
            int len = rtkoutstat(rtk, 2, buff);
            if (len > 0 && len < 131072) fwrite(buff, 1, len, state_dump);
            free(buff);
        }
    }
    if (use_dump) for (i = 0; i < nobs; i++) for (f = 0; f < NFREQ; f++) {
        int sat = obs[i].sat;
        uint8_t code = obs[i].code[f];
        double freq, phase_m, unc;
        char sid[8];
        if (!code || sat < 1 || sat > MAXSAT) continue;
        freq = sat2freq(sat, code, nav);
        phase_m = freq > 0.0 ? obs[i].L[f] * CLIGHT / freq : 0.0;
        unc = raw ? selected_adr_unc(raw, nraw, sat, code, phase_m) : -1.0;
        satno2id(sat, sid);
        fprintf(use_dump,
            "%d,%.9f,%s,%d,%s,%d,%d,%d,%.3f,%.3f,%.6f,%d,%s,%d,%.6f,%.6f,%.6f\n",
            week, tow, sid, sat, code2obs(code), f + 1,
            obs[i].P[f] != 0.0, obs[i].L[f] != 0.0,
            obs[i].SNR[f] * SNR_UNIT, rtk->ssat[sat - 1].azel[1] * R2D,
            unc, nav->osb_valid[sat - 1][code] != 0,
            receiver_code_bias_slot(satsys(sat, NULL), f, code, rtk->opt.nf),
            rtk->ssat[sat - 1].vsat[f] != 0,
            rtk->ssat[sat - 1].resp[f], rtk->ssat[sat - 1].resc[f], elapsed);
    }
    return 1;
}

static int process_txt_epoch(const android_clock_t *clock,
                             const android_meas_t *raw, int nraw,
                             gnss_adapter_t *adapter, rtk_t *rtk,
                             const nav_t *nav, FILE *dump,
                             FILE *result_file, FILE *state_dump, FILE *use_dump,
                             const obs_t *rnx,
                             int *rnx_cursor, compare_t *comparison,
                             perf_t *perf, int signal_mode,
                             const char *exclude_raw_signal, int realtime,
                             gtime_t *first_time, double *first_ms)
{
    obsd_t obs[MAXOBS];
    android_meas_t excluded_raw[RAW_PER_EPOCH];
    android_meas_t common_raw[RAW_PER_EPOCH];
    gnss_adapter_stats_t stats;
    const android_meas_t *adapter_raw = raw;
    int adapter_nraw = nraw;
    double t0 = monotonic_ms();
    int nobs;
    if (realtime) {
        gtime_t epoch_time;
        if (gnss_adapter_clock_time(clock, &epoch_time)) {
            if (!first_time->time) {
                *first_time = epoch_time;
                *first_ms = monotonic_ms();
            }
            /* Release the original Android epoch at its GNSS-time target,
             * then run the unchanged adapter and one rtkpos() call. */
            wait_until(*first_ms + timediff(epoch_time, *first_time) * 1000.0);
            t0 = monotonic_ms();
        }
    }
    if (exclude_raw_signal) {
        int i, kept = 0;
        uint8_t excluded_code = obs2code(exclude_raw_signal + 2);
        for (i = 0; i < nraw; i++) {
            int sat, sys;
            uint8_t code;
            if (gnss_adapter_map_signal(raw + i, &sat, &code, NULL) &&
                code == excluded_code) {
                sys = satsys(sat, NULL);
                if ((exclude_raw_signal[0] == 'G' && sys == SYS_GPS) ||
                    (exclude_raw_signal[0] == 'E' && sys == SYS_GAL) ||
                    (exclude_raw_signal[0] == 'C' && sys == SYS_CMP)) continue;
            }
            excluded_raw[kept++] = raw[i];
        }
        adapter_raw = excluded_raw;
        adapter_nraw = kept;
    }
    if (signal_mode == 2) {
        adapter_nraw = select_common_raw(clock, adapter_raw, adapter_nraw, rnx,
                                         *rnx_cursor, common_raw);
        adapter_raw = common_raw;
    }
    nobs = gnss_adapter_convert(adapter, clock, adapter_raw, adapter_nraw,
                                obs, MAXOBS, &stats);
    double adapter_ms = monotonic_ms() - t0;
    if (signal_mode == 2)
        nobs = filter_common_only(obs, nobs, rnx, *rnx_cursor);
    add_adr_stats(perf, raw, nraw);
    perf->raw_meas += stats.input_meas;
    perf->unsupported += stats.unsupported;
    perf->duplicates += stats.duplicate;
    perf->phase_rejected_unc += stats.phase_rejected_unc;
    perf->phase_lli_slip += stats.phase_lli_slip;
    perf->lli_android += stats.lli_android;
    perf->lli_clock += stats.lli_clock;
    return solve_obs_epoch(obs, nobs, adapter_ms, rtk, nav, dump,
                           result_file, state_dump, use_dump, rnx,
                           rnx_cursor, comparison, perf, raw, nraw, adapter, realtime,
                           first_time, first_ms);
}

int main(int argc, char **argv)
{
    replay_args_t args;
    nav_t *nav = NULL;
    nav_t *rinex_nav = NULL; /* OBS headers must not mutate solver products */
    rtk_t *rtk = NULL;
    prcopt_t opt;
    obs_t rinex_obs = {0};
    sta_t sta = {0};
    gnss_adapter_t adapter;
    uint8_t rinex_codes[3][NFREQ] = {{0}};
    clock_history_t history = {0};
    perf_t perf = {0};
    compare_t comparison = {0};
    android_clock_t current_clock = {0}, last_clock = {0}, parsed_clock;
    android_meas_t raw[RAW_PER_EPOCH], parsed;
    gtime_t first_time = {0};
    double first_ms = 0.0, p95 = 0.0;
    FILE *input = NULL, *dump = NULL, *result_file = NULL;
    FILE *state_dump = NULL, *use_dump = NULL;
    char line[LINE_SIZE];
    int nraw = 0, epochs = 0, rnx_cursor = 0, rc = 1, fatal = 0;

    if (!parse_args(argc, argv, &args)) { usage(argv[0]); return 2; }
    if (args.txt && !exists(args.txt)) return 2;
    nav = (nav_t *)calloc(1, sizeof(nav_t));
    rtk = (rtk_t *)calloc(1, sizeof(rtk_t));
    if (!nav || !rtk) { fprintf(stderr, "Out of memory\n"); goto done; }
    configure_ppp(&opt);
    traceopen(args.trace);
    tracelevel(args.trace_level);
    if (!load_products(&args, nav, &opt)) {
        fprintf(stderr, "Product loading failed; no PPP run attempted.\n"); goto done;
    }
    if (!args.rinex_source && !scan_signal_products(args.txt, nav)) {
        fprintf(stderr, "SIGNAL_PRODUCT scan failed.\n"); goto done;
    }
    if (args.rinex) {
        rinex_nav = (nav_t *)calloc(1, sizeof(nav_t));
        if (!rinex_nav) goto done;
        if (!exists(args.rinex) || !readrnx(args.rinex, 1, "", &rinex_obs, rinex_nav, &sta)) {
            fprintf(stderr, "RINEX comparison load failed.\n"); goto done;
        }
        sortobs(&rinex_obs);
        if (args.signal_mode == 1 &&
            !learn_rinex_codes(&rinex_obs, rinex_codes)) goto done;
    }
    if (!args.rinex_source) {
        input = fopen(args.txt, "rb");
        if (!input) goto done;
    }
    if (args.dump) {
        dump = fopen(args.dump, "wb");
        if (!dump) goto done;
        fprintf(dump, "week,tow,sat,signal,slot,P_m,L_cycle,D_Hz,CN0_dBHz,LLI,code\n");
    }
    if (args.result_file) {
        result_file = fopen(args.result_file, "wb");
        if (!result_file) goto done;
        fprintf(result_file,
            "week,tow,Q,ns,nobs,x_m,y_m,z_m,lat_deg,lon_deg,h_m,std_x_m,std_y_m,std_z_m,processing_ms\n");
    }
    if (args.state_dump && !(state_dump = fopen(args.state_dump, "wb"))) goto done;
    if (args.use_dump) {
        use_dump = fopen(args.use_dump, "wb");
        if (!use_dump) goto done;
        fprintf(use_dump,
            "week,tow,sat,sat_no,signal,slot,has_P,has_L,cn0_dbhz,elevation_deg,adr_unc_m,osb_valid,rcb_slot,phase_vsat,postfit_code_cache_m,postfit_phase_cache_m,processing_ms\n");
    }
    gnss_adapter_init(&adapter);
    gnss_adapter_set_adr_unc_max(&adapter, args.adr_unc_max, args.adr_unc_off);
    if (args.signal_mode == 1)
        gnss_adapter_set_rinex_codes(&adapter, rinex_codes);
    if (args.signal_mode == 3)
        gnss_adapter_set_ppp_safe(&adapter, 1);
    rtkinit(rtk, &opt);
    puts("week,tow,satellites,rtkpos_return,solution_status,lat_deg,lon_deg,height_m,processing_ms");
    if (args.rinex_source) {
        int i = 0;
        while (i < rinex_obs.n && !fatal) {
            obsd_t obs[MAXOBS];
            gtime_t epoch_time = rinex_obs.data[i].time;
            int nobs = 0, result;
            while (i < rinex_obs.n &&
                   fabs(timediff(rinex_obs.data[i].time, epoch_time)) <= 0.01) {
                if (satsys(rinex_obs.data[i].sat, NULL) & opt.navsys) {
                    if (nobs >= MAXOBS) { fatal = 1; break; }
                    obs[nobs++] = rinex_obs.data[i];
                }
                i++;
            }
            if (fatal) break;
            result = solve_obs_epoch(obs, nobs, 0.0, rtk, nav, dump,
                                     result_file, state_dump, use_dump, NULL,
                                     &rnx_cursor, &comparison, &perf,
                                     NULL, 0, NULL,
                                     args.realtime, &first_time, &first_ms);
            if (result < 0) { fatal = 1; break; }
            epochs += result;
            if (args.max_epochs > 0 && epochs >= args.max_epochs) break;
        }
    }
    else {
        while (fgets(line, sizeof(line), input)) {
            if (!parse_raw(line, &parsed_clock, &parsed, &last_clock)) continue;
            if (nraw && parsed_clock.time_nanos != current_clock.time_nanos) {
                int result = process_txt_epoch(&current_clock, raw, nraw, &adapter, rtk, nav,
                                               dump, result_file, state_dump, use_dump,
                                               args.rinex ? &rinex_obs : NULL, &rnx_cursor,
                                               &comparison, &perf, args.signal_mode,
                                               args.exclude_raw_signal, args.realtime,
                                               &first_time, &first_ms);
                if (result < 0) { fatal = 1; nraw = 0; break; }
                epochs += result;
                nraw = 0;
                if (args.max_epochs > 0 && epochs >= args.max_epochs) break;
            }
            if (!nraw) {
                current_clock = parsed_clock;
                normalize_clock(&history, &current_clock);
            }
            last_clock = parsed_clock;
            if (nraw < RAW_PER_EPOCH) raw[nraw++] = parsed;
        }
        if (nraw && (!args.max_epochs || epochs < args.max_epochs)) {
            int result = process_txt_epoch(&current_clock, raw, nraw, &adapter, rtk, nav,
                                           dump, result_file, state_dump, use_dump,
                                           args.rinex ? &rinex_obs : NULL, &rnx_cursor,
                                           &comparison, &perf, args.signal_mode,
                                           args.exclude_raw_signal, args.realtime,
                                           &first_time, &first_ms);
            if (result < 0) fatal = 1;
            else epochs += result;
        }
    }
    if (perf.n) {
        qsort(perf.ms, perf.n, sizeof(double), cmp_double);
        p95 = perf.ms[(int)ceil(0.95 * perf.n) - 1];
        fprintf(stderr, "PERF,epochs=%d,mean_ms=%.3f,p95_ms=%.3f,max_ms=%.3f,deadline_miss_gt_1000ms=%d,processing_gt_1000ms=%d\n",
                perf.n, perf.sum / perf.n, p95, perf.max, perf.deadline_miss,
                perf.processing_gt_1000ms);
        fprintf(stderr, "STATUS,Q0=%d,Q5=%d,Q6=%d\n",
                perf.status_count[0], perf.status_count[5], perf.status_count[6]);
        fprintf(stderr, "ADAPTER,raw=%d,code=%d,phase=%d,doppler=%d,unsupported=%d,same_slot_discarded=%d\n",
                perf.raw_meas, perf.code_meas, perf.phase_meas,
                perf.doppler_meas, perf.unsupported, perf.duplicates);
        if (!args.rinex_source) {
            fprintf(stderr, "ADR_FILTER,max=%.3f,off=%d,rejected_unc=%d,phase_lli_slip=%d\n",
                    args.adr_unc_max, args.adr_unc_off,
                    perf.phase_rejected_unc, perf.phase_lli_slip);
            print_adr_stats(&perf);
            fprintf(stderr,
                "LLI_CAUSES,android=%d,clock=%d\n",
                perf.lli_android, perf.lli_clock);
        }
    }
    if (args.rinex && !args.rinex_source) {
        fprintf(stderr, "COMPARE,epochs=%d,matched=%d,missing=%d,signals=%d,missing_sat=%d,code_mismatch=%d,time_mismatch=%d,max_time_s=%.6f\n",
                comparison.epochs, comparison.matched_epochs,
                comparison.missing_epochs, comparison.matched_signals,
                comparison.missing_signals, comparison.code_mismatch,
                comparison.time_mismatch, comparison.max_time);
        fprintf(stderr, "COMPARE_CODE_BREAKDOWN,rinex_blank=%d,different_code=%d\n",
                comparison.unrepresented_signals, comparison.different_code);
        fprintf(stderr, "COMPARE_VALUES,P_mean_m=%.6f,P_max_m=%.6f,L_mean_cycle=%.6f,L_max_cycle=%.6f,D_mean_Hz=%.6f,D_max_Hz=%.6f,SNR_mean_dBHz=%.3f,SNR_max_dBHz=%.3f\n",
                comparison.np ? comparison.sum_p / comparison.np : 0.0, comparison.max_p,
                comparison.nl ? comparison.sum_l / comparison.nl : 0.0, comparison.max_l,
                comparison.nd ? comparison.sum_d / comparison.nd : 0.0, comparison.max_d,
                comparison.ns ? comparison.sum_snr / comparison.ns : 0.0, comparison.max_snr);
        fprintf(stderr, "COMPARE_PAIRS,P=%d,L=%d,D=%d,SNR=%d\n",
                comparison.np, comparison.nl, comparison.nd, comparison.ns);
        fprintf(stderr,
            "PHASE_COMPARE,rinex_phase=%d,matched=%d,mean_m=%.6f,rms_m=%.6f,max_abs_m=%.6f\n",
            comparison.rinex_phase, comparison.phase_matched,
            comparison.phase_matched ? comparison.phase_sum / comparison.phase_matched : 0.0,
            comparison.phase_matched ? sqrt(comparison.phase_sq / comparison.phase_matched) : 0.0,
            comparison.phase_max);
        fprintf(stderr,
            "DOPPLER_COMPARE,matched=%d,mean_mps=%.6f,rms_mps=%.6f,max_abs_mps=%.6f\n",
            comparison.strict_d_n,
            comparison.strict_d_n ? comparison.prr_sum / comparison.strict_d_n : 0.0,
            comparison.strict_d_n ? sqrt(comparison.prr_sq / comparison.strict_d_n) : 0.0,
            comparison.prr_max);
        fprintf(stderr,
            "P_COMPARE,matched=%d,mean_m=%.6f,rms_m=%.6f,max_abs_m=%.6f\n",
            comparison.strict_p_n,
            comparison.strict_p_n ? comparison.p_sum / comparison.strict_p_n : 0.0,
            comparison.strict_p_n ? sqrt(comparison.p_sq / comparison.strict_p_n) : 0.0,
            comparison.p_max);
        fprintf(stderr, "LLI_COMPARE,matched=%d,mismatch=%d\n",
                comparison.lli_match, comparison.lli_mismatch);
        fprintf(stderr,
            "PHASE_MISSING,NO_RAW_SIGNAL=%d,UNSUPPORTED_SIGNAL=%d,ADR_NOT_VALID=%d,ADR_ZERO=%d,ADR_NAN=%d,ADR_UNCERTAINTY=%d,SIGNAL_SELECTION_LOST=%d,CODE_MAPPING_DIFFERENT=%d,OTHER=%d\n",
            comparison.phase_missing[MISS_NO_RAW_SIGNAL],
            comparison.phase_missing[MISS_UNSUPPORTED_SIGNAL],
            comparison.phase_missing[MISS_ADR_NOT_VALID],
            comparison.phase_missing[MISS_ADR_ZERO],
            comparison.phase_missing[MISS_ADR_NAN],
            comparison.phase_missing[MISS_ADR_UNCERTAINTY],
            comparison.phase_missing[MISS_SIGNAL_SELECTION_LOST],
            comparison.phase_missing[MISS_CODE_MAPPING_DIFFERENT],
            comparison.phase_missing[MISS_OTHER]);
    }
    rc = perf.n && !fatal ? 0 : 1;
done:
    if (rtk && rtk->x) rtkfree(rtk);
    if (input) fclose(input);
    if (dump) fclose(dump);
    if (result_file) fclose(result_file);
    if (state_dump) fclose(state_dump);
    if (use_dump) fclose(use_dump);
    freeobs(&rinex_obs);
    if (rinex_nav) free_products(rinex_nav);
    if (nav) free_products(nav);
    free(rinex_nav); free(nav); free(rtk); free(perf.ms); free(perf.adr_unc);
    traceclose();
    if (perf.n) fprintf(stderr, "PPP_REJECT,count=%d\n",
                        count_trace_tag(args.trace, "$PPP_REJECT,"));
    return rc;
}
