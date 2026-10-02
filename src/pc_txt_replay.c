/* One observation source only: no RINEX OBS argument is ever constructed. */
#include "pc_txt_replay.h"
#include "gnss_replay.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int pc_txt_replay_observed(const pc_ppp_files_t *files, const prcopt_t *options,
                           const char *trace_path, int realtime, int max_epochs,
                           const gnss_replay_observer_t *observer)
{
    char state[PC_PPP_PATH_CAP], dump[PC_PPP_PATH_CAP], limit[24];
    char *args[48];
    const char *level = getenv("RTK_PPP_TRACE_LEVEL");
    int n = 0, trace_level = level ? atoi(level) : 3;
    if (!files || files->source != PC_SOURCE_TXT) return 2;
    /* The existing replay accepts 2/3. Reject instead of silently changing it. */
    if (trace_level < 2 || trace_level > 3) {
        fprintf(stderr, "TXT replay RTK_PPP_TRACE_LEVEL must be 2 or 3\n");
        return 2;
    }
    snprintf(state, sizeof(state), "%s.stat", files->path[PC_OUTPUT]);
    snprintf(dump, sizeof(dump), "%s.obs.csv", files->path[PC_OUTPUT]);
    snprintf(limit, sizeof(limit), "%d", max_epochs);
    args[n++] = "pc-txt-replay";
#define ARG(key, value) do { args[n++] = (char *)(key); args[n++] = (char *)(value); } while (0)
    ARG("--source", realtime ? "txt-realtime" : "txt-fast");
    ARG("--txt", files->path[PC_OBS]);
    ARG("--nav", files->path[PC_NAV]);
    ARG("--sp3", files->path[PC_SP3]);
    ARG("--clk", files->path[PC_CLK]);
    ARG("--bia", files->path[PC_BIA]);
    ARG("--ionex", files->path[PC_IONEX]);
    ARG("--vmf0", files->path[PC_VMF_FIRST]);
    ARG("--vmf6", files->path[PC_VMF_SECOND]);
    ARG("--orog", files->path[PC_OROGRAPHY]);
    ARG("--atx", files->path[PC_ATX]);
    ARG("--signal-mode", "ppp-safe");
    ARG("--adr-unc-max", "1.0");
    ARG("--trace", trace_path);
    ARG("--trace-level", trace_level == 2 ? "2" : "3");
    ARG("--result-file", files->path[PC_OUTPUT]);
    ARG("--state-dump", state);
    ARG("--obs-dump", dump);
    if (max_epochs) ARG("--max-epochs", limit);
#undef ARG
    args[n] = NULL;
    printf("PPP_SOURCE,TXT,%s,signal_mode=ppp-safe,rinex_obs=NONE\n",
           realtime ? "realtime-1x" : "fast-diagnostic");
    fflush(stdout);
    return gnss_replay_run_observed(n, args, options, observer);
}

int pc_txt_replay(const pc_ppp_files_t *files, const prcopt_t *options,
                  const char *trace_path, int realtime, int max_epochs)
{
    return pc_txt_replay_observed(files, options, trace_path, realtime, max_epochs, NULL);
}
