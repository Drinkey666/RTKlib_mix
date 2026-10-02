#ifndef PC_LIVE_VIEW_H
#define PC_LIVE_VIEW_H
#include "pc_txt_replay.h"
/* Main thread owns the window; one background worker owns the entire filter. */
int pc_live_view_run(const pc_ppp_files_t *files, const prcopt_t *options,
                     const char *trace_path, int realtime, int max_epochs);
#endif
