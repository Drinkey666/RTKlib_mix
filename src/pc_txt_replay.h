#ifndef PC_TXT_REPLAY_H
#define PC_TXT_REPLAY_H
#include "pc_file_picker.h"
#include "rtklib.h"
#include "gnss_replay.h"
int pc_txt_replay(const pc_ppp_files_t *files, const prcopt_t *options,
                  const char *trace_path, int realtime, int max_epochs);
int pc_txt_replay_observed(const pc_ppp_files_t *files, const prcopt_t *options,
                           const char *trace_path, int realtime, int max_epochs,
                           const gnss_replay_observer_t *observer);
#endif
