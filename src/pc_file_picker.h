#ifndef PC_FILE_PICKER_H
#define PC_FILE_PICKER_H
#include <stddef.h>
#define PC_PPP_PATH_CAP 1024
enum {
    PC_OBS, PC_NAV, PC_SP3, PC_CLK, PC_BIA, PC_IONEX,
    PC_VMF_FIRST, PC_VMF_SECOND, PC_OROGRAPHY, PC_ATX, PC_OUTPUT,
    PC_FILE_COUNT
};
enum { PC_SOURCE_RINEX, PC_SOURCE_TXT };
/* PC_OBS is exactly one active input, never an OBS/TXT pair. */
typedef struct {
    char path[PC_FILE_COUNT][PC_PPP_PATH_CAP];
    int source;
} pc_ppp_files_t;
extern const char *pc_file_keys[PC_FILE_COUNT];
/* Windows-only UI; 1=start, 0=cancel, -1=UI error. Never invokes PPP. */
int pc_file_picker(pc_ppp_files_t *files, char *error, size_t capacity);
/* File/path validation only: does not establish product temporal coverage. */
int pc_files_validate(const pc_ppp_files_t *files, char *error, size_t capacity);
#endif
