#ifndef GNSS_SIGNAL_POLICY_H
#define GNSS_SIGNAL_POLICY_H

#include "rtklib.h"
#include <stdint.h>

/* Slot is zero-based, matching obsd_t P/L/D/code arrays. A positive priority
 * permits the signal in the current PPP/RCB model. This policy is independent
 * of a RINEX OBS file and does not alter satellite bias corrections. */
int gnss_ppp_signal_priority(int sys, int slot, uint8_t code);
int gnss_ppp_signal_allowed(int sys, int slot, uint8_t code);
const char *gnss_ppp_signal_reason(int sys, int slot, uint8_t code);

#endif
