#include "gnss_signal_policy.h"
#include <string.h>

int gnss_ppp_signal_priority(int sys, int slot, uint8_t code)
{
    const char *sig = code2obs(code);
    if (!sig || !*sig) return 0;
    if (sys == SYS_GPS) {
        if (slot == 0 && !strcmp(sig, "1C")) return 100;
        if (slot == 2 && !strcmp(sig, "5Q")) return 100;
    }
    else if (sys == SYS_GAL) {
        if (slot == 0 && !strcmp(sig, "1C")) return 100;
        if (slot == 2 && !strcmp(sig, "5Q")) return 100;
    }
    else if (sys == SYS_CMP) {
        if (slot == 0 && !strcmp(sig, "2I")) return 100;
        if (slot == 1 && !strcmp(sig, "7I")) return 100;
        if (slot == 2 && !strcmp(sig, "5P")) return 100;
    }
    return 0;
}

int gnss_ppp_signal_allowed(int sys, int slot, uint8_t code)
{
    return gnss_ppp_signal_priority(sys, slot, code) > 0;
}

const char *gnss_ppp_signal_reason(int sys, int slot, uint8_t code)
{
    const char *sig = code2obs(code);
    if (gnss_ppp_signal_allowed(sys, slot, code)) return "PPP_POLICY_ALLOWED";
    if (sys == SYS_CMP && slot == 2 && sig && !strcmp(sig, "5Q"))
        return "NO_BDS_C5Q_RCB";
    if (sys == SYS_CMP && slot == 0 && sig && !strcmp(sig, "1P"))
        return "BDS_2I_REFERENCE_SLOT_PRIORITY";
    if (sys == SYS_GAL && slot == 1)
        return "GAL_E5B_ISOLATED_IN_CURRENT_PPP";
    return "NOT_IN_CURRENT_PPP_SIGNAL_POLICY";
}
