#ifndef GNSS_RAW_H
#define GNSS_RAW_H

#include <stdint.h>

#define ANDROID_ADR_VALID 1u
#define ANDROID_ADR_RESET 2u
#define ANDROID_ADR_SLIP 4u
#define ANDROID_ADR_HALF_RESOLVED 8u
#define ANDROID_ADR_HALF_REPORTED 16u

/* Android GnssClock/GnssMeasurement values, independent of TXT and RTKLIB.
 * Empty optional TXT fields are represented by has_* rather than zero. */
typedef struct {
    int64_t time_nanos;
    int64_t full_bias_nanos;
    double bias_nanos;
    int has_full_bias;
    int hardware_clock_discontinuity_count;
    int64_t clock_correction_nanos; /* explicit correction for silent clock jumps */
    int clock_discontinuity; /* synthetic discontinuity at this epoch */
} android_clock_t;

typedef struct {
    int constellation_type; /* Android: GPS=1, BeiDou=5, Galileo=6 */
    int svid;
    double time_offset_nanos;
    uint32_t state;
    int64_t received_sv_time_nanos;
    double received_sv_time_uncertainty_nanos;
    double cn0_dbhz;
    double pseudorange_rate_mps;
    double pseudorange_rate_uncertainty_mps;
    uint32_t adr_state;
    double adr_meters;
    double adr_uncertainty_meters;
    double carrier_frequency_hz;
    char code_type[8];
} android_meas_t;

#endif
