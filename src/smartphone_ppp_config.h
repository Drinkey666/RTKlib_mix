/* Canonical smartphone PPP defaults. Vendor this file unchanged for Android. */
#ifndef SMARTPHONE_PPP_CONFIG_H
#define SMARTPHONE_PPP_CONFIG_H
#include "rtklib.h"
#ifdef __cplusplus
extern "C" {
#endif
#define SMARTPHONE_PPP_PROFILE "pc-baseline-20260930-v1"
/* No products, state allocation or environment overrides are performed here. */
void smartphone_ppp_configure(prcopt_t *opt, solopt_t *sol);
/* Returns 0 with an explanation if compile-time layout/options conflict. */
int smartphone_ppp_validate(const prcopt_t *opt, char *error, size_t capacity);
/* Human-readable effective options for console, Android logcat and trace.
 * The signature covers the printed fields, not product identities or all of
 * prcopt_t. Never interpret equal signatures as proof of equal products. */
unsigned int smartphone_ppp_describe(const prcopt_t *opt, char *out, size_t capacity);
void smartphone_ppp_log(const char *entry, const prcopt_t *opt);
#ifdef __cplusplus
}
#endif
#endif
