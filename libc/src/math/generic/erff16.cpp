//===-- Half-precision erf(x) function ----------------------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "src/__support/FPUtil/FPBits.h"
#include "src/__support/FPUtil/PolyEval.h"
#include "src/__support/FPUtil/multiply_add.h"
#include "src/__support/FPUtil/nearest_integer.h"
#include "src/__support/macros/optimization.h"
#include "src/math/generic/explogxf.h"

namespace LIBC_NAMESPACE {

LLVM_LIBC_FUNCTION(float16, erff16, (float16 x)) {
  using FPBits = fputil::FPBits<float16>;
  FPBits xbits(x);
  bool is_neg = xbits.is_neg();

  // Handle special cases
  if (LIBC_UNLIKELY(xbits.is_inf_or_nan())) {
    if (xbits.is_nan())
      return x;
    // erf(+inf) = 1, erf(-inf) = -1
    return is_neg ? FPBits::neg_one().get_val() : FPBits::one().get_val();
  }

  if (LIBC_UNLIKELY(xbits.is_zero()))
    return x; // erf(0) = 0, erf(-0) = -0

  float16 abs_x = xbits.abs();

  /*
   * * |x| in [0, 1.0]:
   * Use polynomial approximation: erf(x) ≈ x * P(x^2)
   *
   * This is based on the Taylor series:
   * erf(x) = (2/sqrt(pi)) * (x - x³/3 + x⁵/10 - x⁷/42 + x⁹/216 - ...)
   *        = x * (2/sqrt(pi)) * (1 - x^2/3 + x⁴/10 - x⁶/42 + ...)
   *
   * We use Chebyshev polynomial approximation for better accuracy
   */
  if (abs_x <= 1.0f16) {
    float16 x_sq = x * x;

    // Polynomial coefficients for P(x^2) optimized for float16 range
    // These are derived from Chebyshev approximation of (2/sqrt(pi)) * erf(x)/x
    constexpr float16 COEFFS[] = {
        0x1.20dd76p+0f16,  // 2/sqrt(pi) ≈ 1.1284 (c0)
        -0x1.812746p-2f16, // -0.3761 (c1 * 2/sqrt(pi) / 3)
        0x1.ce2f21p-4f16,  // 0.1128 (c2)
        -0x1.b52b6p-6f16,  // -0.0268 (c3)
        0x1.7d4932p-9f16,  // 0.0046 (c4)
    };

    // P(x^2) = c0 + c1*x^2 + c2*x⁴ + c3*x⁶ + c4*x⁸
    float16 poly = fputil::polyeval(x2, COEFFS[0], COEFFS[1], COEFFS[2],
                                    COEFFS[3], COEFFS[4]);

    float16 result = x * poly;
    return result;
  }

  /*
   * * |x| in (1.0, 3.0]:
   * Use asymptotic expansion: erfc(x) ≈ exp(-x^2) / (xsqrt(pi)) * (1 +
   * correction) Then erf(x) = 1 - erfc(x)
   *
   * The asymptotic series is:
   * erfc(x) ~ (exp(-x^2) / (x*sqrt(pi)) * (1 - 1/(2x^2) + 3/(4x⁴) - 15/(8x⁶) +
   * ...)
   */
  if (abs_x <= 3.0f16) {
    float16 x2 = abs_x * abs_x;
    float16 inv_x2 = 1.0f16 / x2;

    // Asymptotic series coefficients
    // For float16, we only need first few terms due to limited precision
    static constexpr float16 ASYMP_COEFFS[] = {
        1.0f16,    // 1
        -0.5f16,   // -1/2
        0.75f16,   // 3/4
        -1.875f16, // -15/8
    };

    float16 series = fputil::polyeval(inv_x2, ASYMP_COEFFS[0], ASYMP_COEFFS[1],
                                      ASYMP_COEFFS[2], ASYMP_COEFFS[3]);

    // Compute exp(-x^2) using optimized exponential
    float16 exp_neg_x2 = fputil::exp(-x2);

    // erfc(x) = exp(-x^2) / (xsqrt(pi)) * series
    static constexpr float16 inv_sqrt_pi =
        0x1.20dd76p-1f16; // 1/sqrt(pi) ≈ 0.5642
    float16 erfc_val = (exp_neg_x2 * inv_sqrt_pi * series) / abs_x;

    float16 erf_val = 1.0f16 - erfc_val;
    return is_neg ? -erf_val : erf_val;
  }

  /*
   * * |x| > 3.0:
   *   erf(x) ≈ +/-1, with exponentially small correction
   *   For x > 3.0, erf(x) is very close to 1
   *   The difference from 1 is approximately exp(-x^2) / (xsqrt(pi))
   */
  if (abs_x >= 4.0f16) {
    return is_neg ? -1.0f16 : 1.0f16;
  }

  float16 x2 = abs_x * abs_x;
  float16 exp_neg_x2 = fputil::exp(-x2);
  static constexpr float16 inv_sqrt_pi = 0x1.20dd76p-1f16;
  float16 correction = (exp_neg_x2 * inv_sqrt_pi) / abs_x;

  float16 erf_val = 1.0f16 - correction;
  return is_neg ? -erf_val : erf_val;
}

} // namespace LIBC_NAMESPACE