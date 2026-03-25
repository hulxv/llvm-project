.. _trig_algorithm:

==========================================
Trigonometric and Inverse Trig Algorithms
==========================================

.. default-role:: math

This page describes the algorithms used in LLVM libc for trigonometric
functions (sin, cos, tan) and inverse trigonometric functions (asin, acos,
atan, and their ``*pi`` variants).

For the general framework (two-phase computation, Ziv's test, etc.), see
:ref:`math_algorithms`.

.. contents:: Table of Contents
   :depth: 3
   :local:


Trigonometric Functions
=======================

sin, cos, tan (double)
----------------------

All three functions share a common range reduction strategy and differ only in
how the reduced values are combined.

Range reduction
~~~~~~~~~~~~~~~

The goal is to compute:

.. math::

   k = \operatorname{round}\!\left(\frac{128 \, x}{\pi}\right),
   \qquad
   y = x - k \cdot \frac{\pi}{128}

so that `|y| \leq \pi/256`.

**Small inputs** (`|x| < 2^{16}`):  A double-precision constant
`128/\pi \approx` ``0x1.45f306dc9c883p5`` is used for the initial
multiplication, and a 3-part representation of `\pi/128`:

.. math::

   \frac{\pi}{128} \approx p_0 + p_1 + p_2

with `p_0 =` ``-0x1.921fb544p-6``, `p_1 =` ``-0x1.0b4611a6p-40``,
`p_2 =` ``-0x1.3198a2e037073p-75``, is used for the subtraction
`y = x - k(p_0 + p_1 + p_2)`.  The error is bounded by `2^{-111}`.

**Large inputs** (`|x| \geq 2^{16}`):  An adaptive-precision lookup table
of `128/\pi` with 64 entries of 4 double-precision chunks each is indexed
by the exponent of `x`.  This handles arbitrarily large inputs without loss
of accuracy.

Lookup tables
~~~~~~~~~~~~~

A table of 256 DoubleDouble entries stores `\sin(k\pi/128)` for
`k = 0, 1, \ldots, 255`.  Cosine values are obtained by a shifted index:

.. math::

   \cos\!\left(\frac{k\pi}{128}\right)
   = \sin\!\left(\frac{(k + 64)\pi}{128}\right)

A memory-saving variant uses only 65 entries and exploits the symmetry
`\sin(k\pi/128) = \sin((128 - k)\pi/128)`.

For the accurate phase (Float128), a separate 65-entry table of 128-bit
values is used.

Polynomial evaluation
~~~~~~~~~~~~~~~~~~~~~

After range reduction, `\sin(y)` and `\cos(y)` are evaluated simultaneously
by ``sincos_eval``:

.. math::

   \sin(y) &\approx y - \frac{y^3}{3!} + \frac{y^5}{5!} - \frac{y^7}{7!}
   \quad \text{(degree-7 Taylor)} \\
   \cos(y) &\approx 1 - \frac{y^2}{2!} + \frac{y^4}{4!} - \frac{y^6}{6!}
   + \frac{y^8}{8!}
   \quad \text{(degree-8 Taylor)}

These are evaluated in DoubleDouble arithmetic with errors bounded by
`2^{-72}` for sin and `2^{-81}` for cos.

Combination via addition formulas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**sin(x):**

.. math::

   \sin(x) = \sin(y)\cos\!\left(\frac{k\pi}{128}\right)
   + \cos(y)\sin\!\left(\frac{k\pi}{128}\right)

**cos(x):**

.. math::

   \cos(x) = \cos(y)\cos\!\left(\frac{k\pi}{128}\right)
   - \sin(y)\sin\!\left(\frac{k\pi}{128}\right)

**tan(x):**

.. math::

   \tan(x) =
   \frac{\sin\!\left(\frac{k\pi}{128}\right) + \tan(y)\cos\!\left(\frac{k\pi}{128}\right)}
        {\cos\!\left(\frac{k\pi}{128}\right) - \tan(y)\sin\!\left(\frac{k\pi}{128}\right)}

where `\tan(y)` is approximated by a degree-9 Taylor polynomial:

.. math::

   \tan(y) \approx y + \frac{y^3}{3} + \frac{2y^5}{15}
   + \frac{17y^7}{315} + \frac{62y^9}{2835}


sinf, cosf, tanf (float)
-------------------------

The single-precision variants use a similar strategy but with simpler
range reduction (mod `\pi/32` instead of `\pi/128`), 32-entry lookup
tables, and Sollya-generated minimax polynomials of degree 6–7 for sin/cos
and degree 9 for tan.

Exceptional values are handled via precomputed lookup tables rather than
Ziv's test, since the smaller output space makes exhaustive verification
feasible.


Inverse Trigonometric Functions
===============================

asin, asinpi (double)
---------------------

Domain splitting
~~~~~~~~~~~~~~~~

The domain `[-1, 1]` is split into three regions:

1. **Small inputs** (`|x| < 2^{-26}`):

   .. math::

      \operatorname{asin}(x) \approx x,
      \qquad
      \operatorname{asinpi}(x) \approx \frac{x}{\pi}

   The approximation has relative error less than `2^{-54}`.

2. **|x| < 0.5** — Direct polynomial:

   .. math::

      \operatorname{asin}(x) \approx x \cdot P(x^2)

   where `P` is a minimax polynomial approximating
   `\operatorname{asin}(x)/x` (or `\operatorname{asin}(x)/(\pi x)` for
   asinpi).

   The argument `x^2` is further decomposed as `x^2 = k/64 + \delta` where
   `k = \operatorname{round}(64 \cdot x^2)` and `|\delta| < 1/128`.  A
   64-entry lookup table stores precomputed polynomial values at each `k/64`,
   and a correction polynomial handles `\delta`.

3. **0.5 ≤ |x| < 1** — Range reduction via complementary angle:

   .. math::

      \operatorname{asin}(x)
      = \frac{\pi}{2} - 2\operatorname{asin}\!\left(\sqrt{\frac{1-x}{2}}\right)

   Let `u = (1 - |x|)/2` and `v = \sqrt{u}`.  Then:

   .. math::

      \operatorname{asin}(x) = \frac{\pi}{2} - 2v \cdot P(u)

   and for asinpi:

   .. math::

      \operatorname{asinpi}(x) = \frac{1}{2} - 2v \cdot P_\pi(u)

   The square root `v` is refined using Newton–Raphson correction terms:

   .. math::

      v_{\text{lo}} = \frac{u - v_{\text{hi}}^2}{2 v_{\text{hi}}},
      \qquad
      v_{\text{ll}} = -v_{\text{lo}} \cdot \frac{u - v_{\text{hi}}^2}{4u}

   These extra terms ensure the Float128 fallback path has sufficient
   precision.

4. **|x| = 1**: returns `\pm\pi/2` (asin) or `\pm 0.5` (asinpi) exactly.

5. **|x| > 1**: domain error (``EDOM``, ``FE_INVALID``, returns NaN).

Subnormal handling (asinpi)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For very small inputs (`|x| < 2^{-511}`), computing `x^2` would underflow.
Instead, `x/\pi` is computed directly in Float128 arithmetic.

At the subnormal/normal boundary (`|x/\pi| \approx 2^{-1022}`), IEEE 754
"after rounding" tininess detection requires special care.  The
implementation examines the top 53 bits of the Float128 mantissa to
determine whether rounding would carry the result from subnormal to normal
(`2^{-1022}`), checking each rounding mode separately.  See the source in
``asinpi.h`` for the detailed logic.


acos, acospi (double)
---------------------

These use the same polynomial infrastructure as asin, with different
combination constants:

.. math::

   \operatorname{acos}(x) = \frac{\pi}{2} - \operatorname{asin}(x)

For `0.5 \leq x < 1`:

.. math::

   \operatorname{acos}(x) = 2\operatorname{asin}\!\left(\sqrt{\frac{1-x}{2}}\right)

For `-1 < x \leq -0.5`:

.. math::

   \operatorname{acos}(x) = \pi - 2\operatorname{asin}\!\left(\sqrt{\frac{1+x}{2}}\right)


atan, atan2 (double)
--------------------

Domain splitting
~~~~~~~~~~~~~~~~

1. **|x| < 2^-26**: `\operatorname{atan}(x) \approx x`.

2. **|x| < 1**: Table-based reduction with polynomial correction.
   Let `k = \operatorname{round}(64 \cdot |x|)`.  Then:

   .. math::

      \operatorname{atan}(x)
      = \operatorname{atan}\!\left(\frac{k}{64}\right)
      + \operatorname{atan}\!\left(\frac{|x| - k/64}{1 + k|x|/64}\right)

   A 65-entry lookup table stores `\operatorname{atan}(k/64)` as DoubleDouble
   values.  The correction term is small enough for a degree-13 minimax
   polynomial.

3. **|x| ≥ 1**: Reciprocal transformation.

   .. math::

      \operatorname{atan}(x)
      = \operatorname{sign}(x)\left(\frac{\pi}{2}
      - \operatorname{atan}\!\left(\frac{1}{|x|}\right)\right)

   This maps the input to `[0, 1]`, where the table-based method applies.

4. **|x| ≥ 2^53**:
   `\operatorname{atan}(x) \approx \operatorname{sign}(x) \cdot \pi/2`.


Float variants (asinf, acosf, atanf)
-------------------------------------

The single-precision inverse trig functions use:

- Degree-20 minimax polynomials (asinf, acosf) or degree-10 Taylor series (atanf)
- Double-precision intermediates for accuracy
- Precomputed exceptional value tables instead of Ziv's test
- 16-entry lookup tables (atanf) or no tables (asinf/acosf with direct polynomial)


Source files
============

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - File
     - Contents
   * - ``sin.h``
     - `\sin(x)` (double)
   * - ``cos.h``
     - `\cos(x)` (double)
   * - ``sincos.h``
     - simultaneous `\sin` and `\cos`
   * - ``tan.h``
     - `\tan(x)` (double)
   * - ``sincos_eval.h``
     - shared sin/cos polynomial evaluation
   * - ``range_reduction_double_common.h``
     - lookup tables, constants, Float128 tables
   * - ``range_reduction_double_fma.h``
     - FMA-based range reduction
   * - ``range_reduction_double_nofma.h``
     - non-FMA range reduction fallback
   * - ``asin.h``
     - `\operatorname{asin}(x)` (double)
   * - ``asinpi.h``
     - `\operatorname{asinpi}(x)` (double)
   * - ``acos.h``
     - `\operatorname{acos}(x)` (double)
   * - ``atan.h``
     - `\operatorname{atan}(x)` (double)
   * - ``asin_utils.h``
     - shared polynomials and tables for asin/acos family
   * - ``atan_utils.h``
     - shared polynomials and tables for atan family

All source files are located under ``libc/src/__support/math/``.