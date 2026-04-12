//===-- Unittests for acospi ----------------------------------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "hdr/fenv_macros.h"
#include "src/math/acospi.h"
#include "test/UnitTest/FPMatcher.h"
#include "test/UnitTest/Test.h"

using LlvmLibcAcospiTest = LIBC_NAMESPACE::testing::FPTest<double>;

TEST_F(LlvmLibcAcospiTest, SpecialNumbers) {
  EXPECT_FP_EQ_WITH_EXCEPTION_ALL_ROUNDING(aNaN, LIBC_NAMESPACE::acospi(sNaN),
                                           FE_INVALID);
  EXPECT_FP_EQ_ALL_ROUNDING(aNaN, LIBC_NAMESPACE::acospi(aNaN));
  EXPECT_FP_EQ_ALL_ROUNDING(0.5, LIBC_NAMESPACE::acospi(zero));
  EXPECT_FP_EQ_ALL_ROUNDING(0.5, LIBC_NAMESPACE::acospi(neg_zero));
  EXPECT_FP_EQ_ALL_ROUNDING(aNaN, LIBC_NAMESPACE::acospi(inf));
  EXPECT_FP_EQ_ALL_ROUNDING(aNaN, LIBC_NAMESPACE::acospi(neg_inf));
  EXPECT_FP_EQ_ALL_ROUNDING(aNaN, LIBC_NAMESPACE::acospi(2.0));
  EXPECT_FP_EQ_ALL_ROUNDING(aNaN, LIBC_NAMESPACE::acospi(-2.0));
  EXPECT_FP_EQ(0.0, LIBC_NAMESPACE::acospi(1.0));
  EXPECT_FP_EQ(1.0, LIBC_NAMESPACE::acospi(-1.0));
}

#ifdef LIBC_TEST_FTZ_DAZ

using namespace LIBC_NAMESPACE::testing;

TEST_F(LlvmLibcAcospiTest, FTZMode) {
  ModifyMXCSR mxcsr(FTZ);

  EXPECT_TRUE(0.5 == LIBC_NAMESPACE::acospi(min_denormal));
}

TEST_F(LlvmLibcAcospiTest, DAZMode) {
  ModifyMXCSR mxcsr(DAZ);

  EXPECT_TRUE(0.5 == LIBC_NAMESPACE::acospi(min_denormal));
}

TEST_F(LlvmLibcAcospiTest, FTZDAZMode) {
  ModifyMXCSR mxcsr(FTZ | DAZ);

  EXPECT_TRUE(0.5 == LIBC_NAMESPACE::acospi(min_denormal));
}

#endif
