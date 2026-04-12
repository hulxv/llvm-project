//===-- Double-precision acospi function ----------------------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "src/math/acospi.h"
#include "src/__support/math/acospi.h"

namespace LIBC_NAMESPACE_DECL {

LLVM_LIBC_FUNCTION(double, acospi, (double x)) { return math::acospi(x); }

} // namespace LIBC_NAMESPACE_DECL
