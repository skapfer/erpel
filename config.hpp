#ifndef ERPEL_CONFIG_HPP_INCLUDED
#define ERPEL_CONFIG_HPP_INCLUDED

#ifndef _REENTRANT
#error "_REENTRANT missing"
#endif
#include <pthread.h>

#ifdef __GNUC__
#define does_not_return __attribute__ ((noreturn))
#else
#define does_not_return
#endif // __GNUC__

#ifdef __INTEL_COMPILER_BUILD_DATE
// evaluate operands in unspecified
#pragma warning (disable: 981)
// external function definition with no prior declaration
#pragma warning (disable: 1418)
// value copied to temporary, reference to temporary used
#pragma warning (disable: 383)
// floating pnt cmp unreliable
#pragma warning (disable: 1572)
#endif // __INTEL_COMPILER_BUILD_DATE

#endif // ERPEL_CONFIG_HPP_INCLUDED
