#ifndef MATHCOMPAT_H
#define MATHCOMPAT_H
/*
Szymon Rusinkiewicz
Princeton University

mathcompat.h
Defines math compatibility functions / overloads for pre-C++11 compilers.

This isn't the full C++11 set - it includes:
	acosh, asinh, atanh, cbrt, exp2, expm1, fdim, hypot, log1p, log2,
	round, trunc
plus the following type traits:
	integral_constant, true_type, false_type, is_integral,
	is_floating_point, is_arithmetic, conditional, enable_if,
*/


// Get rid of any bad #defines.
#ifdef hypot
# undef hypot
#endif
#ifdef round
# undef round
#endif
#ifdef trunc
# undef trunc
#endif


// Now define versions for early versions of MSVS
#define inline TRIMESH_INLINE

#if defined(_MSC_VER)

// MSVS documentation implies these were only defined in MSVS 2013, but
// internet comments suggest they were defined in VS 2010.
#if _MSC_VER < 1600
static inline float hypotf(float x, float y)
	{ return _hypotf(x, y); }
static inline double hypot(double x, double y)
	{ return _hypot(x, y); }
static inline long double hypotl(long double x, long double y)
	{ return _hypotl(x, y); }
#endif

// Available since MSVS 2012
#if _MSC_VER < 1700
static inline float acoshf(float x)
	{ return log(x + sqrt(x * x - 1)); }
static inline double acosh(double x)
	{ return log(x + sqrt(x * x - 1)); }
static inline long double acoshl(long double x)
	{ return log(x + sqrt(x * x - 1)); }

static inline float asinhf(float x)
	{ return log(x + sqrt(x * x + 1)); }
static inline double asinh(double x)
	{ return log(x + sqrt(x * x + 1)); }
static inline long double asinhl(long double x)
	{ return log(x + sqrt(x * x + 1)); }

static inline float atanhf(float x)
	{ return (log(1 + x) - log(1 - x)) / 2; }
static inline double atanh(double x)
	{ return (log(1 + x) - log(1 - x)) / 2; }
static inline long double atanhl(long double x)
	{ return (log(1 + x) - log(1 - x)) / 2; }

static inline float cbrtf(float x)
	{ return (x < 0) ? -pow(-x, 1.0f / 3.0f) : pow(x, 1.0f / 3.0f); }
static inline double cbrt(double x)
	{ return (x < 0) ? -pow(-x, 1.0 / 3.0) : pow(x, 1.0 / 3.0); }
static inline long double cbrtl(long double x)
	{ return (x < 0) ? -pow(-x, 1.0L / 3.0L) : pow(x, 1.0L / 3.0L); }

static inline float roundf(float x)
	{ return (x < 0) ? ceil(x - 0.5f) : floor(x + 0.5f); }
static inline double round(double x)
	{ return (x < 0) ? ceil(x - 0.5) : floor(x + 0.5); }
static inline long double roundl(long double x)
	{ return (x < 0) ? ceil(x - 0.5L) : floor(x + 0.5L); }

static inline float truncf(float x)
	{ return (x < 0) ? ceil(x) : floor(x); }
static inline double trunc(double x)
	{ return (x < 0) ? ceil(x) : floor(x); }
static inline long double truncl(long double x)
	{ return (x < 0) ? ceil(x) : floor(x); }
#endif

// Available since MSVS 2013
#if _MSC_VER < 1800
static inline float exp2f(float x)
	{ return exp(0.6931472f * x); }
static inline double exp2(double x)
	{ return exp(0.69314718055994530942 * x); }
static inline long double exp2l(long double x)
	{ return exp(0.693147180559945309417232121458176568L * x); }

static inline float expm1f(float x)
	{ return (fabs(x) < 6e-3f) ? (x + 0.5f * x * x) : (exp(x) - 1); }
static inline double expm1(double x)
	{ return (fabs(x) < 7e-6) ? (x + 0.5 * x * x) : (exp(x) - 1); }
static inline long double expm1l(long double x)
	{ return (fabs(x) < 5e-7L) ? (x + 0.5L * x * x) : (exp(x) - 1); }

static inline float log1pf(float x)
	{ return (fabs(x) < 4e-3f) ? (x - 0.5f * x * x) : log(x + 1); }
static inline double log1p(double x)
	{ return (fabs(x) < 6e-6) ? (x - 0.5 * x * x) : log(x + 1); }
static inline long double log1pl(long double x)
	{ return (fabs(x) < 4e-7L) ? (x - 0.5L * x * x) : log(x + 1); }

static inline float log2f(float x)
	{ return log(1.4426950f * x); }
static inline double log2(double x)
	{ return log(1.4426950408889634074 * x); }
static inline long double log2l(long double x)
	{ return log(1.442695040888963407359924681001892137L * x); }

static inline float fdimf(float x, float y)
	{ return x > y ? x - y : 0.0f; }
static inline double fdim(double x, double y)
	{ return x > y ? x - y : 0.0; }
static inline long double fdiml(long double x, long double y)
	{ return x > y ? x - y : 0.0L; }
#endif

#endif // defined(_MSC_VER)


// Finally, define overloads for pre-C++11 compilers
// (except for libc++, which has these already, even in C++03 mode)
#if (defined(_MSC_VER) && _MSC_VER < 1700) || \
    (!defined(_MSC_VER) && !defined(_LIBCPP_VERSION) && __cplusplus < 201103L)

#define DECLARE_OVERLOADS_ONEARG(name) \
template <class T> \
static inline double name(T x) \
	{ using namespace ::std; return name(double(x)); } \
static inline float name(float x) \
	{ using namespace ::std; return name ## f (x); } \
static inline long double name(long double x) \
	{ using namespace ::std; return name ## l (x); }

#define DECLARE_OVERLOADS_TWOARG(name) \
template <class T> \
static inline double name(T x, T y) \
	{ using namespace ::std; return name(double(x), double(y)); } \
static inline float name(float x, float y) \
	{ using namespace ::std; return name ## f (x, y); } \
static inline long double name(long double x, long double y) \
	{ using namespace ::std; return name ## l (x, y); }

DECLARE_OVERLOADS_ONEARG(acosh)
DECLARE_OVERLOADS_ONEARG(asinh)
DECLARE_OVERLOADS_ONEARG(atanh)
DECLARE_OVERLOADS_ONEARG(cbrt)
DECLARE_OVERLOADS_ONEARG(exp2)
DECLARE_OVERLOADS_ONEARG(expm1)
DECLARE_OVERLOADS_TWOARG(fdim)
DECLARE_OVERLOADS_TWOARG(hypot)
DECLARE_OVERLOADS_ONEARG(log1p)
DECLARE_OVERLOADS_ONEARG(log2)
DECLARE_OVERLOADS_ONEARG(round)
DECLARE_OVERLOADS_ONEARG(trunc)

#undef DECLARE_OVERLOADS_ONEARG
#undef DECLARE_OVERLOADS_TWOARG

// Inject into namespace std
namespace std {
	using ::acosh;
	using ::asinh;
	using ::atanh;
	using ::cbrt;
	using ::exp2;
	using ::expm1;
	using ::fdim;
	using ::hypot;
	using ::log1p;
	using ::log2;
	using ::round;
	using ::trunc;
}
#endif // _MSC_VER < 1800 || __cplusplus < 201103L


// A few type_traits for pre-C++11
#if _MSC_VER >= 1600 || defined(_LIBCPP_VERSION) || __cplusplus >= 201103L
# include <type_traits>
#elif !defined(_LIBCPP_TYPE_TRAITS) && !defined(_TR1_TYPE_TRAITS) && !defined(_GLIBCXX_TYPE_TRAITS)
namespace std {

template<class T, T v> struct integral_constant { static const T value = v; };

typedef integral_constant<bool, true> true_type;
typedef integral_constant<bool, false> false_type;

template <class T> struct is_integral : integral_constant<bool,
	::std::numeric_limits<T>::is_integer> {};
template <class T, int N> struct is_integral<T[N]> : false_type {};
template <class T> struct is_integral<const T> : is_integral<T> {};

template <class T> struct is_floating_point : integral_constant<bool,
	::std::numeric_limits<T>::is_specialized
	&& !::std::numeric_limits<T>::is_integer> {};
template <class T, int N> struct is_floating_point<T[N]> : false_type {};
template <class T> struct is_floating_point<const T> : is_floating_point<T> {};

template <class T> struct is_arithmetic : integral_constant<bool,
	is_integral<T>::value || is_floating_point<T>::value> {};

template<bool B, class T, class F> struct conditional { typedef T type; };
template<class T, class F> struct conditional<false, T, F> { typedef F type; };

template<bool B, class T> struct enable_if {};
template<class T> struct enable_if<true, T> { typedef T type; };

} // namespace std
#endif // _MSC_VER >= 1600 || __cplusplus >= 201103L

#undef inline

#endif
