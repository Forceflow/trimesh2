#ifndef MATHUTIL_H
#define MATHUTIL_H
/*
Szymon Rusinkiewicz
Princeton University

mathutil.h
Miscellaneous math utilities.

Includes the following header files, vaguely portably:
	cstddef, cstdlib, cmath, limits, algorithm, utility
Defines compatibility functions / overloads / classes for pre-C++11 compilers:
	see mathcompat.h for details
Defines float and long double versions of math constants
Defines additional math constants:
	M_2PI, M_PI_3, M_PI_180, M_180_PI
Defines NEXT_MOD3 and PREV_MOD3 macros
Defines likely() and unlikely() macros to help optimize conditionals.
Defines the following utility functions:
	sqr, cube, sgn, sign
Defines the following GLSL-inspired functions:
	radians, degrees, fract, clamp, mix, step, smoothstep
Defines fast, portable random number generators:
	xorshift_rnd, uniform_rnd, normal_rnd
Defines TRIMESH_STATIC_CHECK to do static asserts in pre-11 C++
Defines TRIMESH_DEPRECATED to annotate functions as deprecated in pre-14 C++
*/


// Windows defines min and max as macros, which prevents us from using the
// type-safe versions from std::, as well as interfering with method defns.
// Undefine them, and also define NOMINMAX, preventing future bad definitions.
#ifdef min
# undef min
#endif
#ifdef max
# undef max
#endif
#ifndef NOMINMAX
# define NOMINMAX
#endif


// Some systems require this to get M_PI and friends
#ifndef _USE_MATH_DEFINES
 #define _USE_MATH_DEFINES
#endif


// Standard includes
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>


// Force inlining, even when compiling without optimization.
#if defined(_MSC_VER)
#  define TRIMESH_INLINE __forceinline
#elif defined(__GNUC__) && (__GNUC__ > 3)
#  define TRIMESH_INLINE __inline __attribute__ ((__always_inline__))
#else
#  define TRIMESH_INLINE inline
#endif


// Include math compatibility functions
#include "mathcompat.h"


// Float and long double versions of math constants
#ifndef M_Ef
# define M_Ef           2.7182818f
# define M_LOG2Ef       1.4426950f
# define M_LOG10Ef      0.4342945f
# define M_LN2f         0.6931472f
# define M_LN10f        2.3025851f
# define M_PIf          3.1415927f
# define M_PI_2f        1.5707963f
# define M_PI_4f        0.7853982f
# define M_1_PIf        0.3183099f
# define M_2_PIf        0.6366198f
# define M_2_SQRTPIf    1.1283792f
# define M_SQRT2f       1.4142136f
# define M_SQRT1_2f     0.7071068f
#endif

#ifndef M_El
# define M_El           2.718281828459045235360287471352662498L
# define M_LOG2El       1.442695040888963407359924681001892137L
# define M_LOG10El      0.434294481903251827651128918916605082L
# define M_LN2l         0.693147180559945309417232121458176568L
# define M_LN10l        2.302585092994045684017991454684364208L
# define M_PIl          3.141592653589793238462643383279502884L
# define M_PI_2l        1.570796326794896619231321691639751442L
# define M_PI_4l        0.785398163397448309615660845819875721L
# define M_1_PIl        0.318309886183790671537767526745028724L
# define M_2_PIl        0.636619772367581343075535053490057448L
# define M_2_SQRTPIl    1.128379167095512573896158903121545172L
# define M_SQRT2l       1.414213562373095048801688724209698079L
# define M_SQRT1_2l     0.707106781186547524400844362104849039L
#endif


// Nonstandard math constants

// 2 * pi - not to be confused with M_2_PI, which is 2 / pi
#ifndef M_2PI
# define M_2PIf         6.2831855f
# define M_2PI          6.28318530717958647693
# define M_2PIl         6.283185307179586476925286766559005768L
#endif

// pi / 3
#ifndef M_PI_3
# define M_PI_3f        1.0471976f
# define M_PI_3         1.04719755119659774615
# define M_PI_3l        1.047197551196597746154214461093167628L
#endif

// pi / 180
#ifndef M_PI_180
# define M_PI_180f      0.017453293f
# define M_PI_180       0.0174532925199432957692
# define M_PI_180l      0.01745329251994329576923690768488612713L
#endif

// 180 / pi
#ifndef M_180_PI
# define M_180_PIf      57.2957795f
# define M_180_PI       57.2957795130823208768
# define M_180_PIl      57.29577951308232087679815481410517033L
#endif


// i + 1 and i - 1 modulo 3
// This way of computing it tends to be faster than using %
#define NEXT_MOD3(i) ((i) < 2 ? (i) + 1 : (i) - 2)
#define PREV_MOD3(i) ((i) > 0 ? (i) - 1 : (i) + 2)


// Let gcc optimize conditional branches a bit better...
#ifndef likely
# if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 296)
#  define likely(x)   (__builtin_expect((x), 1))
#  define unlikely(x) (__builtin_expect((x), 0))
# else
#  define likely(x) (x)
#  define unlikely(x) (x)
# endif
#endif


// The remainder is enclosed in namespace trimesh
namespace trimesh {

using ::std::size_t;

#define inline TRIMESH_INLINE


// Square and cube
template <class T>
static inline T sqr(const T &x)
{
	return x * x;
}

template <class T>
static inline T cube(const T &x)
{
	return x * x * x;
}


// Sign of a scalar.  Note that sgn(0) == 0.
template <class T>
static inline T sgn(const T &x)
{
	return (x > 0) ? T(1) : (x < 0) ? T(-1) : T(0);
}


// Alternate spelling of sgn()
template <class T>
static inline T sign(const T &x)
{
	return sgn(x);
}


// Scalar functions based on GLSL - the ones that require vectors are in Vec.h
template <class T>
static inline double radians(const T &x)
{
	return M_PI_180 * x;
}
static inline float radians(const float &x)
{
	return M_PI_180f * x;
}
static inline long double radians(const long double &x)
{
	return M_PI_180l * x;
}

template <class T>
static inline double degrees(const T &x)
{
	return M_180_PI * x;
}
static inline float degrees(const float &x)
{
	return M_180_PIf * x;
}
static inline long double degrees(const long double &x)
{
	return M_180_PIl * x;
}

template <class T>
static inline T fract(const T &x)
{
	return x - floor(x);
}

template <class T, class A, class B>
static inline T clamp(const T &x, const A &a, const B &b)
{
	return (x > T(a)) ? (x < T(b)) ? x : T(b) : T(a);  // returns a if x is NaN
}

template <class T, class A>
static inline T mix(const T &x, const T &y, const A &a)
{
	return (1 - a) * x + a * y;
}

template <class T, class A>
static inline T step(const A &a, const T &x)
{
	return (x < T(a)) ? T(0) : T(1);
}

template <class T, class A, class B>
static inline T smoothstep(const A &a, const B &b, const T &x)
{
	if (unlikely(b <= a))
		return step(x, a);
	T t = (x - a) / (b - a);
	return (t > 0) ? (t < 1) ? (t * t * (3 - 2 * t)) : T(1) : T(0);
}


// Fast, portable pseudorandom number generator using Marsaglia's xorshift.
// Returns random unsigned int.  Pass in 0 to reset.
static inline unsigned xorshift_rnd(unsigned n = ~0u)
{
	static unsigned x = 2463534242u;
	if (likely(n == ~0u)) {
		x ^= x << 13;
		x ^= x >> 17;
		x ^= x << 5;
		return x;
	}
	if (!n)
		x = 2463534242u;
	else
		x = n;
	return x;
}


// Return uniformly-distributed numbers.  Version for integral types.
// If n > 0, returns values in [0, n).
// If n < 0, returns values in [-n , n].  (Note calling convention!)
template <class T>
static inline
typename ::std::enable_if< ::std::is_integral<T>::value, T >::type
uniform_rnd(T n)
{
	if (n > 0) {
		return T(xorshift_rnd() % unsigned(n));
	} else if (n < 0) {
		return T(xorshift_rnd() % unsigned(-2 * n + 1)) + n;
	} else {
		return 0;
	}
}


// Return uniformly-distributed numbers.  Version for floating-point types.
// If n > 0, returns values in (0, n).
// If n < 0, returns values in (-n , n).
template <class T>
static inline
typename ::std::enable_if< ::std::is_floating_point<T>::value, T >::type
uniform_rnd(T n)
{
	if (n > 0) {
		const T scale = T(1) / ::std::numeric_limits<unsigned>::max();
		return T(xorshift_rnd()) * scale * n;
	} else if (n < 0) {
		const T scale = T(2) / ::std::numeric_limits<unsigned>::max();
		return (T(xorshift_rnd()) * scale - T(1)) * n;
	} else {
		return 0;
	}
}


// Without any arguments, return uniform floats in (0,1)
static inline float uniform_rnd()
{
	return uniform_rnd(1.0f);
}


// Returns normal-distributed (Gaussian) random numbers with the given sigma.
// Uses Marsaglia's polar form of the Box-Muller method.
template <class T>
static inline
typename ::std::enable_if< ::std::is_floating_point<T>::value, T >::type
normal_rnd(T sigma)
{
	if (sigma == 0)
		return 0;

	static bool have_saved = false;
	static T saved = 0;
	if (have_saved) {
		have_saved = false;
		return sigma * saved;
	}

	T x, y, r2;
	do {
		x = uniform_rnd(T(-1));
		y = uniform_rnd(T(-1));
		r2 = sqr(x) + sqr(y);
	} while (r2 >= T(1));

	if (r2 == T(0))
		return T(0);

	T g = sqrt(T(-2) * log(r2) / r2);
	saved = x * g;
	have_saved = true;
	return sigma * y * g;
}


// For integral types, call double version
template <class T>
static inline
typename ::std::enable_if< ::std::is_integral<T>::value, double >::type
normal_rnd(T sigma)
{
	return normal_rnd(double(sigma));
}


// Without any arguments, call float version with sigma = 1
static inline float normal_rnd()
{
	return normal_rnd(1.0f);
}


// Boost-like compile-time assertion checking
template <bool X> struct STATIC_ASSERTION_FAILURE;
template <> struct STATIC_ASSERTION_FAILURE<true>
	{ void operator () () {} };
#define TRIMESH_STATIC_CHECK(expr) STATIC_ASSERTION_FAILURE<bool(expr)>()


// Mark a function as deprecated
#if __cplusplus >= 201402L
# define TRIMESH_DEPRECATED [[deprecated]]
#elif defined(__GNUC__)
# define TRIMESH_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
# define TRIMESH_DEPRECATED __declspec(deprecated)
#else
# define TRIMESH_DEPRECATED
#endif

#undef inline

} // namespace trimesh

#endif
