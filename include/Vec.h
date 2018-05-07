#ifndef VEC_H
#define VEC_H
/*
Szymon Rusinkiewicz
Princeton University

Vec.h
Class for a constant-length vector, meant to be generally useful for graphics.
Attempts to be similar to the union of <array>, <valarray>, and GLSL vectors,
where convenient.

Creation:
    vec v1;                     // Initialized to (0, 0, 0)
    vec v2(1.23f);              // Initialized to (1.23f, 1.23f, 1.23f)
    vec v3(1, 2, 3);            // Initialized to (1, 2, 3)
    vec v4(v3);                 // Copy constructor

    float farray[3];
    vec v5 = vec(farray);       // Explicit: "vec v5 = farray" won't work

    SomeOtherVectorType v_other;
    vec v6 = vec(v_other);      // Anything for which operator[] is defined

    point p1, p2, p3;           // Same as vec
    Vec<3,double> vd;           // The "vec" used above is Vec<3,float>
    dvec4 vd42;                 // See typedefs below

    v1 = vec::uniform_rnd(3);   // Vector of random numbers in (0,3)
    v1 = vec::normal_rnd(42);   // Gaussian noise vector with sigma = 42

Assignment:
    v1 = v2;
    v2 = 42;                    // Assigns (42, 42, 42)
    v3.fill(42);                // Assigns (42, 42, 42)
    v4.set(1,2,3);
    v5 = farray;

Access:
    float f = v1[0];            // Subscript - not range checked
    f = get<2>(v1);             // Subscript - range checked at compile time
    f = v1.at(42);              // Subscript - range checked at run time
    f = v1.x + v2.y;            // GLSL-like.  Access to 1 component only
    float *fp = v1;             // Implicit conversion to float *
    ivec3 iv(5,6,7);
    int ind = iv.indexof(6);    // Find index of component - returns -1 if none

Vec, Vec/Vec and Vec/scalar operators:
    v1 = -v2;                   // Unary - and +
    if (!v1) {}                 // Check for all components zero
    v1 += v2;                   // Also -=, *=, /=
    v1 *= 2;                    // Also /=, +=, -=
    v1 = v2 + v3 - v4;          // Also *, / (all componentwise)
    v1 = 3.5f * v2 + v3 / 2;    // Also +, - (scalar can come first or second)
    if (v1 == v2) {}            // Also !=, <, >, <=, >=
    std::set<vec>               // This is why we need operator <

Other Vec/Vec and Vec/scalar functions:
    v1.min(v2);                 // Set v1 to min of v1 and v2 - also max, clamp
    v1.clamp(-1, 1);            // Also min, max
    v1 = min(v2, v3);           // Componentwise min - also max, clamp
    v1 = clamp(v2, 0, 1);       // Componentwise clamp - also min, max
    v1.swap(v2);                // Swap - atomic (OpenMP)
    swap(v1, v2);               // Swap - non-atomic

Dot product.  There are four ways to write the dot product of v1 and v2:
    f = v1 DOT v2
    f = v1 ^ v2
    f = dot(v1, v2)
    f = v1.dot(v2)
  SMR prefers the first of these, but all are implemented to allow for
  the range of personal preferences.

Cross product - only in 3 dimensions.  The possible spellings are:
    v1 = v2 CROSS v3;
    v1 = v2 % v3;
    v1 = cross(v2, v3);
    v1 = v2.cross(v3);

Make a vector unit-length.  This operates in-place (not GLSL compatible!):
    normalize(v1);

and this is the version that leaves the original alone, emulating GLSL:
    vec n = normalized(v1);

Functions on vecs:
    f = len(v1);                // Length - also spelled length()
    f = len2(v1);               // Squared length - also length2()
    f = dist(p1, p2);           // Distance - also distance()
    f = dist2(p1, p2);          // Squared distance - also distance2()
    f = angle(v1, v2);          // Angle between vectors
    f = v1.sum();               // From valarrays - see other functions below
    v1 = sin(v2);               // Componentwise - see list of functions below
    v1 = v2.apply(sin);         // Componentwise - any one-argument function
    v1 = reflect(v2, n);        // Reflected vector - n must be unit-length
    v1 = refract(v2, n, 1.5f);  // Refracted vector - n must be unit-length
    v1 = faceforward(v2,v3,v4); // v2 if (v3 DOT v4) > 0, -v2 otherwise
    v1 = trinorm(p1,p2,p3);     // Normal of triangle (area-weighted)

Input/output:
    cout << v1 << endl;         // iostream output in the form (1, 2, 3)
    cin >> v2;                  // iostream input - see below for input format
*/

#include "mathutil.h"
#include <iterator>
#include <stdexcept>
#include <iostream>


#define inline TRIMESH_INLINE


namespace trimesh {


// Storage for Vecs
template <size_t D, class T>
struct Vec_data {
	T v[D];
};

template <class T>
struct Vec_data<1,T> {
	union {
		T v[1];
		struct { T x; };
		struct { T r; };
		struct { T s; };
	};
};

template <class T>
struct Vec_data<2,T> {
	union {
		T v[2];
		struct { T x, y; };
		struct { T r, g; };
		struct { T s, t; };
	};
};

template <class T>
struct Vec_data<3,T> {
	union {
		T v[3];
		struct { T x, y, z; };
		struct { T r, g, b; };
		struct { T s, t, p; };
	};
};

template <class T>
struct Vec_data<4,T> {
	union {
		T v[4];
		struct { T x, y, z, w; };
		struct { T r, g, b, a; };
		struct { T s, t, p, q; };
	};
};


// Utility class for uninitialized constructor
struct Vec_uninitialized {};
#define VEC_UNINITIALIZED ((::trimesh::Vec_uninitialized *) 0)


// Vec class declaration
template <size_t D, class T = float>
class Vec : public Vec_data<D,T> {
protected:
	// Force dependent name lookup for inherited v
	using Vec_data<D,T>::v;

public:
	// Types
	typedef T value_type;
	typedef value_type *pointer;
	typedef const value_type *const_pointer;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef value_type *iterator;
	typedef const value_type *const_iterator;
	typedef ::std::reverse_iterator<iterator> reverse_iterator;
	typedef ::std::reverse_iterator<const_iterator> const_reverse_iterator;
	typedef ::std::size_t size_type;
	typedef ::std::ptrdiff_t difference_type;

	// A type giving the result of any operation (e.g. length) that
	// must produce a floating-point result.  This is double for
	// integral types, else just T itself.
	typedef typename ::std::conditional< ::std::is_integral<T>::value,
	                                     double, T >::type float_type;

public:
	// Constructor for no arguments - everything initialized to zero
	inline Vec() : Vec_data<D,T>() {}

	// Uninitialized constructor - meant mostly for internal use
	inline explicit Vec(Vec_uninitialized *) {}

	// Constructors for 2 - 4 arguments
	inline Vec(const T &x_, const T &y_)
		{ TRIMESH_STATIC_CHECK(D == 2); v[0] = x_; v[1] = y_; }
	inline Vec(const T &x_, const T &y_, const T &z_)
		{ TRIMESH_STATIC_CHECK(D == 3); v[0] = x_; v[1] = y_; v[2] = z_; }
	inline Vec(const T &x_, const T &y_, const T &z_, const T &w_)
		{ TRIMESH_STATIC_CHECK(D == 4); v[0] = x_; v[1] = y_; v[2] = z_; v[3] = w_; }

	// Constructor for 1 scalar argument, which is duplicated into
	// all components.  Explicit.
	template <class S>
	inline explicit Vec(S x_,
		typename ::std::enable_if< ::std::is_arithmetic<S>::value, void >::type * = 0)
	{
		for (size_type i = 0; i < D; i++)
			v[i] = x_;
	}

	// Constructor for 1 argument that's a pointer, array, or
	// anything else that can be accessed using [].  Explicit.
	template <class S>
	inline explicit Vec(const S &v_,
		typename ::std::enable_if< !::std::is_arithmetic<S>::value, void >::type * = 0)
	{
		for (size_type i = 0; i < D; i++)
			v[i] = v_[i];
	}

	// Return a Vec full of uniformly-distributed random numbers
	static inline Vec<D,T> uniform_rnd(T sigma = 1)
	{
		Vec result(VEC_UNINITIALIZED);
		for (size_type i = 0; i < D; i++)
			result[i] = ::trimesh::uniform_rnd(sigma);
		return result;
	}

	// Return a Vec full of normally-distributed random numbers
	static inline Vec<D,T> normal_rnd(T sigma = 1)
	{
		Vec result(VEC_UNINITIALIZED);
		for (size_type i = 0; i < D; i++)
			result[i] = ::trimesh::normal_rnd(sigma);
		return result;
	}

	// Assignment operator equivalents of the one-parameter constructors
	template <class S>
	inline typename ::std::enable_if< ::std::is_arithmetic<S>::value, Vec & >::type
	operator = (S x_)
	{
		for (size_type i = 0; i < D; i++)
			v[i] = x_;
		return *this;
	}

	template <class S>
	inline typename ::std::enable_if< !::std::is_arithmetic<S>::value, Vec & >::type
	operator = (const S &v_)
	{
		for (size_type i = 0; i < D; i++)
			v[i] = v_[i];
		return *this;
	}

	// Using default copy constructor, assignment operator, and destructor

	// Array reference - no bounds checking
	inline reference operator [] (size_type i)
		{ return v[i]; }
	inline reference operator [] (int i)
		{ return v[i]; }
	inline const_reference operator [] (size_type i) const
		{ return v[i]; }
	inline const_reference operator [] (int i) const
		{ return v[i]; }

	// Array reference with run-time bounds checking
	inline reference at(size_type i)
	{
		if (i >= D)
			throw ::std::out_of_range("Vec::at");
		return v[i];
	}
	inline const_reference at(size_type i) const
	{
		if (i >= D)
			throw ::std::out_of_range("Vec::at");
		return v[i];
	}

	// Other accessors, for compatibility with std::array
	inline reference front()
		{ return v[0]; }
	inline const_reference front() const
		{ return v[0]; }
	inline reference back()
		{ return v[D-1]; }
	inline const_reference back() const
		{ return v[D-1]; }

	// Conversion to pointer
	inline operator T * ()
		{ return v; }
	inline operator const T * ()
		{ return v; }
	inline operator const T * () const
		{ return v; }
	inline pointer data()
		{ return v; }
	inline const_pointer data() const
		{ return v; }

	// Iterators
	inline iterator begin()
		{ return v; }
	inline const_iterator begin() const
		{ return v; }
	inline const_iterator cbegin() const
		{ return v; }
	inline iterator end()
		{ return begin() + D; }
	inline const_iterator end() const
		{ return begin() + D; }
	inline const_iterator cend() const
		{ return begin() + D; }
	inline reverse_iterator rbegin()
		{ return reverse_iterator(end()); }
	inline const_reverse_iterator rbegin() const
		{ return const_reverse_iterator(end()); }
	inline const_reverse_iterator crbegin() const
		{ return const_reverse_iterator(end()); }
	inline reverse_iterator rend()
		{ return reverse_iterator(begin()); }
	inline const_reverse_iterator rend() const
		{ return const_reverse_iterator(begin()); }
	inline const_reverse_iterator crend() const
		{ return const_reverse_iterator(begin()); }

	// Capacity
	inline size_type size() const
		{ return D; }
	inline size_type max_size() const
		{ return D; }

	// empty() - check for all components zero.  Note that this definition
	// of empty() is different from std::array, so it's marked deprecated.
	TRIMESH_DEPRECATED
	inline bool empty() const
		{ return !(*this); }

	// Set all components to zero
	inline void clear()
		{ for (size_type i = 0; i < D; i++) v[i] = 0; }

	// Set all elements to some constant
	inline void fill(const value_type &x_)
	{
		for (size_type i = 0; i < D; i++)
			v[i] = x_;
	}

	inline void set(const value_type &x_)
	{
		for (size_type i = 0; i < D; i++)
			v[i] = x_;
	}

	// Set elements to explicit values (only for dimensions 2-4)
	inline void set(const T &x_, const T &y_)
		{ TRIMESH_STATIC_CHECK(D == 2); v[0] = x_; v[1] = y_; }
	inline void set(const T &x_, const T &y_, const T &z_)
		{ TRIMESH_STATIC_CHECK(D == 3); v[0] = x_; v[1] = y_; v[2] = z_; }
	inline void set(const T &x_, const T &y_, const T &z_, const T &w_)
		{ TRIMESH_STATIC_CHECK(D == 4); v[0] = x_; v[1] = y_; v[2] = z_; v[3] = w_; }

	// Componentwise Vec/Vec member operators.
	// (*= and /= included, since some people actually want to do that...)
	inline Vec &operator += (const Vec &v_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] += v_[i];
		return *this;
	}
	inline Vec &operator -= (const Vec &v_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] -= v_[i];
		return *this;
	}
	inline Vec &operator *= (const Vec &v_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] *= v_[i];
		return *this;
	}
	inline Vec &operator /= (const Vec &v_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] /= v_[i];
		return *this;
	}

	// Vec/scalar member operators.
	// (+= and -= included, since some people actually want to do that...)
	inline Vec &operator += (const T &x_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] += x_;
		return *this;
	}
	inline Vec &operator -= (const T &x_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] -= x_;
		return *this;
	}
	inline Vec &operator *= (const T &x_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] *= x_;
		return *this;
	}
	inline Vec &operator /= (const T &x_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp atomic
			v[i] /= x_;
		return *this;
	}

	// Vec/scalar operators - these are friends so that implicit casting
	// can happen on the scalar
	inline friend const Vec operator + (const T &x, const Vec &v)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = x + v[i];
		return result;
	}
	inline friend const Vec operator + (const Vec &v, const T &x)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = v[i] + x;
		return result;
	}
	inline friend const Vec operator - (const T &x, const Vec &v)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = x - v[i];
		return result;
	}
	inline friend const Vec operator - (const Vec &v, const T &x)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = v[i] - x;
		return result;
	}
	inline friend const Vec operator * (const T &x, const Vec &v)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = x * v[i];
		return result;
	}
	inline friend const Vec operator * (const Vec &v, const T &x)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = v[i] * x;
		return result;
	}
	inline friend const Vec operator / (const T &x, const Vec &v)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = x / v[i];
		return result;
	}
	inline friend const Vec operator / (const Vec &v, const T &x)
	{
		using namespace ::std;
		Vec result(VEC_UNINITIALIZED);
		for (size_t i = 0; i < D; i++)
			result[i] = v[i] / x;
		return result;
	}

	// Comparing Vecs and scalars shouldn't work - too easy to make
	// mistakes.  The TRIMESH_STATIC_CHECKs below must depend on D to make
	// clang not attempt to instantiate them when this file is parsed.
	inline friend void operator == (const Vec &, const T &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator == (const T &, const Vec &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator != (const Vec &, const T &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator != (const T &, const Vec &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator >  (const Vec &, const T &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator >  (const T &, const Vec &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator >= (const Vec &, const T &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator >= (const T &, const Vec &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator <  (const Vec &, const T &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator <  (const T &, const Vec &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator <= (const Vec &, const T &)
		{ TRIMESH_STATIC_CHECK(!D); }
	inline friend void operator <= (const T &, const Vec &)
		{ TRIMESH_STATIC_CHECK(!D); }

	// Outside of class: Vec/Vec operators + - * / % ^ << >> == != < > <= >=

	// Vec/Vec in-place min, max, and clamp
	inline Vec &min(const Vec &v_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp critical
			if (v[i] > v_[i]) v[i] = v_[i];
		return *this;
	}
	inline Vec &max(const Vec &v_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp critical
			if (v[i] < v_[i]) v[i] = v_[i];
		return *this;
	}
	inline Vec &clamp(const Vec &a, const Vec &b)
	{
		for (size_type i = 0; i < D; i++) {
#pragma omp critical
			if (v[i] > b[i])
				v[i] = b[i];
			else if (!(v[i] >= a[i]))
				v[i] = a[i];
		}
		return *this;
	}

	// Vec/scalar in-place min, max, and clamp
	inline Vec &min(const T &x_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp critical
			if (v[i] > x_) v[i] = x_;
		return *this;
	}
	inline Vec &max(const T &x_)
	{
		for (size_type i = 0; i < D; i++)
#pragma omp critical
			if (v[i] < x_) v[i] = x_;
		return *this;
	}
	inline Vec &clamp(const T &a, const T &b)
	{
		for (size_type i = 0; i < D; i++) {
#pragma omp critical
			if (v[i] > b)
				v[i] = b;
			else if (!(v[i] >= a))
				v[i] = a;
		}
		return *this;
	}

	// Swap with another vector.  (Also exists as a global function.)
	inline void swap(Vec &v_)
	{
		using namespace ::std;
#pragma omp critical
		for (size_type i = 0; i < D; i++) swap(v[i], v_[i]);
	}

	// Dot product with another vector (also exists as operator ^)
	inline value_type dot(const Vec &v_) const
	{
		value_type total = v[0] * v_[0];
		for (size_type i = 1; i < D; i++)
			total += v[i] * v_[i];
		return total;
	}

	// Cross product with another vector (also exists as operator %)
	inline Vec<3,T> cross(const Vec<3,T> &v_) const
	{
		TRIMESH_STATIC_CHECK(D == 3);
		return Vec<3,T>(v[1] * v_[2] - v[2] * v_[1],
				v[2] * v_[0] - v[0] * v_[2],
				v[0] * v_[1] - v[1] * v_[0]);
	}

	// Some partial compatibility with std::valarray, plus generalizations
	inline value_type sum() const
	{
		value_type total = v[0];
		for (size_type i = 1; i < D; i++)
			total += v[i];
		return total;
	}
	inline value_type sumabs() const
	{
		using namespace ::std;
		value_type total = abs(v[0]);
		for (size_type i = 1; i < D; i++)
			total += abs(v[i]);
		return total;
	}
	inline value_type sumsqr() const
	{
		value_type total = sqr(v[0]);
		for (size_type i = 1; i < D; i++)
			total += sqr(v[i]);
		return total;
	}
	inline float_type avg() const
		{ return float_type(sum()) / D; }
	inline float_type avgabs() const
		{ return float_type(sumabs()) / D; }
	inline float_type mean() const
		{ return float_type(sum()) / D; }
	inline float_type meanabs() const
		{ return float_type(sumabs()) / D; }
	inline float_type rms() const
		{ using namespace ::std;
		  return sqrt(float_type(sumsqr()) / D); }
	inline value_type product() const
	{
		value_type total = v[0];
		for (size_type i = 1; i < D; i++)
			total *= v[i];
		return total;
	}
	inline value_type min() const
	{
		value_type m = v[0];
		for (size_type i = 1; i < D; i++)
			if (v[i] < m)
				m = v[i];
		return m;
	}
	inline value_type minabs() const
	{
		using namespace ::std;
		value_type m = abs(v[0]);
		for (size_type i = 1; i < D; i++) {
			value_type absvi = abs(v[i]);
			if (absvi < m)
				m = absvi;
		}
		return m;
	}
	inline value_type max() const
	{
		value_type m = v[0];
		for (size_type i = 1; i < D; i++)
			if (v[i] > m)
				m = v[i];
		return m;
	}
	inline value_type maxabs() const
	{
		using namespace ::std;
		value_type m = abs(v[0]);
		for (size_type i = 1; i < D; i++) {
			value_type absvi = abs(v[i]);
			if (absvi > m)
				m = absvi;
		}
		return m;
	}
	inline Vec apply(value_type func(value_type)) const
	{
		Vec result(VEC_UNINITIALIZED);
		for (size_type i = 0; i < D; i++)
			result[i] = func(v[i]);
		return result;
	}
	inline Vec apply(value_type func(const value_type&)) const
	{
		Vec result(VEC_UNINITIALIZED);
		for (size_type i = 0; i < D; i++)
			result[i] = func(v[i]);
		return result;
	}
	inline Vec cshift(int n) const
	{
		Vec result(VEC_UNINITIALIZED);
		if (n < 0)
			n = (n % D) + D;
		for (size_type i = 0; i < D; i++)
			result[i] = v[(i+n)%D];
		return result;
	}
	inline Vec shift(int n) const
	{
		using namespace ::std;
		if (unlikely(abs(n) >= D))
			return Vec();
		Vec result; // Must start as zero, so no VEC_UNINITIALIZED
		size_type start = n < 0 ? -n : 0;
		size_type stop = n > 0 ? D - n : D;
		for (size_type i = start; i < stop; i++)
			result[i] = v[i+n];
		return result;
	}

	// Returns index of first element of the Vec that matches the
	// given value exactly.  Returns -1 if not found.
	inline int indexof(const T &x_) const
	{
		for (size_t i = 0; i < D; i++) {
			if (v[i] == x_)
				return i;
		}
		return -1;
	}
}; // class Vec


// Shorthands for particular flavors of Vecs
typedef Vec<3,float>    vec;
typedef Vec<3,float>  point;
typedef Vec<2,float>   vec2;
typedef Vec<3,float>   vec3;
typedef Vec<4,float>   vec4;
typedef Vec<2,float> point2;
typedef Vec<3,float> point3;
typedef Vec<4,float> point4;
typedef Vec<2,int>    ivec2;
typedef Vec<3,int>    ivec3;
typedef Vec<4,int>    ivec4;
typedef Vec<2,double> dvec2;
typedef Vec<3,double> dvec3;
typedef Vec<4,double> dvec4;


// Nonmember operators that take two Vecs
template <size_t D, class T>
static inline const Vec<D,T> operator + (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	Vec<D,T> result(VEC_UNINITIALIZED);
	for (size_t i = 0; i < D; i++)
		result[i] = v1[i] + v2[i];
	return result;
}

template <size_t D, class T>
static inline const Vec<D,T> operator - (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	Vec<D,T> result(VEC_UNINITIALIZED);
	for (size_t i = 0; i < D; i++)
		result[i] = v1[i] - v2[i];
	return result;
}

template <size_t D, class T>
static inline const Vec<D,T> operator * (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	Vec<D,T> result(VEC_UNINITIALIZED);
	for (size_t i = 0; i < D; i++)
		result[i] = v1[i] * v2[i];
	return result;
}

template <size_t D, class T>
static inline const Vec<D,T> operator / (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	Vec<D,T> result(VEC_UNINITIALIZED);
	for (size_t i = 0; i < D; i++)
		result[i] = v1[i] / v2[i];
	return result;
}


// Dot product
template <size_t D, class T>
static inline const T operator ^ (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	T sum = v1[0] * v2[0];
	for (size_t i = 1; i < D; i++)
		sum += v1[i] * v2[i];
	return sum;
}
#define DOT ^

template <size_t D, class T>
static inline const T dot(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return v1 DOT v2;
}


// Cross product
template <class T>
static inline const Vec<3,T> operator % (const Vec<3,T> &v1, const Vec<3,T> &v2)
{
	return Vec<3,T>(v1[1]*v2[2] - v1[2]*v2[1],
			v1[2]*v2[0] - v1[0]*v2[2],
			v1[0]*v2[1] - v1[1]*v2[0]);
}
#define CROSS %

template <class T>
static inline const Vec<3,T> cross(const Vec<3,T> &v1, const Vec<3,T> &v2)
{
	return v1 CROSS v2;
}


// Component-wise equality and inequality.  These return a single bool,
// unlike valarrays, which return one bool per component.
// (#include the usual caveats about comparing floats for equality...)
template <size_t D, class T>
static inline bool operator == (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	for (size_t i = 0; i < D; i++)
		if (v1[i] != v2[i])
			return false;
	return true;
}

template <size_t D, class T>
static inline bool operator != (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	for (size_t i = 0; i < D; i++)
		if (v1[i] != v2[i])
			return true;
	return false;
}


// Comparison by lexicographical ordering - not necessarily useful on its own,
// but necessary in order to put Vecs in sets, maps, etc.
template <size_t D, class T>
static inline bool operator < (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	for (size_t i = 0; i < D; i++) {
		if (v1[i] < v2[i])
			return true;
		else if (v1[i] != v2[i]) // Equivalent to > but catches NaN
			return false;
	}
	return false;
}

template <size_t D, class T>
static inline bool operator > (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return v2 < v1;
}

template <size_t D, class T>
static inline bool operator <= (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	for (size_t i = 0; i < D; i++) {
		if (v1[i] < v2[i])
			return true;
		else if (v1[i] != v2[i])
			return false;
	}
	return true;
}

template <size_t D, class T>
static inline bool operator >= (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return v2 <= v1;
}


// Unary + and -
template <size_t D, class T>
static inline const Vec<D,T> &operator + (const Vec<D,T> &v)
{
	return v;
}

template <size_t D, class T>
static inline const Vec<D,T> operator - (const Vec<D,T> &v)
{
	Vec<D,T> result(VEC_UNINITIALIZED);
	for (size_t i = 0; i < D; i++)
		result[i] = -v[i];
	return result;
}


// Unary ! - check for all elements zero
template <size_t D, class T>
static inline bool operator ! (const Vec<D,T> &v)
{
	for (size_t i = 0; i < D; i++)
		if (v[i] != 0) return false;
	return true;
}


// iostream output.  Formats result as "(1, 2, 3)"
template <size_t D, class T>
static inline ::std::ostream &operator << (::std::ostream &os, const Vec<D,T> &v)

{
	os << "(";
	for (size_t i = 0; i < D - 1; i++)
		os << v[i] << ", ";
	return os << v[D-1] << ")";
}


// iostream input.  Accepts the vec surrounded by (), [], or nothing,
// with components separated by comma, semicolon, or whitespace
template <size_t D, class T>
static inline ::std::istream &operator >> (::std::istream &is, Vec<D,T> &v)
{
	using namespace ::std;
	char c1 = 0, c2 = 0;

	is >> c1;
	if (c1 != '(' && c1 != '[') {
		c1 = 0;
		is.unget();
	}

	is >> v[0];
	for (size_t i = 1; i < D; i++) {
		is >> ws >> c2;
		if (c2 != ',' && c2 != ';')
			is.unget();
		is >> v[i];
	}

	if (c1) {
		is >> ws >> c2;
		if (c1 == '(' && c2 != ')')
			is.setstate(ios::failbit);
		else if (c1 == '[' && c2 != ']')
			is.setstate(ios::failbit);
	}

	if (!is.good())
		v = Vec<D,T>();
	return is;
}


// Vec functions based on GLSL - the scalar ones are in mathutil.h
template <size_t D, class T>
static inline Vec<D,T> faceforward(const Vec<D,T> &N, const Vec<D,T> &I,
                                   const Vec<D,T> &Nref)
{
	return ((Nref DOT I) < 0) ? N : -N;
}

template <size_t D, class T>
static inline Vec<D,T> reflect(const Vec<D,T> &I, const Vec<D,T> &N)
{
	return I - (2 * (N DOT I)) * N;
}

template <size_t D, class T>
static inline Vec<D,T> refract(const Vec<D,T> &I, const Vec<D,T> &N,
                               const T &eta)
{
	using namespace ::std;
	T NdotI = N DOT I;
	T k = 1 - sqr(eta) * (1 - sqr(NdotI));
	if (unlikely(k < 0))
		return Vec<D,T>();
	else
		return eta * I - (eta * NdotI * sqrt(k)) * N;
}


// Squared length
template <size_t D, class T>
static inline const T len2(const Vec<D,T> &v)
{
	T l2 = sqr(v[0]);
	for (size_t i = 1; i < D; i++)
		l2 += sqr(v[i]);
	return l2;
}


// Length
template <size_t D, class T>
static inline const typename Vec<D,T>::float_type
len(const Vec<D,T> &v)
{
	using namespace ::std;
	return sqrt(len2(v));
}


// Alternate, GLSL-compatible spelling of len2() and len()
template <size_t D, class T>
static inline const T length2(const Vec<D,T> &v)
{
	return len2(v);
}

template <size_t D, class T>
static inline const typename Vec<D,T>::float_type
length(const Vec<D,T> &v)
{
	return len(v);
}


// Squared distance
template <size_t D, class T>
static inline const T dist2(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	T d2 = sqr(v2[0] - v1[0]);
	for (size_t i = 1; i < D; i++)
		d2 += sqr(v2[i] - v1[i]);
	return d2;
}


// Distance
template <size_t D, class T>
static inline const typename Vec<D,T>::float_type
dist(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	using namespace ::std;
	return sqrt(dist2(v1, v2));
}


// Alternate, GLSL-compatible spelling of dist2() and dist()
template <size_t D, class T>
static inline const T distance2(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return dist2(v1, v2);
}

template <size_t D, class T>
static inline const typename Vec<D,T>::float_type
distance(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return dist(v1, v2);
}


// In-place normalization to unit length.  For historical reasons, this is
// incompatible with the GLSL normalize() - that's implemented as normalized()
template <size_t D, class T>
static inline void normalize(Vec<D,T> &v)
{
	T l = len(v);
	if (likely(l > 0)) {
		l = 1 / l;
		for (size_t i = 0; i < D; i++)
			v[i] *= l;
	} else {
		// Make sure we have sane output for length 0 and NaN
		for (size_t i = 0; i < D - 1; i++)
			v[i] = 0;
		v[D-1] = 1;
	}
}


// Returns a normalized vector while leaving the original alone
template <size_t D, class T>
static inline Vec<D,T> normalized(const Vec<D,T> &v)
{
	Vec<D,T> w(v);
	normalize(w);
	return w;
}


// Area-weighted triangle face normal
template <class T>
static inline Vec<3,T> trinorm(const Vec<3,T> &v0, const Vec<3,T> &v1, const Vec<3,T> &v2)
{
	return T(0.5) * ((v1 - v0) CROSS (v2 - v0));
}


// Angle between two vectors
template <size_t D, class T>
static inline const typename Vec<D,T>::float_type
angle(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	using namespace ::std;
	typedef typename Vec<D,T>::float_type FT;

	// Formula from section 12 of
	// http://www.cs.berkeley.edu/~wkahan/Mindless.pdf
	Vec<D,FT> x(v1), y(v2);
	x *= len(v2);
	y *= len(v1);
	return 2 * atan2(len(x-y), len(x+y));
}


} // namespace trimesh


// Generic macros for declaring 1-, 2-, and 3- argument
// componentwise functions on Vecs.
#define VEC_DECLARE_ONEARG(name) \
 template < ::std::size_t D, class T > \
 static inline trimesh::Vec<D,T> name(const trimesh::Vec<D,T> &v) \
 { \
	using namespace ::std; \
	using namespace ::trimesh; \
	Vec<D,T> result(VEC_UNINITIALIZED); \
	for (size_t i = 0; i < D; i++) \
		result[i] = name(v[i]); \
	return result; \
 }

// Vector-scalar, scalar-vector, and componentwise vector-vector versions
#define VEC_DECLARE_TWOARG_VS(name) \
 template < ::std::size_t D, class T > \
 static inline trimesh::Vec<D,T> name(const trimesh::Vec<D,T> &v, const T &a) \
 { \
	using namespace ::std; \
	using namespace ::trimesh; \
	Vec<D,T> result(VEC_UNINITIALIZED); \
	for (size_t i = 0; i < D; i++) \
		result[i] = name(v[i], a); \
	return result; \
 }
#define VEC_DECLARE_TWOARG_SV(name) \
 template < ::std::size_t D, class T > \
 static inline trimesh::Vec<D,T> name(const T &a, const trimesh::Vec<D,T> &v) \
 { \
	using namespace ::std; \
	using namespace ::trimesh; \
	Vec<D,T> result(VEC_UNINITIALIZED); \
	for (size_t i = 0; i < D; i++) \
		result[i] = name(a, v[i]); \
	return result; \
 }
#define VEC_DECLARE_TWOARG_VV(name) \
 template < ::std::size_t D, class T > \
 static inline trimesh::Vec<D,T> name(const trimesh::Vec<D,T> &v, const trimesh::Vec<D,T> &w) \
 { \
	using namespace ::std; \
	using namespace ::trimesh; \
	Vec<D,T> result(VEC_UNINITIALIZED); \
	for (size_t i = 0; i < D; i++) \
		result[i] = name(v[i], w[i]); \
	return result; \
 }

#define VEC_DECLARE_THREEARG_VSS(name) \
 template < ::std::size_t D, class T > \
 static inline trimesh::Vec<D,T> name(const trimesh::Vec<D,T> &v, const T &a, const T &b) \
 { \
	using namespace ::std; \
	using namespace ::trimesh; \
	Vec<D,T> result(VEC_UNINITIALIZED); \
	for (size_t i = 0; i < D; i++) \
		result[i] = name(v[i], a, b); \
	return result; \
 }
#define VEC_DECLARE_THREEARG_SSV(name) \
 template < ::std::size_t D, class T > \
 static inline trimesh::Vec<D,T> name(const T &a, const T &b, const trimesh::Vec<D,T> &v) \
 { \
	using namespace ::std; \
	using namespace ::trimesh; \
	Vec<D,T> result(VEC_UNINITIALIZED); \
	for (size_t i = 0; i < D; i++) \
		result[i] = name(a, b, v[i]); \
	return result; \
 }
#define VEC_DECLARE_THREEARG_VVV(name) \
 template < ::std::size_t D, class T > \
 static inline trimesh::Vec<D,T> name(const trimesh::Vec<D,T> &v, const trimesh::Vec<D,T> &w, const trimesh::Vec<D,T> &x) \
 { \
	using namespace ::std; \
	using namespace ::trimesh; \
	Vec<D,T> result(VEC_UNINITIALIZED); \
	for (size_t i = 0; i < D; i++) \
		result[i] = name(v[i], w[i], x[i]); \
	return result; \
 }


// The following is the list of functions in C89 and C++98,
// PLUS the ones in mathcompat.h (which are POSIX / C99 / C++11)
// MINUS frexp, ldexp, and modf (which have irregular calling conventions).
// They are supposed to be in namespace std, but Visual Studio and some
// older compilers also declare them in the global namespace.
// In the name of compatibility, we put them in both global and std.
VEC_DECLARE_ONEARG(abs)
VEC_DECLARE_ONEARG(acos)
VEC_DECLARE_ONEARG(acosh)
VEC_DECLARE_ONEARG(asin)
VEC_DECLARE_ONEARG(asinh)
VEC_DECLARE_ONEARG(atan)
VEC_DECLARE_ONEARG(atanh)
VEC_DECLARE_TWOARG_VV(atan2)
VEC_DECLARE_ONEARG(cbrt)
VEC_DECLARE_ONEARG(ceil)
VEC_DECLARE_ONEARG(cos)
VEC_DECLARE_ONEARG(cosh)
VEC_DECLARE_ONEARG(exp)
VEC_DECLARE_ONEARG(exp2)
VEC_DECLARE_ONEARG(expm1)
VEC_DECLARE_ONEARG(fabs)
VEC_DECLARE_TWOARG_VS(fdim)
VEC_DECLARE_TWOARG_SV(fdim)
VEC_DECLARE_TWOARG_VV(fdim)
VEC_DECLARE_ONEARG(floor)
VEC_DECLARE_TWOARG_VS(fmod)
VEC_DECLARE_TWOARG_VV(fmod)
VEC_DECLARE_TWOARG_VS(hypot)
VEC_DECLARE_TWOARG_SV(hypot)
VEC_DECLARE_TWOARG_VV(hypot)
VEC_DECLARE_ONEARG(log)
VEC_DECLARE_ONEARG(log10)
VEC_DECLARE_ONEARG(log1p)
VEC_DECLARE_ONEARG(log2)
VEC_DECLARE_TWOARG_VS(pow)
VEC_DECLARE_TWOARG_SV(pow)
VEC_DECLARE_TWOARG_VV(pow)
VEC_DECLARE_ONEARG(round)
VEC_DECLARE_ONEARG(sin)
VEC_DECLARE_ONEARG(sinh)
VEC_DECLARE_ONEARG(sqrt)
VEC_DECLARE_ONEARG(tan)
VEC_DECLARE_ONEARG(tanh)
VEC_DECLARE_ONEARG(trunc)


// Inject into namespace std
namespace std {
	using ::abs;
	using ::acos;
	using ::asin;
	using ::atan;
	using ::atan2;
	using ::cbrt;
	using ::ceil;
	using ::cos;
	using ::cosh;
	using ::exp;
	using ::fabs;
	using ::floor;
	using ::fmod;
	using ::hypot;
	using ::log;
	using ::log10;
	using ::pow;
	using ::round;
	using ::sin;
	using ::sinh;
	using ::sqrt;
	using ::tan;
	using ::tanh;
	using ::trunc;

	// These are only in namespace std.
	VEC_DECLARE_TWOARG_VS(min)
	VEC_DECLARE_TWOARG_SV(min)
	VEC_DECLARE_TWOARG_VV(min)
	VEC_DECLARE_TWOARG_VS(max)
	VEC_DECLARE_TWOARG_SV(max)
	VEC_DECLARE_TWOARG_VV(max)

	// Swap two Vecs.  Not atomic, unlike class method.
	template <size_t D, class T>
	static inline void swap(const ::trimesh::Vec<D,T> &v1, const ::trimesh::Vec<D,T> &v2)
	{
		for (size_t i = 0; i < D; i++)
			swap(v1[i], v2[i]);
	}

	// Get an element with compile-time bounds checking
	template <size_t I, size_t D, class T>
	static inline T &get(::trimesh::Vec<D,T> &v)
	{
		using namespace ::trimesh;
		TRIMESH_STATIC_CHECK(I < D);
		return v[I];
	}
	template <size_t I, size_t D, class T>
	static inline const T &get(const ::trimesh::Vec<D,T> &v)
	{
		using namespace ::trimesh;
		TRIMESH_STATIC_CHECK(I < D);
		return v[I];
	}
} // namespace std


// These are new functions declared in namespace trimesh (in mathutil.h)
namespace trimesh {
	VEC_DECLARE_ONEARG(sqr)
	VEC_DECLARE_ONEARG(cube)
	VEC_DECLARE_ONEARG(sgn)
	VEC_DECLARE_ONEARG(radians)
	VEC_DECLARE_ONEARG(degrees)
	VEC_DECLARE_ONEARG(fract)
	VEC_DECLARE_THREEARG_VSS(clamp)
	VEC_DECLARE_THREEARG_VVV(clamp)
	VEC_DECLARE_TWOARG_SV(step)
	VEC_DECLARE_TWOARG_VV(step)
	VEC_DECLARE_THREEARG_SSV(smoothstep)
	VEC_DECLARE_THREEARG_VVV(smoothstep)
} // namespace trimesh


#undef VEC_DECLARE_ONEARG
#undef VEC_DECLARE_TWOARG_VS
#undef VEC_DECLARE_TWOARG_SV
#undef VEC_DECLARE_TWOARG_VV
#undef VEC_DECLARE_THREEARG_VSS
#undef VEC_DECLARE_THREEARG_SSV
#undef VEC_DECLARE_THREEARG_VVV

#undef inline

#endif
