#ifndef BOX_H
#define BOX_H
/*
Szymon Rusinkiewicz
Princeton University

Box.h
Templated axis-aligned bounding boxes - meant to be used with Vec.h
*/

#include "Vec.h"
#include "strutil.h"
#include <iomanip>
#include <fstream>


#define inline TRIMESH_INLINE


namespace trimesh {

template <size_t D, class T = float>
class Box {
public:
	// Types
	typedef T value_type;
	typedef value_type *pointer;
	typedef const value_type *const_pointer;
	typedef value_type &reference;
	typedef const value_type &const_reference;
	typedef Vec<D,T> point_type;
	typedef typename Vec<D,T>::float_type float_type;

	// Public (!) members
	point_type min, max;
	bool valid;

	// Construct as empty
	inline Box() : valid(false)
		{}

	// Construct from a single point
	inline Box(const point_type &p) : min(p), max(p), valid(true)
		{}

	// Construct from two points
	inline Box(const point_type &p1, const point_type &p2)
		: min(p1), max(p2), valid(true)
	{
		using namespace ::std;
		for (size_t i = 0; i < D; i++) {
			if (min[i] > max[i])
				swap(min[i], max[i]);
		}
	}

	// Mark invalid
	inline void clear()
		{ valid = false; }

	// Return center point, (vector) diagonal, and (scalar) radius
	inline point_type center() const
	{
		if (unlikely(!valid))
			return point_type();
		return float_type(0.5) * (min + max);
	}
	inline point_type size() const
	{
		if (unlikely(!valid))
			return point_type();
		return max - min;
	}
	inline float_type radius() const
	{
		if (unlikely(!valid))
			return 0;
		return float_type(0.5) * dist(min, max);
	}

	// Grow a bounding box to encompass a point or another Box
	inline Box<D,T> &operator += (const point_type &p)
	{
		if (likely(valid)) {
			for (size_t i = 0; i < D; i++) {
				if (p[i] < min[i])
					min[i] = p[i];
				else if (p[i] > max[i])
					max[i] = p[i];
			}
		} else {
			min = p;
			max = p;
			valid = true;
		}
		return *this;
	}
	inline Box<D,T> &operator += (const Box<D,T> &b)
	{
		if (likely(valid)) {
			for (size_t i = 0; i < D; i++) {
				if (b.min[i] < min[i])
					min[i] = b.min[i];
				if (b.max[i] > max[i])
					max[i] = b.max[i];
			}
		} else {
			min = b.min;
			max = b.max;
			valid = true;
		}
		return *this;
	}

	inline friend const Box operator + (const Box &b, const point_type &p)
		{ return Box(b) += p; }
	inline friend const Box operator + (const point_type &p, const Box &b)
		{ return Box(b) += p; }
	inline friend const Box operator + (const Box &b1, const Box &b2)
		{ return Box(b1) += b2; }

	// Does a Box contain, or at least touch, a point?
	inline bool contains(const point_type &p) const
	{
		if (unlikely(!valid))
			return false;
		for (size_t i = 0; i < D; i++) {
			if (p[i] < min[i] || p[i] > max[i])
				return false;
		}
		return true;
	}

	// Does a Box contain another Box?
	inline bool contains(const Box &b) const
	{
		if (unlikely(!valid || !b.valid))
			return false;
		for (size_t i = 0; i < D; i++) {
			if (b.min[i] < min[i] || b.max[i] > max[i])
				return false;
		}
		return true;
	}

	// Do two Boxes intersect, or at least touch?
	inline bool intersects(const Box &b)
	{
		if (unlikely(!valid || !b.valid))
			return false;
		for (size_t i = 0; i < D; i++) {
			if (b.max[i] < min[i] || b.min[i] > max[i])
				return false;
		}
		return true;
	}

	// Read a Box from a file
	inline bool read(const ::std::string &filename)
	{
		using namespace ::std;
		fstream f(filename.c_str());
		Box<D,T> B;
		f >> B;
		f.close();
		if (f.good()) {
			*this = B;
			return true;
		}
		return false;
	}

	// Write a Box to a file
	inline bool write(const ::std::string &filename) const
	{
		using namespace ::std;
		const int digits = 2 + numeric_limits<T>::digits10;
		ofstream f(filename.c_str());
		f << setprecision(digits) << *this;
		f.close();
		return f.good();
	}

	// iostream operators
	inline friend ::std::ostream &operator << (::std::ostream &os, const Box &b)
	{
		using namespace ::std;
		const size_t n = b.min.size();
		for (size_t i = 0; i < n; i++)
			os << b.min[i] << (i == n-1 ? "\n" : " ");
		for (size_t i = 0; i < n; i++)
			os << b.max[i] << (i == n-1 ? "\n" : " ");
		return os;
	}
	inline friend ::std::istream &operator >> (::std::istream &is, Box &b)
	{
		using namespace ::std;
		const size_t n = b.min.size();
		for (size_t i = 0; i < n; i++)
			is >> b.min[i];
		for (size_t i = 0; i < n; i++)
			is >> b.max[i];
		for (size_t i = 0; i < n; i++)
			if (b.min[i] > b.max[i])
				swap(b.min[i], b.max[i]);
		b.valid = is.good();
		return is;
	}
}; // class Box


typedef Box<3,float>    box;
typedef Box<2,float>   box2;
typedef Box<3,float>   box3;
typedef Box<4,float>   box4;
typedef Box<2,int>    ibox2;
typedef Box<3,int>    ibox3;
typedef Box<4,int>    ibox4;
typedef Box<2,double> dbox2;
typedef Box<3,double> dbox3;
typedef Box<4,double> dbox4;


// Equality and inequality.  An invalid box compares unequal to anything.
template <size_t D, class T>
static inline bool operator == (const Box<D,T> &b1, const Box<D,T> &b2)
{
	return b1.valid && b2.valid && b1.min == b2.min && b1.max == b2.max;
}

template <size_t D, class T>
static inline bool operator != (const Box<D,T> &b1, const Box<D,T> &b2)
{
	return !(b1 == b2);
}


// (In-)validity testing
template <size_t D, class T>
static inline bool operator ! (const Box<D,T> &b)
{
	return !b.valid;
}


// Generate a .bbox filename from an input (scan) filename
static inline ::std::string bboxname(const ::std::string &filename)
{
	return replace_ext(filename, "bbox");
}

} // namespace trimesh

#undef inline

#endif
