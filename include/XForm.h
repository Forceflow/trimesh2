#ifndef XFORM_H
#define XFORM_H
/*
Szymon Rusinkiewicz
Princeton University

XForm.h
Transformations (represented internally as column-major 4x4 matrices)

Supports the following operations:
	xform xf, xf2;			// Initialized to the identity
	XForm<float> xf3;		// Just "xform" is XForm<double>
	xf(1,2)=3.0;			// Access by row/column
	xf[4]=5.0;			// Access in column-major order
	xf=xform::trans(u,v,w);		// An xform that translates
	xf=xform::rot(ang,ax);		// An xform rotating by ang around ax
	xf=xform::rot_into(dir1,dir2);	// An xform rotating dir1 into dir2
	xf=xform::scale(s);		// An xform that scales
	xf=xform::scale(sx,sy,sz);	// An xform that scales
	xf=xform::ortho(l,r,b,t,n,f);   // Like GLortho
	xf=xform::frustum(l,r,b,t,n,f);	// Like GLfrustum
	xf=xform::outer(v, w);		// Outer product of two vecs
	xf=xform::fromarray(A);		// From A[3][3] or A[4][4]
	glMultMatrixd(xf);		// Conversion to column-major array
	xf = xf * xf2;			// Matrix-matrix multiplication
	vec v = xf * vec(1,2,3);	// Matrix-vector multiplication
	xf = inv(xf) + transp(xf2);	// Inverse, transpose
	xf2 = rot_only(xf);		// Just the upper 3x3 of xf
	xf2 = trans_only(xf);		// Just the translation of xf
	xf2 = norm_xf(xf);		// Inverse transpose, no translation
	invert(xf);			// Inverts xform in place
	transpose(xf);			// Transposes xform in place
	decompose_rot(xf,ang,ax);	// Decompose to angle and axis
	orthogonalize(xf);		// Makes matrix orthonormal
	xf = mix(xf1, xf2, 0.5);	// Interpolation of rigid-body xforms
	bool ok = xf.read("file.xf");	// Read xform from file
	xf.write("file.xf");		// Write xform to file
	xfname("file.ply")		// Returns string("file.xf")
*/

#include "Vec.h"
#include "lineqn.h"
#include "strutil.h"
#include <iomanip>
#include <fstream>


#define inline TRIMESH_INLINE


namespace trimesh {

template <class T>
class XForm {
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

protected:
	// The internal representation: array.  Column-major (OpenGL) order.
	T m[16];

public:
	// Constructors: defaults to identity
	XForm(const T m0 =1, const T m1 =0, const T m2 =0, const T m3 =0,
	      const T m4 =0, const T m5 =1, const T m6 =0, const T m7 =0,
	      const T m8 =0, const T m9 =0, const T m10=1, const T m11=0,
	      const T m12=0, const T m13=0, const T m14=0, const T m15=1)
	{
		m[0]  = m0;  m[1]  = m1;  m[2]  = m2;  m[3]  = m3;
		m[4]  = m4;  m[5]  = m5;  m[6]  = m6;  m[7]  = m7;
		m[8]  = m8;  m[9]  = m9;  m[10] = m10; m[11] = m11;
		m[12] = m12; m[13] = m13; m[14] = m14; m[15] = m15;
	}

	// Constructor for 1 scalar argument, which is duplicated into
	// all components.  Explicit.
	template <class S>
	explicit XForm(S x,
		typename ::std::enable_if< ::std::is_arithmetic<S>::value, void >::type * = 0)
	{
		for (size_type i = 0; i < 16; i++)
			m[i] = x;
	}

	// Constructor for 1 argument that's a pointer, array, or
	// anything else that can be accessed using [].  Explicit.
	template <class S>
	explicit XForm(const S &x)
	{
		for (size_type i = 0; i < 16; i++)
			m[i] = x[i];
	}

	// operator= equivalent of the one-parameter constructor
	template <class S>
	inline XForm &operator = (const S &x)
	{
		for (size_type i = 0; i < 16; i++)
			m[i] = x[i];
		return *this;
	}

	// Default destructor, copy constructor, assignment operator

	// Array reference - no bounds checking
	inline reference operator [] (size_type i)
		{ return m[i]; }
	inline reference operator [] (int i)
		{ return m[i]; }
	inline const_reference operator [] (size_type i) const
		{ return m[i]; }
	inline const_reference operator [] (int i) const
		{ return m[i]; }

	// Access by row/column
	inline reference operator () (size_type r, size_type c)
		{ return m[r + c * 4]; }
	inline reference operator () (int r, int c)
		{ return m[r + c * 4]; }
	inline const_reference operator () (size_type r, size_type c) const
		{ return m[r + c * 4]; }
	inline const_reference operator () (int r, int c) const
		{ return m[r + c * 4]; }

	// Other accessors, for compatibility with std::array
	inline reference front()
		{ return m[0]; }
	inline const_reference front() const
		{ return m[0]; }
	inline reference back()
		{ return m[15]; }
	inline const_reference back() const
		{ return m[15]; }

	// Conversion to pointer
	inline operator T * ()
		{ return m; }
	inline operator const T * ()
		{ return m; }
	inline operator const T * () const
		{ return m; }
	inline pointer data()
		{ return m; }
	inline const_pointer data() const
		{ return m; }

	// Iterators
	inline iterator begin()
		{ return m; }
	inline const_iterator begin() const
		{ return m; }
	inline const_iterator cbegin() const
		{ return m; }
	inline iterator end()
		{ return begin() + 16; }
	inline const_iterator end() const
		{ return begin() + 16; }
	inline const_iterator cend() const
		{ return begin() + 16; }
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
		{ return 16; }
	inline size_type max_size() const
		{ return 16; }

	// Static members - really just fancy constructors
	static inline XForm<T> identity()
		{ return XForm<T>(); }

	static inline XForm<T> trans(const T &tx, const T &ty, const T &tz)
		{ return XForm<T>(1,0,0,0,0,1,0,0,0,0,1,0,tx,ty,tz,1); }

	template <class S>
	static inline XForm<T> trans(const S &t)
		{ return XForm<T>::trans(t[0], t[1], t[2]); }

	static inline XForm<T> rot(const T &angle,
		const T &rx, const T &ry, const T &rz)
	{
		// Angle in radians, unlike OpenGL
		using namespace ::std;
		T l = sqrt(rx*rx + ry*ry + rz*rz);
		if (unlikely(l == 0))
			return XForm<T>();
		T l1 = 1 / l, x = rx * l1, y = ry * l1, z = rz * l1;
		T s = sin(angle), c = cos(angle);
		T xs = x*s, ys = y*s, zs = z*s, c1 = 1 - c;
		T xx = c1*x*x, yy = c1*y*y, zz = c1*z*z;
		T xy = c1*x*y, xz = c1*x*z, yz = c1*y*z;
		return XForm<T>(xx+c,  xy+zs, xz-ys, 0,
				xy-zs, yy+c,  yz+xs, 0,
				xz+ys, yz-xs, zz+c,  0,
				0, 0, 0, 1);
	}

	template <class S>
	static inline XForm<T> rot(const T &angle, const S &axis)
		{ return XForm<T>::rot(angle, axis[0], axis[1], axis[2]); }

	static inline XForm<T> rot_into(T d1x, T d1y, T d1z, T d2x, T d2y, T d2z)
	{
		using namespace ::std;

		// Normalize input directions
		T l1 = sqrt(d1x*d1x + d1y*d1y + d1z*d1z);
		T l2 = sqrt(d2x*d2x + d2y*d2y + d2z*d2z);
		if (unlikely(l1 == 0 || l2 == 0))
			return XForm<T>();
		T rl1 = 1 / l1; d1x *= rl1; d1y *= rl1; d1z *= rl1;
		T rl2 = 1 / l2; d2x *= rl2; d2y *= rl2; d2z *= rl2;

		// Find cosine of rotation angle
		T c = d1x*d2x + d1y*d2y + d1z*d2z;
		if (c <= -1) {
			// Rotation of 180 degrees.  Find an arbitrary
			// rotation axis perpendicular to d1 and d2
			T ax = d1y, ay = -d1x; // az = 0
			if (ax == 0 && ay == 0)
				ax = 1;
			T rla = 1 / sqrt(ax*ax + ay*ay);
			ax *= rla;
			ay *= rla;

			// Construct rotation matrix
			return XForm<T>(2*ax*ax-1, 2*ax*ay,   0, 0,
					2*ax*ay,   2*ay*ay-1, 0, 0,
					0,         0,        -1, 0,
					0, 0, 0, 1);
		}

		// Else construct rotation axis * sin(angle)
		T sx = d1y * d2z - d1z * d2y;
		T sy = d1z * d2x - d1x * d2z;
		T sz = d1x * d2y - d1y * d2x;
		// We need to multiply by (1 - cos) / sin^2 = 1 / (1 + cos)
		T r1c = 1 / (1 + c);
		return XForm<T>(sx*sx*r1c+c,  sx*sy*r1c+sz, sx*sz*r1c-sy, 0,
				sx*sy*r1c-sz, sy*sy*r1c+c,  sy*sz*r1c+sx, 0,
				sx*sz*r1c+sy, sy*sz*r1c-sx, sz*sz*r1c+c,  0,
				0, 0, 0, 1);
	}

	template <class S>
	static inline XForm<T> rot_into(const S &d1, const S &d2)
		{ return XForm<T>::rot_into(d1[0], d1[1], d1[2], d2[0], d2[1], d2[2]); }

	static inline XForm<T> scale(const T &s)
		{ return XForm<T>(s,0,0,0,0,s,0,0,0,0,s,0,0,0,0,1); }

	static inline XForm<T> scale(const T &sx, const T &sy, const T &sz)
		{ return XForm<T>(sx,0,0,0,0,sy,0,0,0,0,sz,0,0,0,0,1); }

	static inline XForm<T> scale(const T &s, const T &dx, const T &dy, const T &dz)
	{
		T dlen2 = dx*dx + dy*dy + dz*dz;
		T s1 = (s - 1) / dlen2;
		return XForm<T>(1 + s1*dx*dx, s1*dx*dy, s1*dx*dz, 0,
				s1*dx*dy, 1 + s1*dy*dy, s1*dy*dz, 0,
				s1*dx*dz, s1*dy*dz, 1 + s1*dz*dz, 0,
				0, 0, 0, 1);
	}

	template <class S>
	static inline XForm<T> scale(const T &s, const S &dir)
		{ return XForm<T>::scale(s, dir[0], dir[1], dir[2]); }

	// OpenGL-like ortho and frustum
	static inline XForm<T> ortho(const T &l, const T &r, const T &b, const T &t,
		const T &n, const T &f)
	{
		T rrl = 1 / (r - l);
		T rtb = 1 / (t - b);
		T rfn = 1 / (f - n);
		return XForm<T>(2*rrl, 0, 0, 0,
				0, 2*rtb, 0, 0,
				0, 0, -2*rfn, 0,
				-(r+l)*rrl, -(t+b)*rtb, -(f+n)*rfn, 1);
	}

	static inline XForm<T> frustum(const T &l, const T &r, const T &b, const T &t,
		const T &n, const T &f)
	{
		T rrl = 1 / (r - l);
		T rtb = 1 / (t - b);
		T rfn = 1 / (f - n);
		return XForm<T>(2*n*rrl, 0, 0, 0,
				0, 2*n*rtb, 0, 0,
				(r+l)*rrl, (t+b)*rtb, -(f+n)*rfn, -1,
				0, 0, -2*f*n*rfn, 0);
	}

	// Returns y*x^T, thinking of y and x as column 3-vectors
	template <class S>
	static inline XForm<T> outer(const S &y, const S &x)
	{
		XForm<T> result;
		for (size_t i = 0; i < 3; i++)
			for (size_t j = 0; j < 3; j++)
				result[4*i+j] = x[i]*y[j];
		return result;
	}

	// Construct from 3x3 or 4x4 array
	template <class S>
	static inline XForm<T> fromarray(const S x[3][3])
	{
		XForm<T> m;
		m[0]  = x[0][0];  m[1]  = x[1][0];
		m[2]  = x[2][0];  m[3]  = 0;
		m[4]  = x[0][1];  m[5]  = x[1][1];
		m[6]  = x[2][1];  m[7]  = 0;
		m[8]  = x[0][2];  m[9]  = x[1][2];
		m[10] = x[2][2];  m[11] = 0;
		m[12] = 0;        m[13] = 0;
		m[14] = 0;        m[15] = 1;
		return m;
	}

	template <class S>
	static inline XForm<T> fromarray(const S x[4][4])
	{
		XForm<T> m;
		m[0]  = x[0][0];  m[1]  = x[1][0];
		m[2]  = x[2][0];  m[3]  = x[3][0];
		m[4]  = x[0][1];  m[5]  = x[1][1];
		m[6]  = x[2][1];  m[7]  = x[3][1];
		m[8]  = x[0][2];  m[9]  = x[1][2];
		m[10] = x[2][2];  m[11] = x[3][2];
		m[12] = x[0][3];  m[13] = x[1][3];
		m[14] = x[2][3];  m[15] = x[3][3];
		return m;
	}

	// Read an XForm from a file.
	inline bool read(const ::std::string &filename)
	{
		using namespace ::std;
		ifstream f(filename.c_str());
		XForm<T> M;
		f >> M;
		f.close();
		if (f.good()) {
			*this = M;
			return true;
		}
		return false;
	}

	// Write an XForm to a file
	inline bool write(const ::std::string &filename) const
	{
		using namespace ::std;
		const int digits = 2 + numeric_limits<T>::digits10;
		ofstream f(filename.c_str());
		f << setprecision(digits) << *this;
		f.close();
		return f.good();
	}
}; // class XForm


typedef XForm<double>  xform;
typedef XForm<double> dxform;
typedef XForm<float>  fxform;


// Binary operations
template <class T>
static inline XForm<T> operator + (const XForm<T> &xf1, const XForm<T> &xf2)
{
	return XForm<T>(
		xf1[ 0] + xf2[ 0], xf1[ 1] + xf2[ 1],
		xf1[ 2] + xf2[ 2], xf1[ 3] + xf2[ 3],
		xf1[ 4] + xf2[ 4], xf1[ 5] + xf2[ 5],
		xf1[ 6] + xf2[ 6], xf1[ 7] + xf2[ 7],
		xf1[ 8] + xf2[ 8], xf1[ 9] + xf2[ 9],
		xf1[10] + xf2[10], xf1[11] + xf2[11],
		xf1[12] + xf2[12], xf1[13] + xf2[13],
		xf1[14] + xf2[14], xf1[15] + xf2[15]
	);
}

template <class T>
static inline XForm<T> operator - (const XForm<T> &xf1, const XForm<T> &xf2)
{
	return XForm<T>(
		xf1[ 0] - xf2[ 0], xf1[ 1] - xf2[ 1],
		xf1[ 2] - xf2[ 2], xf1[ 3] - xf2[ 3],
		xf1[ 4] - xf2[ 4], xf1[ 5] - xf2[ 5],
		xf1[ 6] - xf2[ 6], xf1[ 7] - xf2[ 7],
		xf1[ 8] - xf2[ 8], xf1[ 9] - xf2[ 9],
		xf1[10] - xf2[10], xf1[11] - xf2[11],
		xf1[12] - xf2[12], xf1[13] - xf2[13],
		xf1[14] - xf2[14], xf1[15] - xf2[15]
	);
}

template <class T>
static inline XForm<T> operator * (const XForm<T> &xf1, const XForm<T> &xf2)
{
	return XForm<T>(
		xf1[ 0]*xf2[ 0]+xf1[ 4]*xf2[ 1]+xf1[ 8]*xf2[ 2]+xf1[12]*xf2[ 3],
		xf1[ 1]*xf2[ 0]+xf1[ 5]*xf2[ 1]+xf1[ 9]*xf2[ 2]+xf1[13]*xf2[ 3],
		xf1[ 2]*xf2[ 0]+xf1[ 6]*xf2[ 1]+xf1[10]*xf2[ 2]+xf1[14]*xf2[ 3],
		xf1[ 3]*xf2[ 0]+xf1[ 7]*xf2[ 1]+xf1[11]*xf2[ 2]+xf1[15]*xf2[ 3],
		xf1[ 0]*xf2[ 4]+xf1[ 4]*xf2[ 5]+xf1[ 8]*xf2[ 6]+xf1[12]*xf2[ 7],
		xf1[ 1]*xf2[ 4]+xf1[ 5]*xf2[ 5]+xf1[ 9]*xf2[ 6]+xf1[13]*xf2[ 7],
		xf1[ 2]*xf2[ 4]+xf1[ 6]*xf2[ 5]+xf1[10]*xf2[ 6]+xf1[14]*xf2[ 7],
		xf1[ 3]*xf2[ 4]+xf1[ 7]*xf2[ 5]+xf1[11]*xf2[ 6]+xf1[15]*xf2[ 7],
		xf1[ 0]*xf2[ 8]+xf1[ 4]*xf2[ 9]+xf1[ 8]*xf2[10]+xf1[12]*xf2[11],
		xf1[ 1]*xf2[ 8]+xf1[ 5]*xf2[ 9]+xf1[ 9]*xf2[10]+xf1[13]*xf2[11],
		xf1[ 2]*xf2[ 8]+xf1[ 6]*xf2[ 9]+xf1[10]*xf2[10]+xf1[14]*xf2[11],
		xf1[ 3]*xf2[ 8]+xf1[ 7]*xf2[ 9]+xf1[11]*xf2[10]+xf1[15]*xf2[11],
		xf1[ 0]*xf2[12]+xf1[ 4]*xf2[13]+xf1[ 8]*xf2[14]+xf1[12]*xf2[15],
		xf1[ 1]*xf2[12]+xf1[ 5]*xf2[13]+xf1[ 9]*xf2[14]+xf1[13]*xf2[15],
		xf1[ 2]*xf2[12]+xf1[ 6]*xf2[13]+xf1[10]*xf2[14]+xf1[14]*xf2[15],
		xf1[ 3]*xf2[12]+xf1[ 7]*xf2[13]+xf1[11]*xf2[14]+xf1[15]*xf2[15]
	);
}


// Matrix-vector multiplication
template <class S, class T>
static inline const S operator * (const XForm<T> &xf, const S &v)
{
	T v0 = v[0], v1 = v[1], v2 = v[2];
	T h = 1 / (xf[3] * v0 + xf[7] * v1 + xf[11] * v2 + xf[15]);

	typedef typename S::value_type Stype;
	return S(Stype(h*(xf[0] * v0 + xf[4] * v1 + xf[8] * v2 + xf[12])),
	         Stype(h*(xf[1] * v0 + xf[5] * v1 + xf[9] * v2 + xf[13])),
	         Stype(h*(xf[2] * v0 + xf[6] * v1 + xf[10] * v2 + xf[14])));
}


// Component-wise equality and inequality (#include the usual caveats
// about comparing floats for equality...)
template <class T>
static inline bool operator == (const XForm<T> &xf1, const XForm<T> &xf2)
{
	for (size_t i = 0; i < 16; i++)
		if (xf1[i] != xf2[i])
			return false;
	return true;
}

template <class T>
static inline bool operator != (const XForm<T> &xf1, const XForm<T> &xf2)
{
	for (size_t i = 0; i < 16; i++)
		if (xf1[i] != xf2[i])
			return true;
	return false;
}


// Inverse
template <class T>
static inline XForm<T> inv(const XForm<T> &xf)
{
	T A[4][4] = { { xf[0], xf[4], xf[8],  xf[12] },
	              { xf[1], xf[5], xf[9],  xf[13] },
	              { xf[2], xf[6], xf[10], xf[14] },
	              { xf[3], xf[7], xf[11], xf[15] } };
	int ind[4];
	bool ok = ludcmp<T,4>(A, ind);
	if (unlikely(!ok))
		return XForm<T>();
	T B[4][4] = { { 1, 0, 0, 0 },
	              { 0, 1, 0, 0 },
	              { 0, 0, 1, 0 },
	              { 0, 0, 0, 1 } };
	for (size_t i = 0; i < 4; i++)
		lubksb<T,4>(A, ind, B[i]);
	return XForm<T>(B[0][0], B[0][1], B[0][2], B[0][3],
			B[1][0], B[1][1], B[1][2], B[1][3],
			B[2][0], B[2][1], B[2][2], B[2][3],
			B[3][0], B[3][1], B[3][2], B[3][3]);
}

template <class T>
static inline void invert(XForm<T> &xf)
{
	xf = inv(xf);
}


// Transpose
template <class T>
static inline XForm<T> transp(const XForm<T> &xf)
{
	return XForm<T>(xf[0], xf[4], xf[8], xf[12],
			xf[1], xf[5], xf[9], xf[13],
			xf[2], xf[6], xf[10], xf[14],
			xf[3], xf[7], xf[11], xf[15]);
}

template <class T>
static inline void transpose(XForm<T> &xf)
{
	xf = transp(xf);
}


// Just the rotational or translational parts
template <class T>
static inline XForm<T> rot_only(const XForm<T> &xf)
{
	return XForm<T>(xf[0], xf[1], xf[2], 0,
			xf[4], xf[5], xf[6], 0,
			xf[8], xf[9], xf[10], 0,
			0, 0, 0, 1);
}

template <class T>
static inline XForm<T> trans_only(const XForm<T> &xf)
{
	return XForm<T>(1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			xf[12], xf[13], xf[14], 1);
}


// Inverse transpose - the transformation to be applied to normals
template <class T>
static inline XForm<T> norm_xf(const XForm<T> &xf)
{
	using namespace ::std;
	XForm<T> M = inv(xf);
	M[12] = M[13] = M[14] = 0;
	swap(M[1], M[4]);
	swap(M[2], M[8]);
	swap(M[6], M[9]);
	return M;
}


// Decompose the rotational part of an xform into angle and axis
template <class T, class S>
static inline void decompose_rot(const XForm<T> &xf, T &angle, S &axis)
{
	using namespace ::std;

	// Initial estimate from antisymmetric component
	T ax = xf[6] - xf[9];
	T ay = xf[8] - xf[2];
	T az = xf[1] - xf[4];
	T l = sqrt(ax*ax + ay*ay + az*az);
	T tr1 = xf[0] + xf[5] + xf[10] - 1.0;
	angle = atan2(l, tr1);

	if (angle < T(M_PI_2l)) {
		axis[0] = ax;
		axis[1] = ay;
		axis[2] = az;
	} else {
		// If rotation is over pi/2, using the symmetric component
		// to find the rotation axis is more accurate than
		// using the antisymmetric part.
		// But we still need the latter to disambiguate flip.
		T tx = abs(2 * xf[0] - tr1);
		T ty = abs(2 * xf[5] - tr1);
		T tz = abs(2 * xf[10] - tr1);
		if (tx > ty && tx > tz) {
			axis[0] = 2 * xf[0] - tr1;
			axis[1] = xf[1] + xf[4];
			axis[2] = xf[8] + xf[2];
			if (xf[6] - xf[9] < 0) {
				axis[0] = -axis[0];
				axis[1] = -axis[1];
				axis[2] = -axis[2];
			}
		} else if (ty > tx && ty > tz) {
			axis[0] = xf[1] + xf[4];
			axis[1] = 2 * xf[5] - tr1;
			axis[2] = xf[6] + xf[9];
			if (xf[8] - xf[2] < 0) {
				axis[0] = -axis[0];
				axis[1] = -axis[1];
				axis[2] = -axis[2];
			}
		} else {
			axis[0] = xf[8] + xf[2];
			axis[1] = xf[6] + xf[9];
			axis[2] = 2 * xf[10] - tr1;
			if (xf[1] - xf[4] < 0) {
				axis[0] = -axis[0];
				axis[1] = -axis[1];
				axis[2] = -axis[2];
			}
		}
	}

	T len_ax = sqrt(sqr(axis[0]) + sqr(axis[1]) + sqr(axis[2]));
	if (likely(len_ax != 0)) {
		T scale = 1 / len_ax;
		axis[0] *= scale;
		axis[1] *= scale;
		axis[2] *= scale;
	} else {
		axis[2] = 1;
	}
}


// Make orthonormal
template <class T>
static inline void orthogonalize(XForm<T> &xf)
{
	if (unlikely(xf[15] == 0)) {	// Yuck.  Doesn't make sense...
		xf[15] = 1;
	} else if (unlikely(xf[15] != 1)) {
		T scale = 1 / xf[15];
		for (size_t i = 0; i < 15; i++)
			xf[i] *= scale;
		xf[15] = 1;
	}

	T angle;
	Vec<3,T> axis;
	decompose_rot(xf, angle, axis);
	Vec<3,T> trans(xf[12], xf[13], xf[14]);

	xf = XForm<T>::rot(angle, axis);
	xf[12] = trans[0];
	xf[13] = trans[1];
	xf[14] = trans[2];
}


// Interpolate between xforms, returning an xform that's a fraction "a" of
// the way from x to y.  Performs screw decomposition, interpolates rotation
// and translation along screw axis, and constructs in-plane translation.
template <class T, class S>
static inline XForm<T> mix(const XForm<T> &x, const XForm<T> &y, const S &a)
{
	XForm<T> x2y = y * inv(x);

	// Find angle/axis of rotation
	T angle;
	Vec<3,T> axis;
	decompose_rot(x2y, angle, axis);

	// Split the translation into two components:
	// along axis, and in the plane perpendicular to axis
	Vec<3,T> trans(x2y[12], x2y[13], x2y[14]);
	Vec<3,T> along = (trans DOT axis) * axis;
	Vec<3,T> in_plane = trans - along;
	Vec<3,T> perp = in_plane CROSS axis;

	// Now build the interpolated transform
	T theta = T(0.5) * angle, a1 = 1 - a;
	T factor = (theta > T(1.0e-8)) ? sin(a * theta) / sin(theta) : a;
	Vec<3,T> trans_out = a * along +
		factor * cos(a1 * theta) * in_plane +
		factor * sin(a1 * theta) * perp;

	return XForm<T>::trans(trans_out) *
	       XForm<T>::rot(a * angle, axis) * x;
}


// iostream operators
template <class T>
static inline ::std::ostream &operator << (::std::ostream &os, const XForm<T> &m)
{
	using namespace ::std;
	for (size_t i = 0; i < 4; i++)
		for (size_t j = 0; j < 4; j++)
			os << m[i+4*j] << (j == 3 ? "\n" : " ");
	return os;
}
template <class T>
static inline ::std::istream &operator >> (::std::istream &is, XForm<T> &m)
{
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 4; j++)
			is >> m[i+4*j];
	}

	if (!is.good()) {
		m = xform::identity();
		return is;
	}

	// 4th row is allowed to fail
	for (size_t j = 0; j < 4; j++)
		is >> m[3+4*j];

	if (!is.good()) {
		m[3] = m[7] = m[11] = 0;
		m[15] = 1;
		is.clear();
	}

	return is;
}


// Generate a .xf filename from an input (scan) filename
static inline ::std::string xfname(const ::std::string &filename)
{
	return replace_ext(filename, "xf");
}


} // namespace trimesh

#undef inline

#endif
