#ifndef ENDIANUTIL_H
#define ENDIANUTIL_H
/*
Szymon Rusinkiewicz
Princeton University

endianutil.h
Endian-ness and byte-swapping
*/

#if defined(__linux__)
# include <endian.h>
# include <byteswap.h>
#elif defined(__APPLE__)
# include <libkern/OSByteOrder.h>
#elif defined(_MSC_VER)
# include <cstdlib>
#endif


namespace trimesh {


// Figure out whether this machine is little- or big-endian
static inline bool we_are_little_endian()
{
#if defined(__linux__)
	return (__BYTE_ORDER == __LITTLE_ENDIAN);
#elif defined(__APPLE__)
	return (OSHostByteOrder() == OSLittleEndian);
#elif defined(_MSC_VER)
	return true;
#else
	// The following appears to be legal according to
	// C99 strict-aliasing rules
	static const int tmp = 1;
	static const bool little_endian = !!(* (const unsigned char *) &tmp);
	return little_endian;
#endif
}

static inline bool we_are_big_endian()
{
	return !we_are_little_endian();
}


// Byte-swap 16-, 32-, and 64-bit quantities
static inline void swap_16(unsigned char *p)
{
#if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 408)
	* (unsigned short *) p = __builtin_bswap16(* (unsigned short *) p);
#elif defined(__linux__)
	* (unsigned short *) p = bswap_16(* (unsigned short *) p);
#elif defined(__APPLE__)
	* (unsigned short *) p = OSSwapInt16(* (unsigned short *) p);
#elif defined(_MSC_VER) && _MSC_VER >= 1020
	* (unsigned short *) p = _byteswap_ushort(* (unsigned short *) p);
#else
	* (unsigned short *) p =
		(((* (unsigned short *) p) & 0xff00u) >> 8) |
		(((* (unsigned short *) p) & 0x00ffu) << 8);
#endif
}

static inline void swap_32(unsigned char *p)
{
#if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 403)
	* (unsigned *) p = __builtin_bswap32(* (unsigned *) p);
#elif defined(__linux__)
	* (unsigned *) p = bswap_32(* (unsigned *) p);
#elif defined(__APPLE__)
	* (unsigned *) p = OSSwapInt32(* (unsigned *) p);
#elif defined(_MSC_VER) && _MSC_VER >= 1020
	* (unsigned *) p = _byteswap_ulong(* (unsigned *) p);
#else
	* (unsigned *) p =
		(((* (unsigned *) p) & 0xff000000u) >> 24) |
		(((* (unsigned *) p) & 0x00ff0000u) >>  8) |
		(((* (unsigned *) p) & 0x0000ff00u) <<  8) |
		(((* (unsigned *) p) & 0x000000ffu) << 24);
#endif
}

static inline void swap_64(unsigned char *p)
{
#if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 403)
	* (unsigned long long *) p =
		__builtin_bswap64(* (unsigned long long *) p);
#elif defined(__linux__)
	* (unsigned long long *) p = bswap_64(* (unsigned long long *) p);
#elif defined(__APPLE__)
	* (unsigned long long *) p = OSSwapInt64(* (unsigned long long *) p);
#elif defined(_MSC_VER) && _MSC_VER >= 1020
	* (unsigned __int64 *) p = _byteswap_uint64(* (unsigned __int64 *) p);
#else
	* (unsigned long long *) p =
		(((* (unsigned long long *) p) & 0xff00000000000000ull) >> 56) |
		(((* (unsigned long long *) p) & 0x00ff000000000000ull) >> 40) |
		(((* (unsigned long long *) p) & 0x0000ff0000000000ull) >> 24) |
		(((* (unsigned long long *) p) & 0x000000ff00000000ull) >>  8) |
		(((* (unsigned long long *) p) & 0x00000000ff000000ull) <<  8) |
		(((* (unsigned long long *) p) & 0x0000000000ff0000ull) << 24) |
		(((* (unsigned long long *) p) & 0x000000000000ff00ull) << 40) |
		(((* (unsigned long long *) p) & 0x00000000000000ffull) << 56);
#endif
}


// Byte swap ints, uints, and floats.  Assumes int is 32bits.
// Going through (unsigned char *) appears to be the legal
// (C99 strict-aliasing compliant) way to do this.
static inline void swap_short(short &x)
{
	swap_16((unsigned char *) &x);
}

static inline void swap_ushort(unsigned short &x)
{
	swap_16((unsigned char *) &x);
}

static inline void swap_int(int &x)
{
	swap_32((unsigned char *) &x);
}

static inline void swap_unsigned(unsigned &x)
{
	swap_32((unsigned char *) &x);
}

static inline void swap_float(float &x)
{
	swap_32((unsigned char *) &x);
}

static inline void swap_double(double &x)
{
	swap_64((unsigned char *) &x);
}


} // namespace trimesh

#endif
