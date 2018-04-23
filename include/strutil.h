#ifndef STRUTIL_H
#define STRUTIL_H
/*
Szymon Rusinkiewicz
Princeton University

strutil.h
Miscellaneous string-manipulation utilities

Usage:
	std::string s("foo.bar");
	std::string s2 = replace_ext(s, "baz");  // "foo.baz"
	std::string s2 = replace_ext(s, "");     // "foo"
	std::string s3 = string_printf("%d", 42);

	// All of the below are case-insensitive
	begins_with("Foobar", "foo")             // true
	ends_with("foobar", "baz")               // false
	contains("fooBar", "bar")                // true
*/


#include <string>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <cstdarg>

#ifdef _MSC_VER
# ifndef strcasecmp
#  define strcasecmp _stricmp
# endif
# ifndef strncasecmp
#  define strncasecmp _strnicmp
# endif
# ifndef strdup
#  define strdup _strdup
# endif
#endif


namespace trimesh {

// Replace the extension of a filename, else add one if none present
static inline ::std::string replace_ext(const ::std::string &filename,
	const ::std::string &ext)
{
	using namespace ::std;
	string x = filename;
	string::size_type dot = x.rfind(".", x.length());
	if (dot != string::npos)
		x.erase(dot);
	if (ext.empty())
		return x;
	else
		return x + string(".") + ext;
}


// sprintf to a c++ string
static inline ::std::string string_printf(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	size_t required = 1 + ::std::vsnprintf(NULL, 0, format, ap);
	va_end(ap);
	char *buf = new char[required];
	va_start(ap, format);
	::std::vsnprintf(buf, required, format, ap);
	va_end(ap);
	::std::string s(buf);
	delete [] buf;
	return s;
}


// Does string s1 begin with / end with / contain s2?  (Case-insensitive)
static inline bool begins_with(const char *s1, const char *s2)
{
	using namespace ::std;
	return !strncasecmp(s1, s2, strlen(s2));
}

static inline bool begins_with(const ::std::string &s1, const ::std::string &s2)
{
	return begins_with(s1.c_str(), s2.c_str());
}

static inline bool ends_with(const char *s1, const char *s2)
{
	using namespace ::std;
	size_t l1 = strlen(s1), l2 = strlen(s2);
	return (l1 >= l2) && !strncasecmp(s1 + l1 - l2, s2, l2);
}

static inline bool ends_with(const ::std::string &s1, const ::std::string &s2)
{
	return ends_with(s1.c_str(), s2.c_str());
}

static inline bool contains(const char *s1, const char *s2)
{
	using namespace ::std;
	size_t len1 = strlen(s1), len2 = strlen(s2);

	// strcasestr is not portable, so simulate it
	char *lower1 = new char[len1 + 1];
	for (size_t i = 0; i < len1; i++)
		lower1[i] = tolower(s1[i]);
	lower1[len1] = '\0';

	char *lower2 = new char[len2 + 1];
	for (size_t i = 0; i < len2; i++)
		lower2[i] = tolower(s2[i]);
	lower2[len2] = '\0';

	bool found = (bool) strstr(lower1, lower2);

	delete [] lower2;
	delete [] lower1;
	return found;
}

static inline bool contains(const ::std::string &s1, const ::std::string &s2)
{
	return contains(s1.c_str(), s2.c_str());
}

} // namespace trimesh

#endif
