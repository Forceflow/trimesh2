/*
Szymon Rusinkiewicz
Princeton University

xf.cc
Simple transformations on .xf files
*/

#include "Vec.h"
#include "XForm.h"
#include "timestamp.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;
using namespace trimesh;


// Is this argument a floating-point number?
static bool isanumber(const char *c)
{
	if (!c || !*c)
		return false;
	char *endptr;
	strtod(c, &endptr);
	return (endptr && *endptr == '\0');
}


// Randomizes random number generator based on timestamp
static void randomize()
{
	union {
		timestamp t;
		struct { unsigned u1, u2, u3, u4; } u;
	} t;
	t.t = now();
	unsigned seed = t.u.u1 ^ t.u.u2 ^ t.u.u3 ^ t.u.u4;
	xorshift_rnd(seed);
	xorshift_rnd();
	xorshift_rnd();
}


// Returns a random unit-length vector on the sphere
static vec rnd_vec()
{
	float phi = uniform_rnd(M_2PIf);
	float z = uniform_rnd(-1.0f);
	float xy = sqrt(1.0f - sqr(z));
	return vec(xy * cos(phi), xy * sin(phi), z);
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s commands...\n", myname);
	fprintf(stderr, "Commands:\n");
	fprintf(stderr, "	-xform file.xf	Transform by the given matrix\n");
	fprintf(stderr, "	file.xf		Transform by the given matrix\n");
	fprintf(stderr, "	-ixform file.xf	Transform by inverse of matrix\n");
	fprintf(stderr, "	inv(file.xf)	Transform by inverse of matrix\n");
	fprintf(stderr, "	-rot r x y z	Rotate r degrees around axis (x,y,z)\n");
	fprintf(stderr, "	-q qr qi qj qk	Rotate by given quaternion\n");
	fprintf(stderr, "	-v qi qj qk iqr	Rotate by given VRIP-style quaternion\n");
	fprintf(stderr, "	-trans x y z	Translate by (x,y,z)\n");
	fprintf(stderr, "	-rnd rot trans	Apply a random transform with maximum magnitude in rot, trans\n");
	fprintf(stderr, "	-inv		Invert the current transformation\n");
	fprintf(stderr, "	-print		Display current transformation on stdout\n");
	fprintf(stderr, "	-o file.xf	Write current transformation to file\n");
	fprintf(stderr, "	-prot		Display current transformation as angle/axis\n");
	fprintf(stderr, "	-pq		Display current transformation as quaternion\n");
	fprintf(stderr, "	-pv		Display current transformation as VRIP-style quaternion\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 2)
		usage(argv[0]);

	xform xf;
	bool autoprint = true;

	for (int i = 1; i < argc; i++) {
		autoprint = true;
		if (!strncmp(argv[i], "inv(", 4) &&
		    argv[i][strlen(argv[i])-1] == ')') {
			int l = strlen(argv[i]);
			char *filename = new char[l-4];
			strncpy(filename, argv[i]+4, l-5);
			filename[l-5] = '\0';
			xform tmp;
			if (!tmp.read(filename)) {
				fprintf(stderr, "Couldn't read %s\n", filename);
				usage(argv[0]);
			}
			xf = xf * inv(tmp);
			delete [] filename;
		} else if (argv[i][0] != '-') {
			xform tmp;
			if (!tmp.read(argv[i])) {
				fprintf(stderr, "Couldn't read %s\n", argv[i]);
				usage(argv[0]);
			}
			xf = xf * tmp;
		} else if (!strcmp(argv[i], "-xf") ||
		           !strcmp(argv[i], "-xform")) {
			i++;
			if (!(i < argc)) {
				fprintf(stderr, "\n-xform requires one argument\n\n");
				usage(argv[0]);
			}
			xform tmp;
			if (!tmp.read(argv[i])) {
				fprintf(stderr, "Couldn't read %s\n", argv[i]);
				usage(argv[0]);
			}
			xf = xf * tmp;
		} else if (!strcmp(argv[i], "-ixf") ||
		           !strcmp(argv[i], "-ixform")) {
			i++;
			if (!(i < argc)) {
				fprintf(stderr, "\n-ixform requires one argument\n\n");
				usage(argv[0]);
			}
			xform tmp;
			if (!tmp.read(argv[i])) {
				fprintf(stderr, "Couldn't read %s\n", argv[i]);
				usage(argv[0]);
			}
			xf = xf * inv(tmp);
		} else if (!strcmp(argv[i], "-rot") ||
		           !strcmp(argv[i], "-rotate")) {
			i += 4;
			if (!(i < argc &&
			      isanumber(argv[i]) && isanumber(argv[i-1]) &&
			      isanumber(argv[i-2]) && isanumber(argv[i-3]))) {
				fprintf(stderr, "\n-rot requires four arguments\n\n");
				usage(argv[0]);
			}
			Vec<3,double> ax(atof(argv[i-2]), atof(argv[i-1]), atof(argv[i]));
			double ang = radians(atof(argv[i-3]));
			xf = xf * xform::rot(ang, ax);
		} else if (!strcmp(argv[i], "-q") ||
		           !strcmp(argv[i], "-quat")) {
			i += 4;
			if (!(i < argc &&
			      isanumber(argv[i]) && isanumber(argv[i-1]) &&
			      isanumber(argv[i-2]) && isanumber(argv[i-3]))) {
				fprintf(stderr, "\n-q requires four arguments\n\n");
				usage(argv[0]);
			}
			Vec<4,double> q(atof(argv[i-3]), atof(argv[i-2]), atof(argv[i-1]), atof(argv[i]));
			normalize(q);
			double c2 = q[0];
			double s2 = sqrt(sqr(q[1])+sqr(q[2])+sqr(q[3]));
			double ang = 2.0 * atan2(s2, c2);
			xf = xf * xform::rot(ang, q[1], q[2], q[3]);
		} else if (!strcmp(argv[i], "-v") ||
		           !strcmp(argv[i], "-vq") ||
		           !strcmp(argv[i], "-vquat") ||
		           !strcmp(argv[i], "-vripquat")) {
			i += 4;
			if (!(i < argc &&
			      isanumber(argv[i]) && isanumber(argv[i-1]) &&
			      isanumber(argv[i-2]) && isanumber(argv[i-3]))) {
				fprintf(stderr, "\n-v requires four arguments\n\n");
				usage(argv[0]);
			}
			Vec<4,double> q(-atof(argv[i]), atof(argv[i-3]), atof(argv[i-2]), atof(argv[i-1]));
			normalize(q);
			double c2 = q[0];
			double s2 = sqrt(sqr(q[1])+sqr(q[2])+sqr(q[3]));
			double ang = 2.0 * atan2(s2, c2);
			xf = xf * xform::rot(ang, q[1], q[2], q[3]);
		} else if (!strcmp(argv[i], "-trans") ||
		           !strcmp(argv[i], "-translate")) {
			i += 3;
			if (!(i < argc && isanumber(argv[i]) &&
			      isanumber(argv[i-1]) && isanumber(argv[i-2]))) {
				fprintf(stderr, "\n-trans requires three arguments\n\n");
				usage(argv[0]);
			}
			Vec<3,double> t(atof(argv[i-2]), atof(argv[i-1]), atof(argv[i]));
			xf = xf * xform::trans(t);
		} else if (!strcmp(argv[i], "-rnd") ||
		           !strcmp(argv[i], "-random")) {
			i += 2;
			if (!(i < argc && isanumber(argv[i]) &&
			      isanumber(argv[i-1]))) {
				fprintf(stderr, "\n-rnd requires two arguments\n\n");
				usage(argv[0]);
			}
			randomize();
			float rotamount = uniform_rnd(atof(argv[i-1]));
			vec rotaxis = rnd_vec();
			float transamount = uniform_rnd(atof(argv[i]));
			vec transdir = rnd_vec();
			xf = xf * xform::trans(transamount * transdir) *
				xform::rot(rotamount, rotaxis);
		} else if (!strcmp(argv[i], "-inv") ||
		           !strcmp(argv[i], "-invert")) {
			invert(xf);
		} else if (!strcmp(argv[i], "-print")) {
			std::cout << xf;
			autoprint = false;
		} else if (!strcmp(argv[i], "-o")) {
			i++;
			if (!(i < argc)) {
				fprintf(stderr, "\n-o requires one argument\n\n");
				usage(argv[0]);
			}
			xf.write(argv[i]);
			autoprint = false;
		} else if (!strcmp(argv[i], "-prot")) {
			double angle;
			dvec3 axis;
			decompose_rot(xf, angle, axis);
			std::cout << "Rotation of " << angle <<
				" (" << degrees(angle) << " degrees) about " <<
				axis << std::endl;
			dvec3 trans(xf[12], xf[13], xf[14]);
			std::cout << "Translation of " << trans << std::endl;
			double along = trans DOT axis;
			double perp = sqrt(len2(trans) - sqr(along));
			std::cout << "  " << along << " along axis and " <<
				perp << " perpendicular" << std::endl;
			autoprint = false;
		} else if (!strcmp(argv[i], "-pq")) {
			double angle;
			dvec3 axis;
			decompose_rot(xf, angle, axis);
			axis *= sin(0.5 * angle);
			std::cout << cos(0.5 * angle) << " " <<
				axis[0] << " " << axis[1] << " " <<
				axis[2] << "  " << xf[12] << " " <<
				xf[13] << " " << xf[14] << std::endl;
			autoprint = false;
		} else if (!strcmp(argv[i], "-pv")) {
			double angle;
			dvec3 axis;
			decompose_rot(xf, angle, axis);
			axis *= sin(-0.5 * angle);
			std::cout << "bmesh " << xf[12] << " " << xf[13] <<
				" " << xf[14] <<  " " << axis[0] << " " <<
				axis[1] << " " << axis[2] << " " <<
				cos(0.5 * angle) << std::endl;
			autoprint = false;
		}
	}

	if (autoprint)
		std::cout << xf;
}
