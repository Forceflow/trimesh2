/*
Szymon Rusinkiewicz
Princeton University

xf.cc
Simple transformations on .xf files
*/

#include "Vec.h"
#include "XForm.h"
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
			double ang = M_PI / 180.0 * atof(argv[i-3]);
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
			Vec<3,double> rotaxis(xf[6] - xf[9], xf[8] - xf[2], xf[1] - xf[4]);
			double l = len(rotaxis);
			if (l)
				rotaxis /= l;
			double tr1 = xf[0] + xf[5] + xf[10] - 1.0;
			double rotamount = atan2(l, tr1);
			double drotamount = rotamount * 180.0 / M_PI;
			std::cout << "Rotation of " << rotamount << 
				" (" << drotamount << " degrees) about " <<
				rotaxis << std::endl;
			Vec<3,double> trans(xf[12], xf[13], xf[14]);
			std::cout << "Translation of " << trans << std::endl;
			double along = trans DOT rotaxis;
			double perp = sqrt(len2(trans) - sqr(along));
			std::cout << "  " << along << " along axis and " <<
				perp << " perpendicular" << std::endl;
			Vec<3,double> center;
			if (l) {
				double cothalf = (2.0 + tr1) / l;
				Vec<3,double> cr = rotaxis CROSS trans;
				center = 0.5 * (cothalf * cr - rotaxis CROSS cr);
			}
			std::cout << "  center = " << center << std::endl;
			autoprint = false;
		} else if (!strcmp(argv[i], "-pq")) {
			Vec<3,double> rotaxis(xf[6] - xf[9], xf[8] - xf[2], xf[1] - xf[4]);
			double l = len(rotaxis);
			if (l)
				rotaxis /= l;
			double rotamount = atan2(l, xf[0] + xf[5] + xf[10] - 1.0);
			rotaxis *= sin(0.5 * rotamount);
			std::cout << cos(0.5 * rotamount) << " " <<
				rotaxis[0] << " " << rotaxis[1] << " " <<
				rotaxis[2] << "  " << xf[12] << " " <<
				xf[13] << " " << xf[14] << std::endl;
			autoprint = false;
		} else if (!strcmp(argv[i], "-pv")) {
			Vec<3,double> rotaxis(xf[6] - xf[9], xf[8] - xf[2], xf[1] - xf[4]);
			double l = len(rotaxis);
			if (l)
				rotaxis /= l;
			double rotamount = atan2(l, xf[0] + xf[5] + xf[10] - 1.0);
			rotaxis *= -sin(0.5 * rotamount);
			std::cout << "bmesh " << xf[12] << " " << xf[13] <<
				" " << xf[14] <<  " " << rotaxis[0] << " " <<
				rotaxis[1] << " " << rotaxis[2] << " " <<
				cos(0.5 * rotamount) << std::endl;
			autoprint = false;
		}
	}

	if (autoprint)
		std::cout << xf;
}

