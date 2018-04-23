/*
Szymon Rusinkiewicz
Princeton University

mesh_make.cc
Create various kinds of meshes for testing...
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;
using namespace trimesh;

#define ATOF(x) ((float) atof(x))

const char *fixed_shape_names[] = {
	"tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron",
	"truncated_tetrahedron", "cuboctahedron", "truncated_cube",
	"truncated_octahedron", "rhombicuboctahedron",
	"truncated_cuboctahedron", "icosidodecahedron",
	"truncated_dodecahedron", "truncated_icosahedron", "snub_cube",
	"rhombicosidodecahedron", "truncated_icosidodecahedron",
	"snub_dodecahedron", "triakis_tetrahedron", "rhombic_dodecahedron",
	"triakis_octahedron", "tetrakis_hexahedron",
	"deltoidal_icositetrahedron", "disdyakis_dodecahedron",
	"rhombic_triacontahedron", "triakis_icosahedron",
	"pentakis_dodecahedron", "pentagonal_icositetrahedron",
	"deltoidal_hexecontahedron", "disdyakis_triacontahedron",
	"pentagonal_hexecontahedron",
};


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s shape outfile\n", myname);
	fprintf(stderr, "Parameterized shapes:\n");
	fprintf(stderr, "	plane m [n]	m x n tessellated square (default n = m)\n");
	fprintf(stderr, "	bump n sigma	n x n tessellated Gaussian bump of width sigma\n");
	fprintf(stderr, "	wave n omega	n x n tessellated sine wave of frequency omega\n");
	fprintf(stderr, "	frac n		n x n fractal landscape\n");
	fprintf(stderr, "	cube n		n x n tessellated cube\n");
	fprintf(stderr, "	disc n m	Circular disc, tessellated with m rings of n points\n");
	fprintf(stderr, "	cyl n m [r]	Cylinder of radius r (default 1)\n");
	fprintf(stderr, "	ccyl n m [r]	Disc-capped cylinder\n");
	fprintf(stderr, "	scyl n m [r]	Hemisphere-capped cylinder\n");
	fprintf(stderr, "	cone n m [r]	Cone\n");
	fprintf(stderr, "	ccone n m [r]	Capped cone\n");
	fprintf(stderr, "	torus n m [r]	Torus of minor radius r (default 0.25)\n");
	fprintf(stderr, "	knot n m [r]	Trefoil knot of minor radius r (default 0.2)\n");
	fprintf(stderr, "	klein n m	Klein bottle\n");
	fprintf(stderr, "	helix n m t [r]	Helix of minor radius r, with t turns\n");
	fprintf(stderr, "	sphere n m	Sphere, tessellated in polar coordinates\n");
	fprintf(stderr, "	platonic n	Platonic solid with n sides\n");
	fprintf(stderr, "	ssphere n m	Sphere, subdivided m times from a Platonic of n sides\n");
	fprintf(stderr, "	sor n curve.txt	Surface of revolution: n copies of curve, rot z axis\n");
	fprintf(stderr, "	teapot n [b]	Teapot - pass b=1 to omit bottom\n");
	fprintf(stderr, "	orgteapot n [b]	Original (taller) teapot\n");
	fprintf(stderr, "\nFixed polyhedra:\n\t");

	int n = sizeof(fixed_shape_names) / sizeof(fixed_shape_names[0]);
	int chars_printed = 0;
	for (int i = 0; i < n; i++) {
		chars_printed += strlen(fixed_shape_names[i]) + 1;
		if (chars_printed > 72) {
			fprintf(stderr, "\n\t");
			chars_printed = strlen(fixed_shape_names[i]) + 1;
		}
		fprintf(stderr, "%s ", fixed_shape_names[i]);
	}
	fprintf(stderr, "\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 3)
		usage(argv[0]);

	TriMesh *mesh = NULL;
	const char *outfilename = argv[2];

	if (!strcmp(argv[1], "plane") && argc > 4) {
		mesh = make_plane(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "plane") && argc > 3) {
		mesh = make_plane(atoi(argv[2]), atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "bump") && argc > 4) {
		mesh = make_bump(atoi(argv[2]), ATOF(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "wave") && argc > 4) {
		mesh = make_wave(atoi(argv[2]), ATOF(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "frac") && argc > 3) {
		mesh = make_frac(atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "cube") && argc > 3) {
		mesh = make_cube(atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "disc") && argc > 4) {
		mesh = make_disc(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "disk") && argc > 4) {
		mesh = make_disc(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "cyl") && argc > 5) {
		mesh = make_cyl(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "cyl") && argc > 4) {
		mesh = make_cyl(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "ccyl") && argc > 5) {
		mesh = make_ccyl(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "ccyl") && argc > 4) {
		mesh = make_ccyl(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "scyl") && argc > 5) {
		mesh = make_scyl(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "scyl") && argc > 4) {
		mesh = make_scyl(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "cone") && argc > 5) {
		mesh = make_cone(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "cone") && argc > 4) {
		mesh = make_cone(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "ccone") && argc > 5) {
		mesh = make_ccone(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "ccone") && argc > 4) {
		mesh = make_ccone(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "torus") && argc > 5) {
		mesh = make_torus(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "torus") && argc > 4) {
		mesh = make_torus(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "knot") && argc > 5) {
		mesh = make_knot(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "knot") && argc > 4) {
		mesh = make_knot(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "klein") && argc > 4) {
		mesh = make_klein(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "helix") && argc > 6) {
		mesh = make_helix(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]), ATOF(argv[5]));
		outfilename = argv[6];
	} else if (!strcmp(argv[1], "helix") && argc > 5) {
		mesh = make_helix(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "sphere") && argc > 4) {
		mesh = make_sphere_polar(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "platonic") && argc > 3) {
		mesh = make_platonic(atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "ssphere") && argc > 4) {
		mesh = make_sphere_subdiv(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "rd") && argc > 2) {
		mesh = make_fixed_shape(SHAPE_RHOMBIC_DODECAHEDRON);
		outfilename = argv[2];
	} else if (!strcmp(argv[1], "rt") && argc > 2) {
		mesh = make_fixed_shape(SHAPE_RHOMBIC_TRIACONTAHEDRON);
		outfilename = argv[2];
	} else if (!strcmp(argv[1], "sor") && argc > 4) {
		FILE *f = fopen(argv[3], "r");
		if (!f) {
			TriMesh::eprintf("Couldn't open %s\n", argv[3]);
			exit(1);
		}
		vector<point> pts;
		float x, y, z;
		while (1) {
			if (fscanf(f, "%f%f%f", &x, &y, &z) != 3)
				break;
			pts.push_back(point(x,y,z));
		}
		fclose(f);
		if (pts.size() < 2) {
			TriMesh::eprintf("Not enough points read from %s\n", argv[3]);
			exit(1);
		}
		mesh = make_surface_of_revolution(atoi(argv[2]), pts);
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "teapot") && argc > 4) {
		mesh = make_teapot(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "teapot") && argc > 3) {
		mesh = make_teapot(atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "orgteapot") && argc > 4) {
		mesh = make_teapot(atoi(argv[2]), atoi(argv[3]), true);
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "orgteapot") && argc > 3) {
		mesh = make_teapot(atoi(argv[2]), false, true);
		outfilename = argv[3];
	} else {
		int n = sizeof(fixed_shape_names) / sizeof(fixed_shape_names[0]);
		for (int i = 0; i < n; i++) {
			if (!strcmp(argv[1], fixed_shape_names[i])) {
				mesh = make_fixed_shape((FixedShape) i);
				break;
			}
		}
	}

	if (mesh) {
		mesh->need_tstrips();
		mesh->write(outfilename);
	} else {
		usage(argv[0]);
	}

}
