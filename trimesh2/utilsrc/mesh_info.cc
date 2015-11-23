/*
Szymon Rusinkiewicz
Princeton University

Benedict Brown
Katholieke Universiteit Leuven

mesh_info.cc
Query various information about meshes
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
using namespace std;
using namespace trimesh;


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s infile ( desired_info | stat_op desired_stat )\n", myname);
	fprintf(stderr, "\nInfo:\n");
	fprintf(stderr, "	faces		Number of faces\n");
	fprintf(stderr, "	vertices	Number of vertices\n");
	fprintf(stderr, "	bbox		Bounding box\n");
	fprintf(stderr, "	csize		Bounding box center and size\n");
	fprintf(stderr, "	bsphere		Bounding sphere\n");
	fprintf(stderr, "	overlap infile2	Overlap area and RMS distance to other mesh\n");
	fprintf(stderr, "\nStatistical operations:\n");
	fprintf(stderr, "	min		Minimum\n");
	fprintf(stderr, "	max		Maximum\n");
	fprintf(stderr, "	mean		Mean\n");
	fprintf(stderr, "	meanabs		Mean of absolute values\n");
	fprintf(stderr, "	rms		Root mean square\n");
	fprintf(stderr, "	median		Median\n");
	fprintf(stderr, "	stdev		Standard deviation from mean\n");
	fprintf(stderr, "	total		Sum of all values\n");
	fprintf(stderr, "\nTypes of statistics:\n");
	fprintf(stderr, "	valence		Vertex valence (number of adjacent vertices)\n");
	fprintf(stderr, "	facearea	Surface area per face\n");
	fprintf(stderr, "	angle		Angle between triangle sides\n");
	fprintf(stderr, "	dihedral	Dihedral angle between faces\n");
	fprintf(stderr, "	edgelen		Edge length\n");
	fprintf(stderr, "	x		Vertex X coordinate\n");
	fprintf(stderr, "	y		Vertex Y coordinate\n");
	fprintf(stderr, "	z		Vertex Z coordinate\n");
	fprintf(stderr, "\n");
	exit(1);
}


int main(int argc, char *argv[])
{
	// Don't clutter the output
	TriMesh::set_verbose(0);

	// Parse command line and read mesh
	if (argc < 3)
		usage(argv[0]);

	TriMesh *mesh = TriMesh::read(argv[1]);
	if (!mesh)
		usage(argv[0]);

	// Display information
	if (argc == 3) {
		if (!strcmp(argv[2], "faces")) {
			mesh->need_faces();
			printf("%d\n", (int) mesh->faces.size());
		} else if (!strcmp(argv[2], "vertices")) {
			printf("%d\n", (int) mesh->vertices.size());
		} else if (!strcmp(argv[2], "bbox")) {
			mesh->need_bbox();
			printf("%g %g %g\n%g %g %g\n",
				mesh->bbox.min[0], mesh->bbox.min[1], mesh->bbox.min[2],
				mesh->bbox.max[0], mesh->bbox.max[1], mesh->bbox.max[2]);
		} else if (!strcmp(argv[2], "csize")) {
			mesh->need_bbox();
			printf("%g %g %g\n%g %g %g\n",
				mesh->bbox.center()[0], mesh->bbox.center()[1], mesh->bbox.center()[2],
				mesh->bbox.size()[0], mesh->bbox.size()[1], mesh->bbox.size()[2]);
		} else if (!strcmp(argv[2], "bsphere")) {
			mesh->need_bsphere();
			printf("%g %g %g\n%g\n",
				mesh->bsphere.center[0], mesh->bsphere.center[1], mesh->bsphere.center[2],
				mesh->bsphere.r);
		} else {
			usage(argv[0]);
		}

		return 0;
	}

	// Overlap calculation
	if (argc == 4 && !strcmp(argv[2], "overlap")) {
		TriMesh *mesh2 = TriMesh::read(argv[3]);
		if (!mesh)
			usage(argv[0]);
		xform xf1, xf2;
		xf1.read(xfname(argv[1]));
		xf2.read(xfname(argv[3]));
		float area = 0.0f, rmsdist = 0.0f;
		find_overlap(mesh, mesh2, xf1, xf2, area, rmsdist);
		printf("%g %g\n", area, rmsdist);

		return 0;
	}

	TriMesh::StatOp op = TriMesh::STAT_MIN;
	if (!strcmp(argv[2], "min"))
		op = TriMesh::STAT_MIN;
	else if (!strcmp(argv[2], "max"))
		op = TriMesh::STAT_MAX;
	else if (!strcmp(argv[2], "mean"))
		op = TriMesh::STAT_MEAN;
	else if (!strcmp(argv[2], "meanabs"))
		op = TriMesh::STAT_MEANABS;
	else if (!strcmp(argv[2], "rms"))
		op = TriMesh::STAT_RMS;
	else if (!strcmp(argv[2], "median"))
		op = TriMesh::STAT_MEDIAN;
	else if (!strcmp(argv[2], "stdev"))
		op = TriMesh::STAT_STDEV;
	else if (!strcmp(argv[2], "total"))
		op = TriMesh::STAT_TOTAL;
	else
		usage(argv[0]);

	TriMesh::StatVal val = TriMesh::STAT_VALENCE;
	if (!strcmp(argv[3], "valence"))
		val = TriMesh::STAT_VALENCE;
	else if (!strcmp(argv[3], "facearea"))
		val = TriMesh::STAT_FACEAREA;
	else if (!strcmp(argv[3], "angle"))
		val = TriMesh::STAT_ANGLE;
	else if (!strcmp(argv[3], "dihedral"))
		val = TriMesh::STAT_DIHEDRAL;
	else if (!strcmp(argv[3], "edgelen"))
		val = TriMesh::STAT_EDGELEN;
	else if (!strcmp(argv[3], "x"))
		val = TriMesh::STAT_X;
	else if (!strcmp(argv[3], "y"))
		val = TriMesh::STAT_Y;
	else if (!strcmp(argv[3], "z"))
		val = TriMesh::STAT_Z;
	else
		usage(argv[0]);

	printf("%g\n", mesh->stat(op, val));

	return 0;
}

