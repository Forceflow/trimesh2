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
	fprintf(stderr, "Usage: %s infile [-noxf] ( desired_info | stat_op desired_stat )\n", myname);
	fprintf(stderr, "\nInfo:\n");
	fprintf(stderr, "	faces		Number of faces\n");
	fprintf(stderr, "	vertices	Number of vertices\n");
	fprintf(stderr, "	bbox		Bounding box\n");
	fprintf(stderr, "	csize		Bounding box center and size\n");
	fprintf(stderr, "	bsphere		Bounding sphere\n");
	fprintf(stderr, "	vert_mean	Mean position (center of mass) of vertices\n");
	fprintf(stderr, "	face_mean	Mean position (center of mass) of faces\n");
	fprintf(stderr, "	vert_stdev	Standard deviation of vertices around mean\n");
	fprintf(stderr, "	face_stdev	Standard deviation of faces around mean\n");
	fprintf(stderr, "	overlap infile2	Overlap area and RMS distance to other mesh\n");
	fprintf(stderr, "	iou infile2	Intersection-over-union area with other mesh\n");
	fprintf(stderr, "\nStatistical operations:\n");
	fprintf(stderr, "	min		Minimum\n");
	fprintf(stderr, "	minabs		Minimum absolute value\n");
	fprintf(stderr, "	max		Maximum\n");
	fprintf(stderr, "	maxabs		Maximum absolute value\n");
	fprintf(stderr, "	sum		Sum of all values\n");
	fprintf(stderr, "	sumabs		Sum of absolute values\n");
	fprintf(stderr, "	sumsqr		Sum of squared values\n");
	fprintf(stderr, "	mean		Mean\n");
	fprintf(stderr, "	meanabs		Mean of absolute values\n");
	fprintf(stderr, "	rms		Root mean square\n");
	fprintf(stderr, "	median		Median\n");
	fprintf(stderr, "	stdev		Standard deviation from mean\n");
	fprintf(stderr, "\nTypes of statistics:\n");
	fprintf(stderr, "	valence		Vertex valence (number of adjacent vertices)\n");
	fprintf(stderr, "	facearea	Surface area per face\n");
	fprintf(stderr, "	angle		Angle between triangle sides\n");
	fprintf(stderr, "	dihedral	Dihedral angle between faces\n");
	fprintf(stderr, "	edgelen		Edge length\n");
	fprintf(stderr, "	x		Vertex X coordinate\n");
	fprintf(stderr, "	y		Vertex Y coordinate\n");
	fprintf(stderr, "	z		Vertex Z coordinate\n");
	fprintf(stderr, "\nAutomatically reads infile.xf unless -noxf is passed\n");
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

	const char *filename = argv[1];
	const char *info_type = argv[2];
	const char *info_param = (argc > 3) ? argv[3] : NULL;
	bool use_xf = true;
	if (!strcmp(argv[2], "-noxf")) {
		if (argc < 4)
			usage(argv[0]);
		use_xf = false;
		info_type = argv[3];
		info_param = (argc > 4) ? argv[4] : NULL;
	}

	TriMesh *mesh = TriMesh::read(filename);
	if (!mesh)
		usage(argv[0]);

	if (use_xf) {
		xform xf;
		if (xf.read(xfname(filename)))
			apply_xform(mesh, xf);
	}

	// Figure out what we want
	if (!strcmp(info_type, "faces")) {
		mesh->need_faces();
		printf("%d\n", (int) mesh->faces.size());
		return 0;
	} else if (!strcmp(info_type, "vertices")) {
		printf("%d\n", (int) mesh->vertices.size());
		return 0;
	} else if (!strcmp(info_type, "bbox")) {
		mesh->need_bbox();
		printf("%g %g %g\n%g %g %g\n",
			mesh->bbox.min[0], mesh->bbox.min[1], mesh->bbox.min[2],
			mesh->bbox.max[0], mesh->bbox.max[1], mesh->bbox.max[2]);
		return 0;
	} else if (!strcmp(info_type, "csize")) {
		mesh->need_bbox();
		printf("%g %g %g\n%g %g %g\n",
			mesh->bbox.center()[0], mesh->bbox.center()[1], mesh->bbox.center()[2],
			mesh->bbox.size()[0], mesh->bbox.size()[1], mesh->bbox.size()[2]);
		return 0;
	} else if (!strcmp(info_type, "bsphere")) {
		mesh->need_bsphere();
		printf("%g %g %g\n%g\n",
			mesh->bsphere.center[0], mesh->bsphere.center[1], mesh->bsphere.center[2],
			mesh->bsphere.r);
		return 0;
	} else if (!strcmp(info_type, "vert_mean")) {
		point p = point_center_of_mass(mesh->vertices);
		printf("%g %g %g\n", p[0], p[1], p[2]);
		return 0;
	} else if (!strcmp(info_type, "face_mean")) {
		point p = mesh_center_of_mass(mesh);
		printf("%g %g %g\n", p[0], p[1], p[2]);
		return 0;
	} else if (!strcmp(info_type, "vert_stdev")) {
		trans(mesh, -point_center_of_mass(mesh->vertices));
		float C[3][3];
		point_covariance(mesh->vertices, C);
		printf("%g\n", sqrt(C[0][0] + C[1][1] + C[2][2]));
		return 0;
	} else if (!strcmp(info_type, "face_stdev")) {
		trans(mesh, -mesh_center_of_mass(mesh));
		float C[3][3];
		mesh_covariance(mesh, C);
		printf("%g\n", sqrt(C[0][0] + C[1][1] + C[2][2]));
		return 0;
	} else if (!strcmp(info_type, "overlap") && info_param) {
		TriMesh *mesh2 = TriMesh::read(info_param);
		if (!mesh2)
			usage(argv[0]);
		xform xf2;
		if (use_xf)
			xf2.read(xfname(argv[3]));
		float area = 0.0f, rmsdist = 0.0f;
		find_overlap(mesh, mesh2, xform(), xf2, area, rmsdist);
		printf("%g %g\n", area, rmsdist);
		return 0;
	} else if (!strcmp(info_type, "iou") && info_param) {
		TriMesh *mesh2 = TriMesh::read(info_param);
		if (!mesh2)
			usage(argv[0]);
		xform xf2;
		if (use_xf)
			xf2.read(xfname(argv[3]));
		printf("%g\n", iou(mesh, mesh2, xform(), xf2));
		return 0;
	}

	// If it wasn't any of those, see whether it's a StatOp
	TriMesh::StatOp op = TriMesh::STAT_MIN;
	if (!strcmp(info_type, "min"))
		op = TriMesh::STAT_MIN;
	else if (!strcmp(info_type, "minabs"))
		op = TriMesh::STAT_MINABS;
	else if (!strcmp(info_type, "max"))
		op = TriMesh::STAT_MAX;
	else if (!strcmp(info_type, "maxabs"))
		op = TriMesh::STAT_MAXABS;
	else if (!strcmp(info_type, "total"))
		op = TriMesh::STAT_SUM;
	else if (!strcmp(info_type, "sum"))
		op = TriMesh::STAT_SUM;
	else if (!strcmp(info_type, "sumabs"))
		op = TriMesh::STAT_SUMABS;
	else if (!strcmp(info_type, "sumsqr"))
		op = TriMesh::STAT_SUMSQR;
	else if (!strcmp(info_type, "mean"))
		op = TriMesh::STAT_MEAN;
	else if (!strcmp(info_type, "meanabs"))
		op = TriMesh::STAT_MEANABS;
	else if (!strcmp(info_type, "rms"))
		op = TriMesh::STAT_RMS;
	else if (!strcmp(info_type, "median"))
		op = TriMesh::STAT_MEDIAN;
	else if (!strcmp(info_type, "stdev"))
		op = TriMesh::STAT_STDEV;
	else
		usage(argv[0]);

	if (!info_param)
		usage(argv[0]);

	TriMesh::StatVal val = TriMesh::STAT_VALENCE;
	if (!strcmp(info_param, "valence"))
		val = TriMesh::STAT_VALENCE;
	else if (!strcmp(info_param, "facearea"))
		val = TriMesh::STAT_FACEAREA;
	else if (!strcmp(info_param, "angle"))
		val = TriMesh::STAT_ANGLE;
	else if (!strcmp(info_param, "dihedral"))
		val = TriMesh::STAT_DIHEDRAL;
	else if (!strcmp(info_param, "edgelen"))
		val = TriMesh::STAT_EDGELEN;
	else if (!strcmp(info_param, "x"))
		val = TriMesh::STAT_X;
	else if (!strcmp(info_param, "y"))
		val = TriMesh::STAT_Y;
	else if (!strcmp(info_param, "z"))
		val = TriMesh::STAT_Z;
	else
		usage(argv[0]);

	printf("%g\n", mesh->stat(op, val));

	return 0;
}
