/*
Szymon Rusinkiewicz
Princeton University

mesh_align.cc
Minimal interface to ICP: register two meshes given an initial guess
for their alignment.
*/


#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "ICP.h"
#include "timestamp.h"
#ifdef _WIN32
# include "wingetopt.h"
#else
# include <unistd.h>
#endif
#include <cstdio>
using namespace std;
using namespace trimesh;


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-options] mesh1.ply mesh2.ply\n", myname);
	fprintf(stderr, "Reads transforms in mesh1.xf and mesh2.xf, updates the latter\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "	-a		Align using affine xform\n");
	fprintf(stderr, "	-b		Bulk mode: overlap checking, write to mesh1--mesh2.xf\n");
	fprintf(stderr, "	-r		Align using rigid-body transform (default)\n");
	fprintf(stderr, "	-s		Align using rigid + isotropic scale\n");
	fprintf(stderr, "	-t		Align using translation only\n");
	fprintf(stderr, "	-v		Verbose\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	ICP_xform_type xform_type = ICP_RIGID;
	int verbose = 0;
	bool bulkmode = false;

	int c;
	while ((c = getopt(argc, argv, "harstvb")) != EOF) {
		switch (c) {
			case 'a': xform_type = ICP_AFFINE; break;
			case 'r': xform_type = ICP_RIGID; break;
			case 's': xform_type = ICP_SIMILARITY; break;
			case 't': xform_type = ICP_TRANSLATION; break;
			case 'v': verbose = 2; break;
			case 'b': bulkmode = true; break;
			default: usage(argv[0]);
		}
	}

	TriMesh::set_verbose(verbose);

	if (argc - optind < 2)
		usage(argv[0]);
	const char *filename1 = argv[optind], *filename2 = argv[optind+1];

	TriMesh *mesh1 = TriMesh::read(filename1);
	if (!mesh1)
		usage(argv[0]);
	TriMesh *mesh2 = TriMesh::read(filename2);
	if (!mesh2)
		usage(argv[0]);

	xform xf1;
	string xffilename1 = xfname(filename1);
	xf1.read(xffilename1);

	xform xf2;
	string xffilename2 = xfname(filename2);
	xf2.read(xffilename2);

	timestamp t = now();
	KDtree *kd1 = new KDtree(mesh1->vertices);
	KDtree *kd2 = new KDtree(mesh2->vertices);
	if (verbose > 1) {
		TriMesh::dprintf("Building KDtrees: %.3f msec.\n",
			(now() - t) * 1000.0f);
	}
	vector<float> weights1, weights2;

	if (bulkmode) {
		float area1 = mesh1->stat(TriMesh::STAT_SUM, TriMesh::STAT_FACEAREA);
		float area2 = mesh2->stat(TriMesh::STAT_SUM, TriMesh::STAT_FACEAREA);
		float overlap_area, overlap_dist;
		find_overlap(mesh1, mesh2, xf1, xf2, kd1, kd2,
			overlap_area, overlap_dist);
		float frac_overlap = overlap_area / min(area1, area2);
		if (frac_overlap < 0.1f) {
			TriMesh::eprintf("Insufficient overlap\n");
			exit(1);
		} else {
			TriMesh::dprintf("%.1f%% overlap\n",
				frac_overlap * 100.0f);
		}
	}

	float err = ICP(mesh1, mesh2, xf1, xf2, kd1, kd2, weights1, weights2,
			0.0f, verbose, xform_type);

	if (err < 0.0f) {
		TriMesh::eprintf("ICP failed\n");
		exit(1);
	}

	TriMesh::eprintf("ICP succeeded - distance = %f\n", err);
	if (bulkmode) {
		string xffilename12 = filename1;
		size_t dot = xffilename12.rfind(".", xffilename12.length());
		if (dot != string::npos)
			xffilename12.erase(dot);
		xffilename12 += string("--") + replace_ext(filename2, "xf");
		xform xf12 = inv(xf2) * xf1;
		xf12.write(xffilename12);
	} else {
		xf2.write(xffilename2);
	}
}
