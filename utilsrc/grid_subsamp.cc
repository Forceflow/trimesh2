/*
Szymon Rusinkiewicz
Princeton University

grid_subsamp.cc
Subsample a mesh grid.
*/

#include "TriMesh.h"
#include <cstdio>
#include <cstdlib>
using namespace std;
using namespace trimesh;


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s in.ply subsamp out.ply\n", myname);
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 4)
		usage(argv[0]);

	const char *infilename = argv[1], *outfilename = argv[3];
	int subsamp = atoi(argv[2]);
	if (subsamp < 2) {
		fprintf(stderr, "subsamp must be >= 2\n");
		usage(argv[0]);
	}

	TriMesh *mesh = TriMesh::read(infilename);
	if (!mesh)
		usage(argv[0]);
	if (mesh->grid.empty()) {
		fprintf(stderr, "No grid found in %s\n", infilename);
		usage(argv[0]);
	}

	TriMesh *outmesh = new TriMesh;
	outmesh->grid_width = mesh->grid_width / subsamp;
	outmesh->grid_height = mesh->grid_height / subsamp;
	if (outmesh->grid_width == 0 || outmesh->grid_height == 0) {
		fprintf(stderr, "Resized size is 0\n");
		usage(argv[0]);
	}

	int n = outmesh->grid_width * outmesh->grid_height;
	outmesh->grid.resize(n, TriMesh::GRID_INVALID);
	for (int i = 0; i < n; i++) {
		int x = i % outmesh->grid_width;
		int y = i / outmesh->grid_width;
		int ind = (subsamp * x) + (subsamp * y) * mesh->grid_width;
		int old_vert = mesh->grid[ind];
		if (old_vert < 0)
			continue;
		outmesh->grid[i] = int(outmesh->vertices.size());
		outmesh->vertices.push_back(mesh->vertices[old_vert]);
		if (!mesh->normals.empty())
			outmesh->normals.push_back(mesh->normals[old_vert]);
		if (!mesh->confidences.empty())
			outmesh->confidences.push_back(mesh->confidences[old_vert]);
	}

	outmesh->write(outfilename);
}
