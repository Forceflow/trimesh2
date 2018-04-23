/*
Szymon Rusinkiewicz
Princeton University

mesh_cat.cc
Concatenate meshes together
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;
using namespace trimesh;


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s infiles... [-share tol] -o outfile\n", myname);
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 4)
		usage(argv[0]);

	vector<TriMesh *> meshes;
	const char *outfile = NULL;
	float tol = -1.0f;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-o") == 0 && i < argc-1) {
			outfile = argv[i+1];
			i++;
			continue;
		}
		if (strncmp(argv[i], "-share", 6) == 0 && i < argc-1) {
			tol = atof(argv[i+1]);
			i++;
			continue;
		}
		TriMesh *m = TriMesh::read(argv[i]);
		if (!m) {
			fprintf(stderr, "Couldn't read file %s\n", argv[i]);
			continue;
		}
		meshes.push_back(m);
	}


	if (outfile)
		join(meshes, tol)->write(outfile);
	else
		fprintf(stderr, "No output file specified\n");
}
