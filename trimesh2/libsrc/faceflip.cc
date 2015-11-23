/*
Szymon Rusinkiewicz
Princeton University

faceflip.cc
Flip the order of vertices in each face.  Turns the mesh inside out.
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#define dprintf TriMesh::dprintf
using namespace std;


namespace trimesh {

void faceflip(TriMesh *mesh)
{
	bool had_tstrips = !mesh->tstrips.empty();
	mesh->need_faces();
	mesh->tstrips.clear();

	dprintf("Flipping faces... ");
	int nf = mesh->faces.size();
#pragma omp parallel for
	for (int i = 0; i < nf; i++)
		swap(mesh->faces[i][0], mesh->faces[i][2]);
	dprintf("Done.\n");

	if (had_tstrips)
		mesh->need_tstrips();
}

}; // namespace trimesh
