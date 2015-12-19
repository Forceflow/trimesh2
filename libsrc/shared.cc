/*
Szymon Rusinkiewicz
Princeton University

shared.cc
Find separate mesh vertices that should be "shared": they lie on separate
connected components, but they are very close to each other.
*/


#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <vector>
using namespace std;


namespace trimesh {

// Merge vertices within tol
void shared(TriMesh *mesh, float tol)
{
	int nv = mesh->vertices.size();
	if (nv < 2)
		return;
	if (mesh->faces.empty())
		return;
	mesh->tstrips.clear();
	mesh->need_neighbors();
	mesh->need_adjacentfaces();

	// We only merge vertices on different connected components.
	// First we find those components.
	vector<int> comps, compsizes;
	find_comps(mesh, comps, compsizes, true);

	// Find boundary vertices
	vector<bool> bdy(nv);
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		bdy[i] = mesh->is_bdy(i);
	}

	vector<int> remap(nv);
	float tol2 = sqr(tol);

	// TODO: use KD tree to avoid O(n^2)
	int next = 0;
	for (int i = 0; i < nv; i++) {
		remap[i] = next++;
		if (!bdy[i] || mesh->adjacentfaces[i].empty())
			continue;
		for (int j = 0; j < i; j++) {
			if (!bdy[j] || mesh->adjacentfaces[j].empty())
				continue;
			if (comps[mesh->adjacentfaces[i][0]] ==
			    comps[mesh->adjacentfaces[j][0]])
				continue;
			if (dist2(mesh->vertices[i], mesh->vertices[j]) > tol2)
				continue;
			remap[i] = remap[j];
			next--;
			break;
		}
	}

	mesh->adjacentfaces.clear();
	mesh->neighbors.clear();
	remap_verts(mesh, remap);
	remove_unused_vertices(mesh);
	orient(mesh);
}

}; // namespace trimesh
