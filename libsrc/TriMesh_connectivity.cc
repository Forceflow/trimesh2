/*
Szymon Rusinkiewicz
Princeton University

TriMesh_connectivity.cc
Manipulate data structures that describe connectivity between faces and verts.
*/


#include "TriMesh.h"
using namespace std;


namespace trimesh {

// Find the direct neighbors of each vertex
void TriMesh::need_neighbors()
{
	if (!neighbors.empty())
		return;

	need_faces();
	if (faces.empty())
		return;

	dprintf("Finding vertex neighbors... ");
	int nv = vertices.size(), nf = faces.size();

	vector<int> numneighbors(nv);
	for (int i = 0; i < nf; i++) {
		numneighbors[faces[i][0]]++;
		numneighbors[faces[i][1]]++;
		numneighbors[faces[i][2]]++;
	}

	neighbors.resize(nv);
	for (int i = 0; i < nv; i++)
		neighbors[i].reserve(numneighbors[i]);

	for (int i = 0; i < nf; i++) {
		for (int j = 0; j < 3; j++) {
			vector<int> &me = neighbors[faces[i][j]];
			int n1 = faces[i][NEXT_MOD3(j)];
			int n2 = faces[i][PREV_MOD3(j)];
			if (find(me.begin(), me.end(), n1) == me.end())
				me.push_back(n1);
			if (find(me.begin(), me.end(), n2) == me.end())
				me.push_back(n2);
		}
	}

	dprintf("Done.\n");
}


// Find the faces touching each vertex
void TriMesh::need_adjacentfaces()
{
	if (!adjacentfaces.empty())
		return;

	need_faces();
	if (faces.empty())
		return;

	dprintf("Finding vertex to triangle maps... ");
	int nv = vertices.size(), nf = faces.size();

	vector<int> numadjacentfaces(nv);
	for (int i = 0; i < nf; i++) {
		numadjacentfaces[faces[i][0]]++;
		numadjacentfaces[faces[i][1]]++;
		numadjacentfaces[faces[i][2]]++;
	}

	adjacentfaces.resize(vertices.size());
	for (int i = 0; i < nv; i++)
		adjacentfaces[i].reserve(numadjacentfaces[i]);

	for (int i = 0; i < nf; i++) {
		for (int j = 0; j < 3; j++)
			adjacentfaces[faces[i][j]].push_back(i);
	}

	dprintf("Done.\n");
}


// Find the face across each edge from each other face (-1 on boundary)
// If topology is bad, not necessarily what one would expect...
void TriMesh::need_across_edge()
{
	if (!across_edge.empty())
		return;

	need_adjacentfaces();
	if (adjacentfaces.empty())
		return;

	dprintf("Finding across-edge maps... ");

	int nf = faces.size();
	across_edge.resize(nf, Face(-1,-1,-1));

#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
		for (int j = 0; j < 3; j++) {
			int v1 = faces[i][NEXT_MOD3(j)];
			int v2 = faces[i][PREV_MOD3(j)];
			const vector<int> &a1 = adjacentfaces[v1];
			for (size_t k1 = 0; k1 < a1.size(); k1++) {
				int other = a1[k1];
				if (other == i)
					continue;
				int v2_in_other = faces[other].indexof(v2);
				if (v2_in_other < 0)
					continue;
				if (faces[other][NEXT_MOD3(v2_in_other)] != v1)
					continue;
				across_edge[i][j] = other;
				break;
			}
		}
	}

	dprintf("Done.\n");
}

} // namespace trimesh
