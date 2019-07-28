/*
Szymon Rusinkiewicz
Princeton University

mesh_check.cc
Get rid of certain kinds of topological weirdness in meshes.
Specifically, this checks for and eliminates the following:
1. Faces with invalid indices
2. Degenerate faces
3. Duplicate faces
4. Unreferenced vertices (except for pure point clouds with no faces at all)
5. Vertices with weird (NaN or really huge) coordinates
6. Double-sided fins (pairs of mirror-image faces)
7. Single-sided fins (single triangles sticking out of a manifold)
8. More than two faces at an edge, or two faces that have the same vertex
ordering at an edge.

Things it does NOT check for, but (in the ideal world) should:
 - Two manifolds touching at a vertex
 - Nonsimple boundaries

Return code is:
 0 if everything went well, mesh was OK
 1 if we couldn't read the mesh
 2 if we found problems in the mesh and modified it
*/


#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
using namespace std;
using namespace trimesh;


// Vertices bigger than this are flagged as invalid
// Change this if you actually deal with values of this order of magnitude...
#define MAX_VERTEX_COORD 1.0e10f


// Check whether vertex i has any invalid coordinates
bool bad_vertex(const TriMesh *themesh, int i)
{
	const float &x = themesh->vertices[i][0];
	const float &y = themesh->vertices[i][1];
	const float &z = themesh->vertices[i][2];

	// Remember, comparing a NaN with anything yields false
	if ((x > -MAX_VERTEX_COORD) && (x < MAX_VERTEX_COORD) &&
	    (y > -MAX_VERTEX_COORD) && (y < MAX_VERTEX_COORD) &&
	    (z > -MAX_VERTEX_COORD) && (z < MAX_VERTEX_COORD))
		return false;
	else
		return true;
}


// Check whether face i has any invalid vertices or is degenerate
bool bad_face(const TriMesh *themesh, int i)
{
	const int &i0 = themesh->faces[i][0];
	const int &i1 = themesh->faces[i][1];
	const int &i2 = themesh->faces[i][2];
	const int nv = themesh->vertices.size();

	if ((i0 >= 0) && (i0 < nv) &&
	    (i1 >= 0) && (i1 < nv) &&
	    (i2 >= 0) && (i2 < nv) &&
	    (i0 != i1) && (i1 != i2) && (i2 != i0))
		return false;
	else
		return true;
}


// Sorts the vertex indices of a face, and returns whether we changed the
// current cyclic ordering.
bool sort_indices_in_face(TriMesh::Face &f)
{
	bool flipped = false;
	if (f[1] < f[0]) {
		swap(f[0], f[1]);
		flipped = !flipped;
	}
	if (f[2] < f[0]) {
		swap(f[0], f[2]);
		flipped = !flipped;
	}
	if (f[2] < f[1]) {
		swap(f[1], f[2]);
		flipped = !flipped;
	}
	return flipped;
}


// Define an ordering on two faces
class OrderFace {
private:
	const TriMesh *themesh;
public:
	OrderFace(const TriMesh *_themesh) : themesh(_themesh)
		{}
	bool operator() (int i1, int i2) const
	{
		const TriMesh::Face &f1 = themesh->faces[i1];
		const TriMesh::Face &f2 = themesh->faces[i2];
		if (f1[0] < f2[0]) return true;
		if (f1[0] > f2[0]) return false;
		if (f1[1] < f2[1]) return true;
		if (f1[1] > f2[1]) return false;
		if (f1[2] < f2[2]) return true;
		return false;
	}
};


// Checks the connectivity at an edge.  Counts the number of polygons that
// meet at that edge, and complains if there is more than one forwards and one
// backwards.  Returns true iff this is a weird edge.
bool bad_edge(const TriMesh *themesh,
              int v1, int v2,
              const vector<bool> &bad_faces)
{
	int numforwards = 0, numbackwards = 0;

	int naf1 = themesh->adjacentfaces[v1].size();
	for (int i = 0; i < naf1; i++) {
		int ind = themesh->adjacentfaces[v1][i];
		if (bad_faces[ind])
			continue;
		const TriMesh::Face &f = themesh->faces[ind];
		for (int j = 0; j < 3; j++) {
			if (f[j] != v1)
				continue;
			if (f[NEXT_MOD3(j)] == v2)
				numforwards++;
			else if (f[PREV_MOD3(j)] == v2)
				numbackwards++;
		}
	}

	if (numforwards > 1 || numbackwards > 1) {
		printf("Edge (%d, %d) has %d forwards and %d backwards faces.\n",
			v1, v2, numforwards, numbackwards);
		return true;
	}
	return false;
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s in.ply out.ply\n", myname);
	exit(1);
}


int main(int argc, char *argv[])
{
	// Deal with command-line options
	if (argc < 3 ||
	    !strncmp(argv[1], "-h", 2) ||
	    !strncmp(argv[1], "--h", 3)) {
		usage(argv[0]);
	}
	const char *infilename = argv[1];
	const char *outfilename = argv[2];
	if (!infilename || !outfilename)
		usage(argv[0]);

	// Read in file
	TriMesh *in = TriMesh::read(infilename);
	if (!in) {
		fprintf(stderr, "Couldn't open file %s\n", infilename);
		exit(1);
	}

	// Initial sanity checks
	int nv = in->vertices.size();
	if (!nv) {
		printf("Input has no vertices!\n");
		exit(0);
	}

	bool had_tstrips = !in->tstrips.empty();
	in->need_faces();
	in->tstrips.clear();
	int nf = in->faces.size();
	if (!nf)
		printf("Warning: input has no faces!\n");


	// Get rid of vertices with weird coordinates
	bool had_problems = false;
	int num_bad_verts = 0;
	vector<bool> bad_verts(nv, false);
	printf("Looking for bad vertices... "); fflush(stdout);
	for (int i = 0; i < nv; i++) {
		if (bad_vertex(in, i)) {
			bad_verts[i] = true;
			num_bad_verts++;
		}
	}
	printf("Found %d funny vertices... Done.\n", num_bad_verts);
	if (num_bad_verts) {
		remove_vertices(in, bad_verts);
		nv = in->vertices.size();
		had_problems = true;
	}

	if (!nf) {
		in->write(outfilename);
		exit(had_problems ? 2 : 0);
	}

	// Get rid of faces with invalid vertices and degenerate faces
	int num_bad_faces = 0;
	vector<bool> bad_faces(nf, false);

	printf("Looking for bad faces... "); fflush(stdout);
	for (int i = 0; i < nf; i++) {
		if (bad_face(in, i)) {
			bad_faces[i] = true;
			num_bad_faces++;
		}
	}
	printf("Found %d bogus triangles... Done.\n", num_bad_faces);


	// To make it easier to detect duplicate faces, we're now going to
	// sort the vertex indices of each face.  When we're doing that,
	// though, we need to record whether we changed the current cyclic
	// ordering.
	printf("Looking for duplicate faces... "); fflush(stdout);
	vector<bool> flipped(nf);
	for (int i = 0; i < nf; i++) {
		if (bad_faces[i])
			continue;
		flipped[i] = sort_indices_in_face(in->faces[i]);
	}

	// OK, now sort the faces themselves
	// (actually, we just sort a list of "pointers"...)
	vector<int> facelist;
	facelist.reserve(nf - num_bad_faces);
	for (int i = 0; i < nf; i++) {
		if (!bad_faces[i])
			facelist.push_back(i);
	}
	sort(facelist.begin(), facelist.end(), OrderFace(in));

	// Go through and look for duplicates
	// We look for identical faces (which will be next to each other after
	// the sort), and accumulate +1 for each non-flipped face and -1 for
	// each flipped face.  At the end of this, if the sum is:
	//  = 0 : all the faces get deleted
	//  > 0 : we leave exactly one non-flipped face
	//  < 0 : we leave exactly one flipped face
	int s = facelist.size();
	int nuked = 0;
	for (int i = 0; i < s; i++) {
		int start = i;
		int tmp = flipped[facelist[start]] ? -1 : 1;
		while ((i+1 < s) &&
		       (in->faces[facelist[start]][2] == in->faces[facelist[i+1]][2]) &&
		       (in->faces[facelist[start]][1] == in->faces[facelist[i+1]][1]) &&
		       (in->faces[facelist[start]][0] == in->faces[facelist[i+1]][0])) {
			i++;
			tmp += flipped[facelist[i]] ? -1 : 1;
			bad_faces[facelist[i]] = true;
			nuked++;
		}
		if (tmp == 0) {
			bad_faces[facelist[start]] = true;
			nuked++;
		} else if (tmp < 0) {
			// Restore this face to its proper ordering
			swap(in->faces[facelist[start]][1],
			     in->faces[facelist[start]][2]);
		}
	}
	printf("Found %d duplicates and fins... Done.\n", nuked);

	// Actually remove the bad triangles
	if (num_bad_faces || nuked) {
		remove_faces(in, bad_faces);
		nf = in->faces.size();
		had_problems = true;
	}

	// Check for broken topology
	for (int i = 0; i < nf; i++)
		bad_faces[i] = false;
	num_bad_faces = 0;

	for (int i = 0; i < nv; i++)
		bad_verts[i] = false;
	num_bad_verts = 0;

	in->need_adjacentfaces();

	printf("Checking for non-manifold topology... "); fflush(stdout);
	vector<int> suspect_faces;
	for (int i = 0; i < nf; i++) {
		if (bad_edge(in, in->faces[i][0], in->faces[i][1], bad_faces)) {
			if (in->adjacentfaces[in->faces[i][2]].size() == 1) {
				bad_faces[i] = true;
				num_bad_faces++;
			} else {
				suspect_faces.push_back(i);
			}
		} else if (bad_edge(in, in->faces[i][1], in->faces[i][2], bad_faces)) {
			if (in->adjacentfaces[in->faces[i][0]].size() == 1) {
				bad_faces[i] = true;
				num_bad_faces++;
			} else {
				suspect_faces.push_back(i);
			}
		} else if (bad_edge(in, in->faces[i][2], in->faces[i][0], bad_faces)) {
			if (in->adjacentfaces[in->faces[i][1]].size() == 1) {
				bad_faces[i] = true;
				num_bad_faces++;
			} else {
				suspect_faces.push_back(i);
			}
		}
	}

	if (num_bad_faces) {
		printf("%d one-sided fins... ", num_bad_faces); fflush(stdout);
	}

	// If any of the faces are still bad, things are not going
	// well.  Nuke all affected vertices.
	for (size_t i = 0; i < suspect_faces.size(); i++) {
		int f = suspect_faces[i];
		if (bad_edge(in, in->faces[f][0], in->faces[f][1], bad_faces) ||
		    bad_edge(in, in->faces[f][1], in->faces[f][2], bad_faces) ||
		    bad_edge(in, in->faces[f][2], in->faces[f][0], bad_faces) ) {
			bad_verts[in->faces[f][0]] = true;
			bad_verts[in->faces[f][1]] = true;
			bad_verts[in->faces[f][2]] = true;
			num_bad_verts++;
		}
	}
	printf("Done.\n");

	in->adjacentfaces.clear();
	if (num_bad_faces) {
		remove_faces(in, bad_faces);
		had_problems = true;
	}
	if (num_bad_verts) {
		remove_vertices(in, bad_verts);
		had_problems = true;
	}

	// Remove unused vertices
	nv = in->vertices.size();
	remove_unused_vertices(in);
	int new_nv = in->vertices.size();
	if (new_nv != nv)
		had_problems = true;

	// Write it out
	if (had_tstrips) {
		in->need_tstrips();
		in->faces.clear();
	}
	in->write(outfilename);

	// Return code
	exit(had_problems ? 2 : 0);
}
