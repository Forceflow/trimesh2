/*
Szymon Rusinkiewicz
Princeton University

conn_comps.cc
Determine the connected components of a mesh, and perform some basic
manipulations on them.  utilsrc/mesh_cc is a front-end for this code.
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include <vector>
#include <stack>
using namespace std;


#define NO_COMP -1
#define FOR_EACH_ADJACENT_FACE(mesh,v,f) \
	for (size_t f_ind = 0, f = mesh->adjacentfaces[v][0]; \
	     (f_ind < mesh->adjacentfaces[v].size()) && \
	     ((f = mesh->adjacentfaces[v][f_ind]) || 1); \
	     f_ind++)


namespace trimesh {

// Helper class for comparing two integers by finding the elements at those
// indices within some array and comparing them
template <class Array>
class CompareArrayElements {
private:
	const Array &a;
public:
	CompareArrayElements(const Array &_a) : a(_a)
		{}
	bool operator () (int i1, int i2) const
	{
		return (a[i1] > a[i2]);
	}
};


// Are two faces connected along an edge (or vertex)?
static bool connected(const TriMesh *mesh, int f1, int f2, bool conn_vert)
{
	int f10=mesh->faces[f1][0], f11=mesh->faces[f1][1], f12=mesh->faces[f1][2];
	int f20=mesh->faces[f2][0], f21=mesh->faces[f2][1], f22=mesh->faces[f2][2];

	if (conn_vert)
		return f10 == f20 || f10 == f21 || f10 == f22 ||
		       f11 == f20 || f11 == f21 || f11 == f22 ||
		       f12 == f20 || f12 == f21 || f12 == f22;
	else
		return (f10 == f20 && (f11 == f22 || f12 == f21)) ||
		       (f10 == f21 && (f11 == f20 || f12 == f22)) ||
		       (f10 == f22 && (f11 == f21 || f12 == f20)) ||
		       (f11 == f20 && f12 == f22) ||
		       (f11 == f21 && f12 == f20) ||
		       (f11 == f22 && f12 == f21);
}


// Helper function for find_comps, below.  Finds and marks all the faces
// connected to f.
static void find_connected(const TriMesh *mesh,
			   vector<int> &comps, vector<int> &compsizes,
			   int f, int whichcomponent, bool conn_vert)
{
	stack<int> s;
	s.push(f);
	while (!s.empty()) {
		int currface = s.top();
		s.pop();
		for (int i = 0; i < 3; i++) {
			int vert = mesh->faces[currface][i];
			FOR_EACH_ADJACENT_FACE(mesh, vert, adjface) {
				if (comps[adjface] != NO_COMP ||
				    !connected(mesh, adjface, currface, conn_vert))
					continue;
				comps[adjface] = whichcomponent;
				compsizes[whichcomponent]++;
				s.push(adjface);
			}
		}
	}
}


// Helper function for find_comps, below.  Sorts the connected components
// from largest to smallest.  Renumbers the elements of compsizes to
// reflect this new numbering.
static void sort_comps(vector<int> &comps, vector<int> &compsizes)
{
	vector<int> comp_pointers(compsizes.size());
	for (size_t i = 0; i < comp_pointers.size(); i++)
		comp_pointers[i] = i;

	sort(comp_pointers.begin(), comp_pointers.end(),
	     CompareArrayElements< vector<int> >(compsizes));

	vector<int> remap_table(comp_pointers.size());
	for (size_t i = 0; i < comp_pointers.size(); i++)
		remap_table[comp_pointers[i]] = i;
	for (size_t i = 0; i < comps.size(); i++)
		comps[i] = remap_table[comps[i]];

	vector<int> newcompsizes(compsizes.size());
	for (size_t i = 0; i < compsizes.size(); i++)
		newcompsizes[i] = compsizes[comp_pointers[i]];
	compsizes = newcompsizes;
}


// Find the connected components of TriMesh "in".
// Considers components to be connected if they touch at a vertex if
//  conn_vert == true, else they need to touch at an edge.
// Outputs:
//  comps is a vector that gives a mapping from each face to its
//   associated connected component.
//  compsizes holds the size of each connected component.
// Connected components are sorted from largest to smallest.
void find_comps(TriMesh *mesh, vector<int> &comps, vector<int> &compsizes,
		bool conn_vert /* = false */)
{
	if (mesh->vertices.empty())
		return;
	mesh->need_faces();
	if (mesh->faces.empty())
		return;
	mesh->need_adjacentfaces();

	int nf = mesh->faces.size();
	comps.clear();
	comps.reserve(nf);
	comps.resize(nf, NO_COMP);
	compsizes.clear();

	for (int i = 0; i < nf; i++) {
		if (comps[i] != NO_COMP)
			continue;
		int comp = compsizes.size();
		comps[i] = comp;
		compsizes.push_back(1);
		find_connected(mesh, comps, compsizes, i, comp, conn_vert);
	}

	if (compsizes.size() > 1)
		sort_comps(comps, compsizes);
}


// Select a particular connected component, and delete all other vertices from
// the mesh.
void select_comp(TriMesh *mesh, const vector<int> &comps, int whichcc)
{
	int numfaces = mesh->faces.size();
	vector<bool> toremove(numfaces, false);
	for (int i = 0; i < numfaces; i++) {
		if (comps[i] != whichcc)
			toremove[i] = true;
	}

	remove_faces(mesh, toremove);
	remove_unused_vertices(mesh);
}


// Select the connected components no smaller than min_size (but no more than
// total_largest components), and delete all other vertices from the mesh.
// Updates comps and compsizes.
void select_big_comps(TriMesh *mesh,
		      const vector<int> &comps, const vector<int> &compsizes,
		      int min_size,
		      int total_largest /* = std::numeric_limits<int>::max() */)
{
	int ncomp = compsizes.size();
	int keep_last = min(ncomp - 1, total_largest - 1);
	while (keep_last > -1 && compsizes[keep_last] < min_size)
		keep_last--;

	int numfaces = mesh->faces.size();
	vector<bool> toremove(numfaces, false);
	for (int i = 0; i < numfaces; i++) {
		if (comps[i] > keep_last)
			toremove[i] = true;
	}

	remove_faces(mesh, toremove);
	remove_unused_vertices(mesh);
}


// Select the connected components no bigger than max_size (but no more than
// total_smallest components), and delete all other vertices from the mesh.
void select_small_comps(TriMesh *mesh,
			const vector<int> &comps, const vector<int> &compsizes,
			int max_size,
			int total_smallest /* = std::numeric_limits<int>::max() */)
{
	int ncomp = compsizes.size();
	int keep_first = max(0, ncomp - total_smallest);
	while (keep_first < ncomp && compsizes[keep_first] > max_size)
		keep_first++;

	int numfaces = mesh->faces.size();
	vector<bool> toremove(numfaces, false);
	for (int i = 0; i < numfaces; i++) {
		if (comps[i] < keep_first)
			toremove[i] = true;
	}

	remove_faces(mesh, toremove);
	remove_unused_vertices(mesh);
}

}; // namespace trimesh
