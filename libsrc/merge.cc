/*
Szymon Rusinkiewicz
Princeton University

merge.cc
Routines for merging meshes and "sharing" vertices.
*/


#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "KDtree.h"
using namespace std;


namespace trimesh {

// Find separate mesh vertices that should be "shared": they lie on separate
// connected components, but they are within "tol" of each other.
void shared(TriMesh *mesh, float tol)
{
	int nv = mesh->vertices.size();
	if (nv < 2)
		return;
	mesh->need_faces();
	if (mesh->faces.empty())
		return;
	mesh->clear_tstrips();
	mesh->clear_grid();
	mesh->need_neighbors();
	mesh->need_adjacentfaces();

	// We only merge vertices on different connected components.
	// First we find those components.
	vector<int> comps, compsizes;
	find_comps(mesh, comps, compsizes, true);
	int ncomps = compsizes.size();

	// We only merge vertices on boundaries.  Find them.
	vector<bool> bdy(nv);
#pragma omp parallel for
	for (int i = 0; i < nv; i++)
		bdy[i] = mesh->is_bdy(i);

	// KD trees for closest-point acceleration
	vector<KDtree *> kd_trees(ncomps);
	vector< vector<const float *> > comp_points(ncomps);
	for (int i = 0; i < nv; i++) {
		if (!bdy[i] || mesh->adjacentfaces[i].empty())
			continue;
		int i_comp = comps[mesh->adjacentfaces[i][0]];
		comp_points[i_comp].push_back(&mesh->vertices[i][0]);
	}
	for (int i = 0; i < ncomps; i++)
		kd_trees[i] = new KDtree(comp_points[i]);

	// Mapping of which vertex is merged with which other
	vector<int> remap;
	remap.reserve(nv);

	// Find shared vertices
	float tol2 = sqr(tol);
	for (int i = 0; i < nv; i++) {
		remap.push_back(i);
		if (!bdy[i] || mesh->adjacentfaces[i].empty())
			continue;
		int i_comp = comps[mesh->adjacentfaces[i][0]];
		for (int j = 0; j < i_comp; j++) {
			const float *match = kd_trees[j]->
				closest_to_pt(mesh->vertices[i], tol2);
			if (!match)
				continue;
			int match_ind = (match - &mesh->vertices[0][0]) / 3;
			remap[i] = match_ind;
			break;
		}
	}
	mesh->clear_adjacentfaces();
	mesh->clear_neighbors();
	while (!kd_trees.empty()) {
		delete kd_trees.back();
		kd_trees.pop_back();
	}

	// Adjust remapping table to use lowest-numbered vertex in each group
	for (int i = 0; i < nv; i++) {
		if (remap[i] == i)
			continue;
		int curr = i, next = remap[curr], lowest = i;
		while (curr != next) {
			curr = next;
			next = remap[curr];
			if (curr < lowest)
				lowest = curr;
		}
		curr = i;
		next = remap[curr];
		remap[curr] = lowest;
		while (curr != next) {
			curr = next;
			next = remap[curr];
			remap[curr] = lowest;
		}
	}

	// Compact remapping table
	int next = 0;
	for (int i = 0; i < nv; i++) {
		if (remap[i] == i)
			remap[i] = next++;
		else
			remap[i] = remap[remap[i]];
	}

	// Remap everything
	remap_verts(mesh, remap);
}


// Join multiple meshes together, possibly sharing vertices within tol.
// If tol < 0, don't share vertices.
TriMesh *join(const vector<TriMesh *> &meshes, float tol)
{
	// Pre-allocate storage
	int out_nv = 0, out_nf = 0;
	bool have_colors = false, have_confidences = false, have_normals = false;
	bool have_tstrips = false;
	int nmeshes = meshes.size();
	for (int i = 0; i < nmeshes; i++) {
		out_nv += meshes[i]->vertices.size();
		meshes[i]->need_faces();
		out_nf += meshes[i]->faces.size();
		if (!meshes[i]->colors.empty())
			have_colors = true;
		if (!meshes[i]->confidences.empty())
			have_confidences = true;
		if (!meshes[i]->normals.empty())
			have_normals = true;
		if (!meshes[i]->tstrips.empty())
			have_tstrips = true;
	}

	TriMesh *outmesh = new TriMesh;
	outmesh->vertices.reserve(out_nv);
	outmesh->faces.reserve(out_nf);
	if (have_colors)
		outmesh->colors.reserve(out_nv);
	if (have_confidences)
		outmesh->confidences.reserve(out_nv);
	if (have_normals)
		outmesh->normals.reserve(out_nv);

	vector<KDtree *> kd_trees;
	vector<int> remap;
	if (tol >= 0.0f)
		remap.resize(out_nv);
	float tol2 = sqr(tol);

	// Now loop over all meshes
	for (int i = 0; i < nmeshes; i++) {
		TriMesh *m = meshes[i];
		int nv = m->vertices.size();
		int onv = outmesh->vertices.size();

		// Vertices and vertex properties
		outmesh->vertices.insert(outmesh->vertices.end(),
		                         m->vertices.begin(),
		                         m->vertices.end());

		if (have_colors && !m->colors.empty()) {
			outmesh->colors.insert(outmesh->colors.end(),
			                       m->colors.begin(),
			                       m->colors.end());
		} else if (have_colors && m->colors.empty()) {
			outmesh->colors.resize(onv + nv, Color::white());
		}

		if (have_confidences && !m->confidences.empty()) {
			outmesh->confidences.insert(outmesh->confidences.end(),
			                            m->confidences.begin(),
			                            m->confidences.end());
		} else if (have_confidences && m->confidences.empty()) {
			outmesh->confidences.resize(onv + nv, 1.0f);
		}

		if (have_normals) {
			if (m->normals.empty())
				m->need_normals();
			outmesh->normals.insert(outmesh->normals.end(),
						m->normals.begin(),
						m->normals.end());
		}

		// Append and renumber faces
		int onf = outmesh->faces.size();
		outmesh->faces.insert(outmesh->faces.end(),
		                      m->faces.begin(),
		                      m->faces.end());
		for (size_t j = onf; j < outmesh->faces.size(); j++) {
			outmesh->faces[j][0] += onv;
			outmesh->faces[j][1] += onv;
			outmesh->faces[j][2] += onv;
		}

		// Share vertices
		if (tol < 0.0f)
			continue;

		m->need_neighbors();
		m->need_adjacentfaces();
		vector<const float *> bdy_pts;
		for (size_t j = 0; j < m->vertices.size(); j++) {
			remap[onv + j] = onv + j;
			if (!m->is_bdy(j))
				continue;
			const point &p = outmesh->vertices[onv + j];
			bdy_pts.push_back(&p[0]);
			for (size_t k = 0; k < kd_trees.size(); k++) {
				const float *match = kd_trees[k]->closest_to_pt(p, tol2);
				if (!match)
					continue;
				int match_ind = (match - &outmesh->vertices[0][0]) / 3;
				remap[onv + j] = match_ind;
				bdy_pts.pop_back();
				break;
			}
		}
		if (i < nmeshes - 1)
			kd_trees.push_back(new KDtree(bdy_pts));
	}

	while (!kd_trees.empty()) {
		delete kd_trees.back();
		kd_trees.pop_back();
	}

	if (!remap.empty())
		remap_verts(outmesh, remap);

	if (have_tstrips)
		outmesh->need_tstrips();

	return outmesh;
}

} // namespace trimesh
