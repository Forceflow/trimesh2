/*
Szymon Rusinkiewicz
Princeton University

overlap.cc
Compute overlap area and mesh-to-mesh distance for two meshes
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "KDtree.h"
using namespace std;


namespace trimesh {

// Find the overlap area and RMS distance from mesh1 to mesh2.  Used by
// find_overlap in both directions, below
static void find_overlap_onedir(TriMesh *mesh1, TriMesh *mesh2,
				const xform &xf1, const xform &xf2,
				const KDtree *kd2,
				float &overlap_area, float &nonoverlap_area,
				float &rmsdist)
{
	overlap_area = 0.0f;
	rmsdist = 0.0f;
	float area_considered = 0.0f;

	xform xf12 = inv(xf2) * xf1;
	xform xf12r = norm_xf(xf12);
	int nv = mesh1->vertices.size();
	int nsamp = min(nv, 10000);

	for (int i = 0; i < nsamp; i++) {
		int ind = int((float) i / nsamp * nv);
		ind = clamp(ind, 0, nv-1);
		float this_area = mesh1->pointareas[ind];
		area_considered += this_area;
		point p = xf12 * mesh1->vertices[ind];
		const float *q = kd2->closest_to_pt(p);
		if (!q)
			continue;
		int ind2 = (q - (const float *) &(mesh2->vertices[0][0])) / 3;
		if (mesh2->is_bdy(ind2))
			continue;
		if (((xf12r * mesh1->normals[ind]) DOT mesh2->normals[ind2])
				<= 0.0f)
			continue;
		overlap_area += this_area;
		rmsdist += this_area *
			sqr((p - point(q)) DOT mesh2->normals[ind2]);
	}

	float total_area = mesh1->stat(TriMesh::STAT_SUM,
		TriMesh::STAT_FACEAREA);

	// Put the "m" and "s" in rms
	if (overlap_area) {
		rmsdist /= overlap_area;
		rmsdist = sqrt(rmsdist);
	}

	// Correct overlap_area for area_considered
	overlap_area *= total_area / area_considered;
	nonoverlap_area = total_area - overlap_area;
}


// Find overlap area and RMS distance between mesh1 and mesh2.
// rmsdist is unchanged if area returned as zero
void find_overlap(TriMesh *mesh1, TriMesh *mesh2,
                  const xform &xf1, const xform &xf2,
                  const KDtree *kd1, const KDtree *kd2,
                  float &area, float &rmsdist)
{
	mesh1->need_normals();
	mesh1->need_neighbors();
	mesh1->need_adjacentfaces();
	mesh1->need_pointareas();
	mesh2->need_normals();
	mesh2->need_neighbors();
	mesh2->need_adjacentfaces();
	mesh2->need_pointareas();

	float area1, area2, non1, non2, rmsdist1, rmsdist2;

	TriMesh::dprintf("Finding overlap 1->2... ");
	find_overlap_onedir(mesh1, mesh2, xf1, xf2, kd2, area1, non1, rmsdist1);
	TriMesh::dprintf("area = %g, RMS distance = %g\n", area1, rmsdist1);
	TriMesh::dprintf("Finding overlap 2->1... ");
	find_overlap_onedir(mesh2, mesh1, xf2, xf1, kd1, area2, non2, rmsdist2);
	TriMesh::dprintf("area = %g, RMS distance = %g\n", area2, rmsdist2);
	area = 0.5f * (area1 + area2);
	if (area)
		rmsdist = 0.5f * (rmsdist1 + rmsdist2);
}


// Return intersection-over-union between mesh1 and mesh2.
float iou(TriMesh *mesh1, TriMesh *mesh2,
          const xform &xf1, const xform &xf2,
	  const KDtree *kd1, const KDtree *kd2)
{
	mesh1->need_normals();
	mesh1->need_neighbors();
	mesh1->need_adjacentfaces();
	mesh1->need_pointareas();
	mesh2->need_normals();
	mesh2->need_neighbors();
	mesh2->need_adjacentfaces();
	mesh2->need_pointareas();

	float area1, area2, non1, non2, rmsdist1, rmsdist2;

	TriMesh::dprintf("Finding overlap 1->2... ");
	find_overlap_onedir(mesh1, mesh2, xf1, xf2, kd2, area1, non1, rmsdist1);
	TriMesh::dprintf("area = %g, RMS distance = %g\n", area1, rmsdist1);
	TriMesh::dprintf("Finding overlap 2->1... ");
	find_overlap_onedir(mesh2, mesh1, xf2, xf1, kd1, area2, non2, rmsdist2);
	TriMesh::dprintf("area = %g, RMS distance = %g\n", area2, rmsdist2);
	float int_area = 0.5f * (area1 + area2);
	float union_area = area1 + non1 + area2 + non2 - int_area;
	return int_area / union_area;
}


// Easy-to-use interfaces
void find_overlap(TriMesh *mesh1, TriMesh *mesh2, float &area, float &rmsdist)
{
	KDtree *kd1 = new KDtree(mesh1->vertices);
	KDtree *kd2 = new KDtree(mesh2->vertices);
	find_overlap(mesh1, mesh2, xform(), xform(), kd1, kd2, area, rmsdist);
	delete kd2;
	delete kd1;
}

void find_overlap(TriMesh *mesh1, TriMesh *mesh2,
                  const xform &xf1, const xform &xf2,
                  float &area, float &rmsdist)
{
	KDtree *kd1 = new KDtree(mesh1->vertices);
	KDtree *kd2 = new KDtree(mesh2->vertices);
	find_overlap(mesh1, mesh2, xf1, xf2, kd1, kd2, area, rmsdist);
	delete kd2;
	delete kd1;
}

float iou(TriMesh *mesh1, TriMesh *mesh2)
{
	KDtree *kd1 = new KDtree(mesh1->vertices);
	KDtree *kd2 = new KDtree(mesh2->vertices);
	float ret = iou(mesh1, mesh2, xform(), xform(), kd1, kd2);
	delete kd2;
	delete kd1;
	return ret;
}

float iou(TriMesh *mesh1, TriMesh *mesh2, const xform &xf1, const xform &xf2)
{
	KDtree *kd1 = new KDtree(mesh1->vertices);
	KDtree *kd2 = new KDtree(mesh2->vertices);
	float ret = iou(mesh1, mesh2, xf1, xf2, kd1, kd2);
	delete kd2;
	delete kd1;
	return ret;
}

} // namespace trimesh
