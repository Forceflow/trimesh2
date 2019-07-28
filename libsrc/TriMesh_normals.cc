/*
Szymon Rusinkiewicz
Princeton University

TriMesh_normals.cc
Compute per-vertex normals for TriMeshes

For meshes, uses average of per-face normals, weighted according to:
  Max, N.
  "Weights for Computing Vertex Normals from Facet Normals,"
  Journal of Graphics Tools, Vol. 4, No. 2, 1999.
unless need_normals(true) is called.

For raw point clouds, fits plane to k nearest neighbors.
*/

#include "TriMesh.h"
#include "KDtree.h"
#include "lineqn.h"
using namespace std;


namespace trimesh {

// Compute from tstrips, Max-weighted
static void normals_from_tstrips_Max(vector<int> &tstrips,
	vector<point> &vertices, vector<vec> &normals)
{
	if (tstrips.back() < 0) {
		bool flip = false;
		for (size_t i = 2; i < tstrips.size(); i++) {
			if (tstrips[i] < 0) {
				flip = false;
				i += 2;
				continue;
			}
			const point &p0 = vertices[tstrips[i-2]];
			const point &p1 = vertices[tstrips[i-1]];
			const point &p2 = vertices[tstrips[i  ]];
			vec a = p0 - p1, b = p1 - p2, c = p2 - p0;
			float l2a = len2(a), l2b = len2(b), l2c = len2(c);
			if (!l2a || !l2b || !l2c) {
				flip = !flip;
				continue;
			}
			vec facenormal = flip ? (b CROSS a) : (a CROSS b);
			normals[tstrips[i-2]] += facenormal * (1.0f / (l2a * l2c));
			normals[tstrips[i-1]] += facenormal * (1.0f / (l2b * l2a));
			normals[tstrips[i  ]] += facenormal * (1.0f / (l2c * l2b));
			flip = !flip;
		}
	} else {
		const int *t = &tstrips[0], *end = t + tstrips.size();
		while (likely(t < end)) {
			int striplen = *t - 2;
			t += 3;
			bool flip = false;
			for (int i = 0; i < striplen; i++, t++) {
				const point &p0 = vertices[*(t-2)];
				const point &p1 = vertices[*(t-1)];
				const point &p2 = vertices[* t   ];
				vec a = p0 - p1, b = p1 - p2, c = p2 - p0;
				float l2a = len2(a), l2b = len2(b), l2c = len2(c);
				if (!l2a || !l2b || !l2c) {
					flip = !flip;
					continue;
				}
				vec facenormal = flip ? (b CROSS a) : (a CROSS b);
				normals[*(t-2)] += facenormal * (1.0f / (l2a * l2c));
				normals[*(t-1)] += facenormal * (1.0f / (l2b * l2a));
				normals[* t   ] += facenormal * (1.0f / (l2c * l2b));
				flip = !flip;
			}
		}
	}
}


// Compute from tstrips, area-weighted
static void normals_from_tstrips_area(vector<int> &tstrips,
	vector<point> &vertices, vector<vec> &normals)
{
	if (tstrips.back() < 0) {
		bool flip = false;
		for (size_t i = 2; i < tstrips.size(); i++) {
			if (tstrips[i] < 0) {
				flip = false;
				i += 2;
				continue;
			}
			const point &p0 = vertices[tstrips[i-2]];
			const point &p1 = vertices[tstrips[i-1]];
			const point &p2 = vertices[tstrips[i  ]];
			vec a = p0 - p1, b = p1 - p2;
			vec facenormal = flip ? (b CROSS a) : (a CROSS b);
			normals[tstrips[i-2]] += facenormal;
			normals[tstrips[i-1]] += facenormal;
			normals[tstrips[i  ]] += facenormal;
			flip = !flip;
		}
	} else {
		const int *t = &tstrips[0], *end = t + tstrips.size();
		while (likely(t < end)) {
			int striplen = *t - 2;
			t += 3;
			bool flip = false;
			for (int i = 0; i < striplen; i++, t++) {
				const point &p0 = vertices[*(t-2)];
				const point &p1 = vertices[*(t-1)];
				const point &p2 = vertices[* t   ];
				vec a = p0 - p1, b = p1 - p2;
				vec facenormal = flip ? (b CROSS a) : (a CROSS b);
				normals[*(t-2)] += facenormal;
				normals[*(t-1)] += facenormal;
				normals[* t   ] += facenormal;
				flip = !flip;
			}
		}
	}
}


// Compute from faces, Max-weighted
static void normals_from_faces_Max(vector<TriMesh::Face> &faces,
	vector<point> &vertices, vector<vec> &normals)
{
	int nf = faces.size();
#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
		const point &p0 = vertices[faces[i][0]];
		const point &p1 = vertices[faces[i][1]];
		const point &p2 = vertices[faces[i][2]];
		vec a = p0 - p1, b = p1 - p2, c = p2 - p0;
		float l2a = len2(a), l2b = len2(b), l2c = len2(c);
		if (!l2a || !l2b || !l2c)
			continue;
		vec facenormal = a CROSS b;
		normals[faces[i][0]] += facenormal * (1.0f / (l2a * l2c));
		normals[faces[i][1]] += facenormal * (1.0f / (l2b * l2a));
		normals[faces[i][2]] += facenormal * (1.0f / (l2c * l2b));
	}
}


// Compute from faces, area-weighted
static void normals_from_faces_area(vector<TriMesh::Face> &faces,
	vector<point> &vertices, vector<vec> &normals)
{
	int nf = faces.size();
#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
		const point &p0 = vertices[faces[i][0]];
		const point &p1 = vertices[faces[i][1]];
		const point &p2 = vertices[faces[i][2]];
		vec a = p0 - p1, b = p1 - p2;
		vec facenormal = a CROSS b;
		normals[faces[i][0]] += facenormal;
		normals[faces[i][1]] += facenormal;
		normals[faces[i][2]] += facenormal;
	}
}


// Compute from points, fitting plane to k-nn
static void normals_from_points(vector<point> &vertices, vector<vec> &normals)
{
	const int k = 6;
	const vec ref(0, 0, 1);
	const float approx_eps = 0.05f;
	KDtree kd(vertices);
	int nv = vertices.size();
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		vector<const float *> knn;
		kd.find_k_closest_to_pt(knn, k, vertices[i], approx_eps);
		int actual_k = knn.size();
		if (actual_k < 3) {
			TriMesh::dprintf("Warning: not enough points for vertex %d\n", i);
			normals[i] = ref;
			continue;
		}
		// Compute covariance
		float C[3][3] = { { 0 } };
		// The KDtree does not return vertices[i] itself, so
		// these are all valid neighbors
		for (int j = 0; j < actual_k; j++) {
			vec d = point(knn[j]) - vertices[i];
			for (int l = 0; l < 3; l++)
				for (int m = 0; m < 3; m++)
					C[l][m] += d[l] * d[m];
		}
		float e[3];
		eigdc<float,3>(C, e);
		normals[i].set(C[0][0], C[1][0], C[2][0]);
		if ((normals[i] DOT ref) < 0.0f)
			normals[i] = -normals[i];
	}
}


// Compute per-vertex normals
void TriMesh::need_normals(bool simple_area_weighted /* = false */)
{
	// Nothing to do if we already have normals
	int nv = vertices.size();
	if (!nv || int(normals.size()) == nv)
		return;

	dprintf("Computing normals... ");
	normals.clear();
	normals.resize(nv);

	// TODO: direct handling of grids
	if (!tstrips.empty()) {
		if (simple_area_weighted)
			normals_from_tstrips_area(tstrips, vertices, normals);
		else
			normals_from_tstrips_Max(tstrips, vertices, normals);
	} else if (need_faces(), !faces.empty()) {
		if (simple_area_weighted)
			normals_from_faces_area(faces, vertices, normals);
		else
			normals_from_faces_Max(faces, vertices, normals);
	} else {
		normals_from_points(vertices, normals);
	}

	// Make them all unit-length
#pragma omp parallel for
	for (int i = 0; i < nv; i++)
		normalize(normals[i]);

	dprintf("Done.\n");
}

} // namespace trimesh
