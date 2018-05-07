/*
Szymon Rusinkiewicz
Princeton University

mesh_shade.cc
Apply procedural shaders to a mesh
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "strutil.h"
#include "KDtree.h"
#include "lineqn.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <algorithm>
using namespace std;
using namespace trimesh;

#define BIGNUM 3.3e33f
#define ATOF(x) ((float) atof(x))


// Apply a solid color to the mesh
void solidcolor(TriMesh *mesh, const char *col)
{
	unsigned c;
	sscanf(col, "%x", &c);
	int r = (c >> 16) & 0xff;
	int g = (c >> 8)  & 0xff;
	int b =  c        & 0xff;
	Color cc = Color(r,g,b);
	int nv = mesh->vertices.size();
	for (int i = 0; i < nv; i++)
		mesh->colors[i] = cc;
}


// Color based on normals
void colorbynormals(TriMesh *mesh)
{
	mesh->need_normals();
	int nv = mesh->vertices.size();
	for (int i = 0; i < nv; i++) {
		mesh->colors[i] = Color(0.5f, 0.5f, 0.5f) +
			0.5f * mesh->normals[i];
	}
}


// Color based on confidences
void colorbyconfidences(TriMesh *mesh, float conf_scale)
{
	int nv = mesh->vertices.size();
	int nc = mesh->confidences.size();
	if (nc != nv)
		return;
	for (int i = 0; i < nv; i++) {
		float c = conf_scale * mesh->confidences[i];
		mesh->colors[i] = Color(0.5f + 0.5f * c, c, c);
	}
}


// Compute a "typical scale" for the mesh: computed as 1% of
// the reciprocal of the 10-th percentile curvature
float typical_scale(TriMesh *mesh)
{
	const float frac = 0.1f;
	const float mult = 0.01f;

	mesh->need_faces();
	if (mesh->faces.empty()) {
		mesh->need_bsphere();
		float f = mult * mesh->bsphere.r;
		TriMesh::dprintf("Typical scale = %f\n", f);
		return f;
	}

	mesh->need_curvatures();
	int nv = mesh->curv1.size();
	int nsamp = min(nv, 500);

	vector<float> samples;
	samples.reserve(nsamp * 2);

	for (int i = 0; i < nsamp; i++) {
		int ind = uniform_rnd(nv);
		samples.push_back(fabs(mesh->curv1[ind]));
		samples.push_back(fabs(mesh->curv2[ind]));
	}

	int which = int(frac * samples.size());
	nth_element(samples.begin(), samples.begin() + which, samples.end());

	float f;
	if (samples[which] == 0.0f) {
		mesh->need_bsphere();
		f = mult * mesh->bsphere.r;
	} else {
		f = mult / samples[which];
	}
	TriMesh::dprintf("Typical scale = %f\n", f);
	return f;
}


// Color based on curvature
void colorbycurv(TriMesh *mesh, const char *scale, const char *smooth)
{
	mesh->need_curvatures();
	float smoothsigma = ATOF(smooth);
	if (smoothsigma > 0.0f) {
		smoothsigma *= mesh->feature_size();
		diffuse_curv(mesh, smoothsigma);
	}
	float cscale = 10.0f * ATOF(scale) * typical_scale(mesh);
	TriMesh::dprintf("Using scale = %f\n", cscale);
	cscale = sqr(cscale);

	int nv = mesh->vertices.size();
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		float H = 0.5f * (mesh->curv1[i] + mesh->curv2[i]);
		float K = mesh->curv1[i] * mesh->curv2[i];
		float h = 4.0f / 3.0f * fabs(atan2(H*H-K,H*H*sgn(H)));
		float s = M_2_PIf * atan((2.0f*H*H-K)*cscale);
		mesh->colors[i] = Color::hsv(h,s,1.0);
	}
}


// Color based on curvature.  Similar to above, but uses a grayscale mapping.
void gcolorbycurv(TriMesh *mesh, const char *scale, const char *smooth)
{
	mesh->need_curvatures();
	float smoothsigma = ATOF(smooth);
	if (smoothsigma > 0.0f) {
		smoothsigma *= mesh->feature_size();
		diffuse_curv(mesh, smoothsigma);
	}
	float cscale = 10.0f * ATOF(scale) * typical_scale(mesh);
	TriMesh::dprintf("Using scale = %f\n", cscale);

	int nv = mesh->vertices.size();
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		float H = 0.5f * (mesh->curv1[i] + mesh->curv2[i]);
		mesh->colors[i] = Color(atan(H*cscale) / M_PIf + 0.5f);
	}
}


// Accessibility shading
void acc(TriMesh *mesh, const char *maxsize_, const char *offset_)
{
	mesh->need_normals();
	float ts = typical_scale(mesh);
	float maxsize = ATOF(maxsize_) * ts;
	float offset = ATOF(offset_) * ts;
	TriMesh::dprintf("Using maxsize = %f, offset = %f\n", maxsize, offset);

	KDtree *kd = new KDtree(mesh->vertices);
	int nv = mesh->vertices.size();
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		const vec &n = mesh->normals[i];
		point p = mesh->vertices[i] + offset * n;
		float tmin = 0, tmax = maxsize;
		for (int iter = 0; iter < 8; iter++) {
			float tmid = 0.5f * (tmin + tmax);
			point q = p + tmid * n;
			const float *qq = kd->closest_to_pt(q, sqr(tmid));
			if (qq)
				tmax = tmid;
			else
				tmin = tmid;
		}
		mesh->colors[i] = Color(0.5f * (tmin + tmax) / maxsize);
	}
	delete kd;
}


// Color by distance to bdy
void bdyshade(TriMesh *mesh, const char *nedges_)
{
	int nedges = atoi(nedges_) + 1;
	int nv = mesh->vertices.size();
	mesh->need_neighbors();
	mesh->flags.resize(nv);
	for (int i = 0; i < nv; i++)
		mesh->flags[i] = mesh->is_bdy(i) ? 0 : nedges;
	for (int iter = 1; iter < nedges; iter++) {
		for (int i = 0; i < nv; i++) {
			for (size_t j = 0; j < mesh->neighbors[i].size(); j++) {
				int n = mesh->neighbors[i][j];
				if (mesh->flags[n] + 1 < mesh->flags[i])
					mesh->flags[i] = mesh->flags[n] + 1;
			}
		}
	}
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		float c = (float) mesh->flags[i] / nedges;
		mesh->colors[i] = Color(1.0f, c, c);
	}
}


// Helper for dist2mesh:
// Find closest point to p on segment from v0 to v1
point closest_on_segment(const point &v0, const point &v1, const point &p)
{
	vec v01 = v1 - v0;
	float d = (p - v0) DOT v01;
	d /= len2(v01);
	if (d < 0.0f)
		d = 0.0f;
	else if (d > 1.0f)
		d = 1.0f;
	return v0 + d * v01;
}


// Helper for dist2mesh:
// Find closest point to p on face i of mesh
point closest_on_face(const TriMesh *mesh, int i, const point &p)
{
	const TriMesh::Face &f = mesh->faces[i];
	const point &v0 = mesh->vertices[f[0]];
	const point &v1 = mesh->vertices[f[1]];
	const point &v2 = mesh->vertices[f[2]];
	vec a = v1 - v0, b = v2 - v0, p1 = p - v0, n = a CROSS b;

	float A[3][3] = { { a[0], b[0], n[0] },
	                  { a[1], b[1], n[1] },
	                  { a[2], b[2], n[2] } };
	float x[3] = { p1[0], p1[1], p1[2] };
	int indx[3];
	ludcmp<float,3>(A, indx);
	lubksb<float,3>(A, indx, x);

	if (x[0] >= 0.0f && x[1] >= 0.0f && x[0] + x[1] <= 1.0f)
		return v0 + x[0] * a + x[1] * b;

	point c01 = closest_on_segment(v0, v1, p);
	point c12 = closest_on_segment(v1, v2, p);
	point c20 = closest_on_segment(v2, v0, p);
	float d01 = dist2(c01, p);
	float d12 = dist2(c12, p);
	float d20 = dist2(c20, p);
	if (d01 < d12) {
		if (d01 < d20) return c01; else return c20;
	} else {
		if (d12 < d20) return c12; else return c20;
	}
}


// Helper for dist2mesh:
// Find (good approximation to) closest point on mesh to p.
// Finds closest vertex, then checks all faces that touch it.
bool find_closest_pt(const TriMesh *mesh, const KDtree *kd, const point &p,
                     float maxdist2, point &pmatch)
{
	// const float *match = kd->closest_to_pt(p, maxdist2);
	// The closest vertex might be much further away than the closest
	// point on the surface.  So, we need to be conservative here.
	const float *match = kd->closest_to_pt(p, 100.0f * maxdist2);
	if (!match)
		return false;
	int ind = (match - (const float *) &(mesh->vertices[0][0])) / 3;
	int nv = mesh->vertices.size();
	if (ind < 0 || ind >= nv)
		return false;

	const vector<int> &a = mesh->adjacentfaces[ind];
	if (a.empty()) {
		pmatch = mesh->vertices[ind];
		return true;
	}

	float closest_dist2 = maxdist2;
	for (size_t i = 0; i < a.size(); i++) {
		point c = closest_on_face(mesh, a[i], p);
		float this_dist2 = dist2(c, p);
		if (this_dist2 < closest_dist2) {
			closest_dist2 = this_dist2;
			pmatch = c;
		}
	}
	return (closest_dist2 != maxdist2);
}


// Color by distance to another mesh
void dist2mesh(TriMesh *mesh, const char *filename, const char *maxdist_)
{
	TriMesh *othermesh = TriMesh::read(filename);
	if (!othermesh) {
		TriMesh::eprintf("Couldn't read %s\n", filename);
		exit(1);
	}
	othermesh->need_adjacentfaces();
	KDtree *kd = new KDtree(othermesh->vertices);

	float maxdist = ATOF(maxdist_);
	float maxdist2 = sqr(maxdist);

	int nv = mesh->vertices.size();
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		const point &p = mesh->vertices[i];
		point pmatch;
		float d = maxdist;
		if (find_closest_pt(othermesh, kd, p, maxdist2, pmatch))
			d = dist(p, pmatch);
		d /= maxdist;
		float H = 4.0f * (1.0f - d);
		float S = 0.7f + 0.3f * d;
		float V = 0.7f + 0.3f * d;
		mesh->colors[i] = Color::hsv(H,S,V);
	}
	delete kd;
	delete othermesh;
}


// Color by distance to given vertex
void findvert(TriMesh *mesh, const char *v_, const char *nedges_)
{
	int v = atoi(v_);
	int nedges = atoi(nedges_) + 1;
	int nv = mesh->vertices.size();
	mesh->need_neighbors();
	mesh->flags.resize(nv);
	for (int i = 0; i < nv; i++)
		mesh->flags[i] = (i == v) ? 0 : nedges;
	for (int iter = 1; iter < nedges; iter++) {
		for (int i = 0; i < nv; i++) {
			for (size_t j = 0; j < mesh->neighbors[i].size(); j++) {
				int n = mesh->neighbors[i][j];
				if (mesh->flags[n] + 1 < mesh->flags[i])
					mesh->flags[i] = mesh->flags[n] + 1;
			}
		}
	}
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		float c = (float) mesh->flags[i] / nedges;
		mesh->colors[i] = Color(1.0f, c, c);
	}
}


// Adjust colors
void remapcolor(TriMesh *mesh, const char *scale_, const char *off_,
		const char *gamma_)
{
	float scale = ATOF(scale_);
	float off = ATOF(off_);
	float gamma = 1.0f / ATOF(gamma_);
	int nv = mesh->vertices.size();
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		Color &c = mesh->colors[i];
		c[0] = pow(c[0] * scale + off, gamma);
		c[1] = pow(c[1] * scale + off, gamma);
		c[2] = pow(c[2] * scale + off, gamma);
	}
}


// Color based on depth in direction (x,y,z)
// To find range, eliminates percentage p of points
void colorbydepth(TriMesh *mesh, const char *x_, const char *y_,
                  const char *z_, const char *p_)
{
	vec dir(ATOF(x_), ATOF(y_), ATOF(z_));
	float p = ATOF(p_);

	int nv = mesh->vertices.size();
	vector<float> depths(nv);
#pragma omp parallel for
	for (int i = 0; i < nv; i++)
		depths[i] = mesh->vertices[i] DOT dir;

	float mind, maxd;
	if (p > 0.0f) {
		int which = int(p * nv);
		nth_element(depths.begin(), depths.begin() + which, depths.end());
		mind = depths[which];

		which = nv - 1 - which;
		nth_element(depths.begin(), depths.begin() + which, depths.end());
		maxd = depths[which];
	} else {
		mind = *min_element(depths.begin(), depths.end());
		maxd = *max_element(depths.begin(), depths.end());
	}
	float mult = 1.0f / (maxd - mind);
	for (int i = 0; i < nv; i++) {
		float d = mesh->vertices[i] DOT dir;
		mesh->colors[i] = Color(mult * (d - mind));
	}
}


// Color based on logical position in grid
void colorbygridpos(TriMesh *mesh)
{
	int w = mesh->grid_width, h = mesh->grid_height;
	if (w <= 0 || h <= 0 || (int) mesh->grid.size() != w * h) {
		TriMesh::dprintf("No grid!\n");
		return;
	}

	mesh->colors.clear();
	mesh->colors.resize(mesh->vertices.size());
	for (int i = 0; i < w * h; i++) {
		int ind = mesh->grid[i];
		if (ind < 0)
			continue;
		if (mesh->colors[ind] != Color(0,0,0))
			TriMesh::dprintf("Multiple references to vertex %d\n", ind);
		float red = float((i % w) + 1) / w;
		float green = float((i / w) + 1) / h;
		mesh->colors[ind] = Color(red, green, 0.0f);
	}
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s model shader [options] outfile\n", myname);
	fprintf(stderr, "Shaders:\n");
	fprintf(stderr, "	color rrggbb	Solid color (specified in hex)\n");
	fprintf(stderr, "	normals		Colored based on normals (r = nx, etc.)\n");
	fprintf(stderr, "	confidences [s]	Color by confidences with optional scale\n");
	fprintf(stderr, "	curv sc sm	Colored based on curvature (args = scale, smoothing)\n");
	fprintf(stderr, "	gcurv sc sm	Grayscale based on curvature (args = scale, smoothing)\n");
	fprintf(stderr, "	acc max off	Accessibility (args = maximum size, offset)\n");
	fprintf(stderr, "	bdy max		Distance to boundary (arg = maximum # edges)\n");
	fprintf(stderr, "	dist m.ply max	Distance to another mesh (args = mesh, max distance)\n");
	fprintf(stderr, "	findvert v max	Distance to vertex (args = vert #, max distance)\n");
	fprintf(stderr, "	remap s o g	Remap existing colors (args = scale, offset, gamma)\n");
	fprintf(stderr, "	depth x y z p	Color by depth in dir (x,y,z); clamp fraction p\n");
	fprintf(stderr, "	gridpos		Color by row/col position in a grid (grids only!)\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 4)
		usage(argv[0]);
	const char *infilename = argv[1];

	TriMesh *themesh = TriMesh::read(infilename);
	if (!themesh)
		usage(argv[0]);
	themesh->colors.resize(themesh->vertices.size());

	const char *shader = argv[2];
	const char *outfilename = argv[3];
	if (begins_with(shader, "col")) {
		if (argc < 5) {
			TriMesh::eprintf("\n\"color\" needs one argument\n\n");
			usage(argv[0]);
		}
		solidcolor(themesh, argv[3]);
		outfilename = argv[4];
	} else if (begins_with(shader, "norm")) {
		colorbynormals(themesh);
	} else if (begins_with(shader, "conf")) {
		float conf_scale = 1.0f;
		if (argc > 4) {
			conf_scale = atof(argv[3]);
			outfilename = argv[4];
		}
		colorbyconfidences(themesh, conf_scale);
	} else if (begins_with(shader, "curv")) {
		if (argc < 6) {
			TriMesh::eprintf("\n\"curv\" needs two arguments\n\n");
			usage(argv[0]);
		}
		colorbycurv(themesh, argv[3], argv[4]);
		outfilename = argv[5];
	} else if (begins_with(shader, "gcurv")) {
		if (argc < 6) {
			TriMesh::eprintf("\n\"gcurv\" needs two arguments\n\n");
			usage(argv[0]);
		}
		gcolorbycurv(themesh, argv[3], argv[4]);
		outfilename = argv[5];
	} else if (begins_with(shader, "acc")) {
		if (argc < 6) {
			TriMesh::eprintf("\n\"acc\" needs two arguments\n\n");
			usage(argv[0]);
		}
		acc(themesh, argv[3], argv[4]);
		outfilename = argv[5];
	} else if (begins_with(shader, "bdy")) {
		if (argc < 5) {
			TriMesh::eprintf("\n\"bdy\" needs one argument\n\n");
			usage(argv[0]);
		}
		bdyshade(themesh, argv[3]);
		outfilename = argv[4];
	} else if (begins_with(shader, "dist")) {
		if (argc < 6) {
			TriMesh::eprintf("\n\"dist\" needs two arguments\n\n");
			usage(argv[0]);
		}
		dist2mesh(themesh, argv[3], argv[4]);
		outfilename = argv[5];
	} else if (begins_with(shader, "findv")) {
		if (argc < 6) {
			TriMesh::eprintf("\n\"findvert\" needs two arguments\n\n");
			usage(argv[0]);
		}
		findvert(themesh, argv[3], argv[4]);
		outfilename = argv[5];
	} else if (begins_with(shader, "remap")) {
		if (argc < 7) {
			TriMesh::eprintf("\n\"remap\" needs three arguments\n\n");
			usage(argv[0]);
		}
		remapcolor(themesh, argv[3], argv[4], argv[5]);
		outfilename = argv[6];
	} else if (begins_with(shader, "depth")) {
		if (argc < 8) {
			TriMesh::eprintf("\n\"depth\" needs four arguments\n\n");
			usage(argv[0]);
		}
		colorbydepth(themesh, argv[3], argv[4], argv[5], argv[6]);
		outfilename = argv[7];
	} else if (begins_with(shader, "grid")) {
		colorbygridpos(themesh);
	} else {
		TriMesh::eprintf("\nUnknown shader [%s]\n\n", shader);
		usage(argv[0]);
	}

	themesh->write(outfilename);
}
