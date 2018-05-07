/*
Szymon Rusinkiewicz
Princeton University

diffuse.cc
Smoothing of meshes and per-vertex fields
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "timestamp.h"
using namespace std;
#define dprintf TriMesh::dprintf


namespace trimesh {

// Approximation to Gaussian...  Used in filtering
static inline float wt(const point &p1, const point &p2, float invsigma2)
{
	float d2 = invsigma2 * dist2(p1, p2);
	//return (d2 >= 4.0f) ? 0.0f : 1.0f - d2 * (0.5f - d2 * 0.0625f);
	return (d2 >= 6.25f) ? 0.0f : 1.0f - d2 * (0.32f - d2 * 0.0256f);
	//return (d2 >= 9.0f) ? 0.0f : exp(-0.5f*d2);
	//return (d2 >= 25.0f) ? 0.0f : exp(-0.5f*d2);
}
static inline float wt(const TriMesh *themesh, int v1, int v2, float invsigma2)
{
	return wt(themesh->vertices[v1], themesh->vertices[v2], invsigma2);
}


// Functor classes for adding scalar, vector, or tensor fields on the surface
template <class T>
struct AccumVec {
	const vector<T> &field;
	AccumVec(const vector<T> &field_) : field(field_)
		{}
	inline void operator() (const TriMesh *, int /* v0 */, T &f,
				float w, int v) const
	{
		f += w * field[v];
	}
};

struct AccumCurv {
	inline void operator() (const TriMesh *themesh, int v0, vec &c,
				float w, int v) const
	{
		vec ncurv;
		proj_curv(themesh->pdir1[v], themesh->pdir2[v],
		          themesh->curv1[v], 0, themesh->curv2[v],
		          themesh->pdir1[v0], themesh->pdir2[v0],
		          ncurv[0], ncurv[1], ncurv[2]);
		c += w * ncurv;
	}
};

struct AccumDCurv {
	inline void operator() (const TriMesh *themesh, int v0, Vec<4> &d,
				float w, int v) const
	{
		Vec<4> ndcurv;
		proj_dcurv(themesh->pdir1[v], themesh->pdir2[v],
		           themesh->dcurv[v],
		           themesh->pdir1[v0], themesh->pdir2[v0],
		           ndcurv);
		d += w * ndcurv;
	}
};


// Diffuse a vector field at 1 vertex, weighted by
// a Gaussian of width 1/sqrt(invsigma2)
template <class ACCUM, class T>
static void diffuse_vert_field(TriMesh *themesh,
                               vector<unsigned> &flags, unsigned &flag_curr,
                               const ACCUM &accum, int v, float invsigma2,
                               T &flt)
{
	if (themesh->neighbors[v].empty()) {
		flt = T();
		accum(themesh, v, flt, 1.0f, v);
		return;
	}

	flt = T();
	accum(themesh, v, flt, themesh->pointareas[v], v);
	float sum_w = themesh->pointareas[v];
	const vec &nv = themesh->normals[v];

	flag_curr++;
	flags[v] = flag_curr;
	vector<int> boundary = themesh->neighbors[v];
	while (!boundary.empty()) {
		int n = boundary.back();
		boundary.pop_back();
		if (flags[n] == flag_curr)
			continue;
		flags[n] = flag_curr;
		if ((nv DOT themesh->normals[n]) <= 0.0f)
			continue;
		// Gaussian weight
		float w = wt(themesh, n, v, invsigma2);
		if (w == 0.0f)
			continue;
		// Downweight things pointing in different directions
		w *= nv DOT themesh->normals[n];
		// Surface area "belonging" to each point
		w *= themesh->pointareas[n];
		// Accumulate weight times field at neighbor
		accum(themesh, v, flt, w, n);
		sum_w += w;
		for (size_t i = 0; i < themesh->neighbors[n].size(); i++) {
			int nn = themesh->neighbors[n][i];
			if (flags[nn] == flag_curr)
				continue;
			boundary.push_back(nn);
		}
	}
	if (sum_w != 0.0f) {
		flt /= sum_w;
	} else {
		flt = T();
		accum(themesh, v, flt, 1.0f, v);
	}
}


// Smooth the mesh geometry.
// XXX - this is perhaps not a great way to do this,
// but it seems to work better than most other things I've tried...
void smooth_mesh(TriMesh *themesh, float sigma)
{
	themesh->need_faces();
	diffuse_normals(themesh, 0.5f * sigma);
	int nv = themesh->vertices.size();
	int nf = themesh->faces.size();

	dprintf("\rSmoothing... ");
	timestamp t = now();

	float invsigma2 = 1.0f / sqr(sigma);

	vector<point> dflt(nv), dflt2(nv);
#pragma omp parallel
	{
		// Thread-local flags
		vector<unsigned> flags(nv);
		unsigned flag_curr = 0;

		// Main filtering step
#pragma omp for
		for (int i = 0; i < nv; i++) {
			diffuse_vert_field(themesh, flags, flag_curr,
				AccumVec<vec>(themesh->vertices),
				i, invsigma2, dflt[i]);
			// Just keep the displacement
			dflt[i] -= themesh->vertices[i];
		}

		// Slightly better small-neighborhood approximation
#pragma omp for
		for (int i = 0; i < nf; i++) {
			point c = (themesh->vertices[themesh->faces[i][0]] +
			           themesh->vertices[themesh->faces[i][1]] +
			           themesh->vertices[themesh->faces[i][2]])
				* (1.0f / 3.0f);
			for (int j = 0; j < 3; j++) {
				int v = themesh->faces[i][j];
				vec d = 0.5f * (c - themesh->vertices[v]);
				dflt[v] += themesh->cornerareas[i][j] /
				           themesh->pointareas[themesh->faces[i][j]] *
				           exp(-0.5f * invsigma2 * len2(d)) * d;
			}
		}

		// Filter displacement field
#pragma omp for
		for (int i = 0; i < nv; i++) {
			diffuse_vert_field(themesh, flags, flag_curr,
				AccumVec<point>(dflt),
				i, invsigma2, dflt2[i]);
		}

		// Update vertex positions
#pragma omp for
		for (int i = 0; i < nv; i++)
			themesh->vertices[i] += dflt[i] - dflt2[i]; // second Laplacian
	} // #pragma omp parallel

	dprintf("Done.  Filtering took %f sec.\n", now() - t);
}


// Filter a vertex using the method of [Jones et al. 2003]
static void jones_filter(TriMesh *themesh,
                         vector<unsigned> &flags, unsigned &flag_curr,
                         int v,
                         float invsigma2_1, float invsigma2_2,
                         vector<point> &oldverts)
{
	const point p = oldverts[v];
	const vec norm = themesh->normals[v];
	point &flt = themesh->vertices[v];

	flt.clear();
	float sum_w = 0.0f;

	flag_curr++;
	vector<int> boundary;
	boundary.push_back(v);
	while (!boundary.empty()) {
		int n = boundary.back();
		boundary.pop_back();
		if (flags[n] == flag_curr)
			continue;
		flags[n] = flag_curr;

		const point &q = oldverts[n];
		float w = wt(p, q, invsigma2_1);
		if (w == 0.0f)
			continue;

		point prediction = q + norm * ((p - q) DOT norm);
		w *= wt(prediction, q, invsigma2_2);
		if (w == 0.0f)
			continue;
		w *= themesh->pointareas[v];
		flt += w * prediction;
		sum_w += w;

		for (size_t i = 0; i < themesh->neighbors[n].size(); i++) {
			int nn = themesh->neighbors[n][i];
			if (flags[nn] == flag_curr)
				continue;
			boundary.push_back(nn);
		}
	}
	flt *= 1.0f / sum_w;
}


// Bilateral smoothing using the method of [Jones et al. 2003]
void bilateral_smooth_mesh(TriMesh *themesh, float sigma1, float sigma2)
{
	bool had_normals = !themesh->normals.empty();
	themesh->need_normals();
	themesh->need_pointareas();
	themesh->need_neighbors();
	int nv = themesh->vertices.size();

	diffuse_normals(themesh, 0.5f * sigma1);

	dprintf("\rSmoothing... ");
	timestamp t = now();

	float invsigma2_1 = 1.0f / sqr(sigma1);
	float invsigma2_2 = 1.0f / sqr(sigma2);

	vector<point> oldverts(themesh->vertices);
#pragma omp parallel
	{
		// Thread-local flags
		vector<unsigned> flags(nv);
		unsigned flag_curr = 0;

#pragma omp for
		for (int i = 0; i < nv; i++)
			jones_filter(themesh, flags, flag_curr,
				i, invsigma2_1, invsigma2_2, oldverts);
	}

	dprintf("Done.  Filtering took %f sec.\n", now() - t);
	themesh->normals.clear();
	if (had_normals)
		themesh->need_normals();
}


// Diffuse an arbitrary per-vertex vector field
template <class T>
void diffuse_vector(TriMesh *themesh, std::vector<T> &field, float sigma)
{
	themesh->need_normals();
	themesh->need_pointareas();
	themesh->need_neighbors();
	int nv = themesh->vertices.size();

	dprintf("\rSmoothing vector field... ");
	timestamp t = now();

	float invsigma2 = 1.0f / sqr(sigma);

	vector<T> flt(nv);
	AccumVec<T> a(field);
#pragma omp parallel
	{
		// Thread-local flags
		vector<unsigned> flags(nv);
		unsigned flag_curr = 0;

#pragma omp for
		for (int i = 0; i < nv; i++)
			diffuse_vert_field(themesh, flags, flag_curr,
				a, i, invsigma2, flt[i]);
	} // #pragma omp parallel

	field = flt;

	dprintf("Done.  Filtering took %f sec.\n", now() - t);
}


// Diffuse the normals across the mesh
void diffuse_normals(TriMesh *themesh, float sigma)
{
	themesh->need_normals();
	themesh->need_pointareas();
	themesh->need_neighbors();
	int nv = themesh->vertices.size();

	dprintf("\rSmoothing normals... ");
	timestamp t = now();

	float invsigma2 = 1.0f / sqr(sigma);

	vector<vec> nflt(nv);
	AccumVec<vec> a(themesh->normals);
#pragma omp parallel
	{
		// Thread-local flags
		vector<unsigned> flags(nv);
		unsigned flag_curr = 0;

#pragma omp for
		for (int i = 0; i < nv; i++) {
			diffuse_vert_field(themesh, flags, flag_curr,
				a, i, invsigma2, nflt[i]);
			normalize(nflt[i]);
		}
	} // #pragma omp parallel

	themesh->normals = nflt;

	dprintf("Done.  Filtering took %f sec.\n", now() - t);
}


// Diffuse the curvatures across the mesh
void diffuse_curv(TriMesh *themesh, float sigma)
{
	themesh->need_normals();
	themesh->need_pointareas();
	themesh->need_curvatures();
	themesh->need_neighbors();
	int nv = themesh->vertices.size();

	dprintf("\rSmoothing curvatures... ");
	timestamp t = now();

	float invsigma2 = 1.0f / sqr(sigma);

	vector<vec> cflt(nv);
#pragma omp parallel
	{
		// Thread-local flags
		vector<unsigned> flags(nv);
		unsigned flag_curr = 0;

#pragma omp for
		for (int i = 0; i < nv; i++)
			diffuse_vert_field(themesh, flags, flag_curr,
				AccumCurv(), i, invsigma2, cflt[i]);

#pragma omp for
		for (int i = 0; i < nv; i++)
			diagonalize_curv(themesh->pdir1[i], themesh->pdir2[i],
			                 cflt[i][0], cflt[i][1], cflt[i][2],
			                 themesh->normals[i],
			                 themesh->pdir1[i], themesh->pdir2[i],
			                 themesh->curv1[i], themesh->curv2[i]);
	} // #pragma omp parallel

	dprintf("Done.  Filtering took %f sec.\n", now() - t);
}


// Diffuse the curvature derivatives across the mesh
void diffuse_dcurv(TriMesh *themesh, float sigma)
{
	themesh->need_normals();
	themesh->need_pointareas();
	themesh->need_curvatures();
	themesh->need_dcurv();
	themesh->need_neighbors();
	int nv = themesh->vertices.size();

	dprintf("\rSmoothing curvature derivatives... ");
	timestamp t = now();

	float invsigma2 = 1.0f / sqr(sigma);

	vector< Vec<4> > dflt(nv);
#pragma omp parallel
	{
		// Thread-local flags
		vector<unsigned> flags(nv);
		unsigned flag_curr = 0;

#pragma omp for
		for (int i = 0; i < nv; i++)
			diffuse_vert_field(themesh, flags, flag_curr,
				AccumDCurv(), i, invsigma2, dflt[i]);
	} // #pragma omp parallel

	themesh->dcurv = dflt;
	dprintf("Done.  Filtering took %f sec.\n", now() - t);
}


// Instantiate a bunch of diffuse_vector forms
template void diffuse_vector< float >(TriMesh *, vector< float > &, float);
template void diffuse_vector< Vec<2,float> >(TriMesh *, vector< Vec<2,float> > &, float);
template void diffuse_vector< Vec<3,float> >(TriMesh *, vector< Vec<3,float> > &, float);
template void diffuse_vector< Vec<4,float> >(TriMesh *, vector< Vec<4,float> > &, float);

} // namespace trimesh
