/*
Szymon Rusinkiewicz
Princeton University

ICP.cc
Iterative Closest Point alignment using covariance-weighted sampling,
adaptive outlier rejection, and symmetric point-to-plane minimization.
*/

#include "ICP.h"
#include "KDtree.h"
#include "timestamp.h"
#include "lineqn.h"
using namespace std;


#define INITIAL_ITERS 2
#define MAX_ITERS 199
#define TERMINATION_ITER_THRESH 20
#define FINAL_ITERS 2
#define MIN_PAIRS 10
#define DESIRED_PAIRS 2000
#define DESIRED_PAIRS_FINAL 5000
#define CDF_UPDATE_INTERVAL 20
#define REJECT_BDY false
#define USE_NORMCOMPAT true
#define REGULARIZATION 0.005f
#define MEDIAN_TO_SIGMA 1.4826f
#define DIST_THRESH_MULT_START 4.0f
#define DIST_THRESH_MULT_FINAL 2.0f
#define NORMDOT_THRESH_START 0.5f
#define NORMDOT_THRESH_FINAL 0.9f
#define THRESH_RATE_OF_CHANGE 0.025f
#define dprintf TriMesh::dprintf


namespace trimesh {


// Type of ICP constraint
enum ICP_iter_type {
	ICP_POINT_TO_POINT, ICP_POINT_TO_PLANE
};


// A pair of points, with associated normals
struct PtPair {
	vec p1, n1, p2, n2;
	PtPair(const point &p1_, const vec &n1_,
	       const point &p2_, const vec &n2_) :
			p1(p1_), n1(n1_), p2(p2_), n2(n2_)
		{}
};


// A class for evaluating compatibility of normals during KDtree searches.
// This is a simplified version that just checks for positive dot product.
class NormCompat : public KDtree::CompatFunc {
private:
	TriMesh *m;
	const vec &n;

public:
	NormCompat(const vec &n_, TriMesh *m_): m(m_), n(n_)
		{}
	virtual bool operator () (const float *p) const
	{
		int idx = (const point *) p - &(m->vertices[0]);
		return (n DOT m->normals[idx]) > 0.0f;
	}
};


// A spatial grid datastructure for fast overlap computation
class Grid {
private:
	enum { GRID_SHIFT = 4, GRID_MAX = (1 << GRID_SHIFT) - 1 };
	float xmin, xmax, ymin, ymax, zmin, zmax, scale;
	vector<char> g;
public:
	inline float bbox_size() const
	{
		return dist(point(xmin, ymin, zmin), point(xmax, ymax, zmax));
	}
	inline bool valid(const point &p) const
	{
		return p[0] >= xmin && p[1] >= ymin && p[2] >= zmin &&
		       p[0] <= xmax && p[1] <= ymax && p[2] <= zmax;
	}
	inline int ind(int xcell, int ycell, int zcell) const
	{
		xcell = clamp(xcell, 0, int(GRID_MAX));
		ycell = clamp(ycell, 0, int(GRID_MAX));
		zcell = clamp(zcell, 0, int(GRID_MAX));
		return (xcell << (2 * GRID_SHIFT)) +
		       (ycell << GRID_SHIFT) +
		        zcell;
	}
	inline int ind(const point &p) const
	{
		return ind(int(scale * (p[0] - xmin)),
		           int(scale * (p[1] - ymin)),
		           int(scale * (p[2] - zmin)));
	}
	inline bool overlaps(const point &p) const
	{
		return valid(p) && g[ind(p)];
	}
	Grid(const vector<point> &pts);
};


// Compute a Grid from a list of points
Grid::Grid(const vector<point> &pts)
{
	size_t gsize = (1 << (3 * GRID_SHIFT));
	g.resize(gsize);
	if (pts.empty()) {
		xmin = xmax = ymin = ymax = zmin = zmax = scale = 0.0f;
		return;
	}

	// Find bounding box of pts
	xmin = xmax = pts[0][0];
	ymin = ymax = pts[0][1];
	zmin = zmax = pts[0][2];
	size_t npts = pts.size();
	for (size_t i = 1; i < npts; i++) {
		if      (pts[i][0] < xmin) xmin = pts[i][0];
		else if (pts[i][0] > xmax) xmax = pts[i][0];
		if      (pts[i][1] < ymin) ymin = pts[i][1];
		else if (pts[i][1] > ymax) ymax = pts[i][1];
		if      (pts[i][2] < zmin) zmin = pts[i][2];
		else if (pts[i][2] > zmax) zmax = pts[i][2];
	}
	scale = 1.0f / max(max(xmax - xmin, ymax - ymin), zmax - zmin);
	scale *= float(1 << GRID_SHIFT);

	// Set grid cells of pts
	vector<char> tmpg(1 << (3 * GRID_SHIFT));
	for (size_t i = 0; i < npts; i++)
		tmpg[ind(pts[i])] = 1;

	// Dilate
	const int xoff[27] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
	                        0,  0,  0,  0,  0,  0,  0,  0,  0,
	                        1,  1,  1,  1,  1,  1,  1,  1,  1 };
	const int yoff[27] = { -1, -1, -1,  0,  0,  0,  1,  1,  1,
	                       -1, -1, -1,  0,  0,  0,  1,  1,  1,
	                       -1, -1, -1,  0,  0,  0,  1,  1,  1 };
	const int zoff[27] = { -1,  0,  1, -1,  0,  1, -1,  0,  1,
	                       -1,  0,  1, -1,  0,  1, -1,  0,  1,
	                       -1,  0,  1, -1,  0,  1, -1,  0,  1 };
#pragma omp parallel for
	for (size_t i = 0; i < gsize; i++) {
		int x = i >> (2 * GRID_SHIFT);
		int y = (i >> GRID_SHIFT) & GRID_MAX;
		int z = i & GRID_MAX;
		for (int j = 0; j < 27; j++) {
			// ind() clamps to the edge of the grid
			if (tmpg[ind(x + xoff[j], y + yoff[j], z + zoff[j])]) {
				g[i] = 1;
				break;
			}
		}
	}
}


// Return Grid for a mesh
Grid *make_grid(TriMesh *mesh)
{
	return new Grid(mesh->vertices);
}


// Determine which points on mesh1 and mesh2 overlap the other,
// filling in o1 and o2.  Also fills in maxdist, if it is <= 0 on input.
void compute_overlaps(TriMesh *mesh1, TriMesh *mesh2,
                      const xform &xf1, const xform &xf2,
                      const KDtree *kd1, const KDtree *kd2,
                      const Grid *g1, const Grid *g2,
                      vector<float> &o1, vector<float> &o2,
                      float &maxdist, int verbose)
{
	timestamp t = now();
	if (maxdist <= 0.0f) {
		if (g1 && g2) {
			maxdist = min(g1->bbox_size(), g2->bbox_size());
		} else {
			mesh1->need_bbox();
			mesh2->need_bbox();
			maxdist = min(len(mesh1->bbox.size()),
			              len(mesh2->bbox.size()));
		}
	}

	const size_t nv1 = mesh1->vertices.size(), nv2 = mesh2->vertices.size();
	o1.clear(); o1.resize(nv1);
	o2.clear(); o2.resize(nv2);

	xform xf12 = inv(xf2) * xf1;
	xform xf21 = inv(xf1) * xf2;
	float maxdist2 = sqr(maxdist);

#pragma omp parallel
	{
#pragma omp for nowait
		for (size_t i = 0; i < nv1; i++) {
			point p = xf12 * mesh1->vertices[i];
			if (g2 && !g2->overlaps(p))
				continue;
			if (kd2) {
				const float *match = kd2->closest_to_pt(p, maxdist2);
				if (!match)
					continue;
			}
			o1[i] = 1;
		}

#pragma omp for
		for (size_t i = 0; i < nv2; i++) {
			point p = xf21 * mesh2->vertices[i];
			if (g1 && !g1->overlaps(p))
				continue;
			if (kd1) {
				const float *match = kd1->closest_to_pt(p, maxdist2);
				if (!match)
					continue;
			}
			o2[i] = 1;
		}
	} // omp parallel

	if (verbose > 1) {
		dprintf("Computed overlaps in %.3f msec.\n\n",
			(now() - t) * 1000.0f);
	}
}


// Select a number of points and find correspondences
static void select_and_match(TriMesh *mesh1, TriMesh *mesh2,
                             const xform &xf1, const xform &xf2,
                             const KDtree *kd2,
                             const vector<float> &sampcdf1, float cdfincr,
                             float maxdist, bool do_flip,
                             vector<PtPair> &pairs)
{
	xform nxf1 = norm_xf(xf1);
	xform nxf2 = norm_xf(xf2);
	xform xf12 = inv(xf2) * xf1;
	xform nxf12 = norm_xf(xf12);
	bool is_pointcloud1 = (mesh1->faces.empty() &&
		mesh1->tstrips.empty() && mesh1->grid.empty());
	bool is_pointcloud2 = (mesh2->faces.empty() &&
		mesh2->tstrips.empty() && mesh2->grid.empty());
	bool is_pointcloud = (is_pointcloud1 || is_pointcloud2);
	float maxdist2 = sqr(maxdist);

	size_t i = 0, n = sampcdf1.size();
	float cdfval = uniform_rnd(cdfincr);
	while (cdfval < 1.0f) {
		if (sampcdf1[i] <= cdfval) {
			// Find next sample point on mesh1 using the CDF.
			// We're looking for the first point such that
			// its CDF value is > cdfval.  The algorithm is a
			// variant of binary search.  We start with a
			// step size of 1 and repeatedly double it until
			// we overshoot...
			size_t step = 1;
			while (i + step < n && sampcdf1[i+step] <= cdfval) {
				i = i + step;
				step <<= 1;
			}
			// ... then repeatedly halve step size until it is
			// back to 0.  Invariants: step is a power of 2 and
			// sampcdf1[i] <= cdfval and sampcdf1[i+step] > cdfval
			step >>= 1;
			while (step) {
				if (i + step < n && sampcdf1[i+step] <= cdfval)
					i = i + step;
				step >>= 1;
			}
			// i ended up being the *last* location at which
			// sampcdf1[i] <= cdfval, so we increment it to find
			// the *first* location where sampcdf1[i] > cdfval
			i++;
		}

		// Set up cdfval for the next iteration
		cdfval += cdfincr;

		// Transform into coords of mesh2 and match
		point p12 = xf12 * mesh1->vertices[i];

		const float *match;
		if (USE_NORMCOMPAT && !is_pointcloud) {
			vec n12 = nxf12 * mesh1->normals[i];
			NormCompat nc(n12, mesh2);
			match = kd2->closest_to_pt(p12, maxdist2, &nc);
		} else {
			match = kd2->closest_to_pt(p12, maxdist2);
		}
		if (!match)
			continue;

		// Reject boundary matches, if desired
		int match_ind = (const point *) match - &(mesh2->vertices[0]);
		if (REJECT_BDY && !is_pointcloud && mesh2->is_bdy(match_ind))
			continue;

		// Project both points into world coords and save
		point p1 = xf1 * mesh1->vertices[i];
		vec n1 = nxf1 * mesh1->normals[i];
		point p2 = xf2 * mesh2->vertices[match_ind];
		vec n2 = nxf2 * mesh2->normals[match_ind];
		if ((n1 DOT n2) < 0.0f)
			n2 = -n2;

		if (do_flip)
			pairs.push_back(PtPair(p2, n2, p1, n1));
		else
			pairs.push_back(PtPair(p1, n1, p2, n2));
	}
}


// Determinant of a 3x3 matrix
static float det(const float (&A)[3][3])
{
	return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
	       A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +
	       A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}


// Do rigid-body point-to-point alignment.  (This is done in the early stages
// of registration to improve stability.)
static void align_pt2pt(const vector<PtPair> &pairs,
			const point &centroid1, const point &centroid2,
			xform &alignxf)
{
	size_t n = pairs.size();

	float A[3][3] = { { 0 } };
	for (size_t i = 0; i < n; i++) {
		vec v1 = pairs[i].p1 - centroid1;
		vec v2 = pairs[i].p2 - centroid2;
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				A[j][k] += v1[j] * v2[k];
	}
	float s[3], V[3][3];
	svd<float,3,3>(A, s, V);
	if ((det(A) < 0.0f) ^ (det(V) < 0.0f)) {
		V[2][0] = -V[2][0];
		V[2][1] = -V[2][1];
		V[2][2] = -V[2][2];
	}
	alignxf = xform::trans(centroid1) *
	          xform::fromarray(A) * transp(xform::fromarray(V)) *
	          xform::trans(-centroid2);
}


// Do symmetric point-to-plane alignment, returning alignxf
// as well as eigenvectors and inverse eigenvalues
static void align_pt2pl(const vector<PtPair> &pairs, float scale,
                        const point &centroid1, const point &centroid2,
                        float (&evec)[6][6], float (&einv)[6],
                        xform &alignxf)
{
	size_t n = pairs.size();

	float b[6] = { 0 };
	for (size_t i = 0; i < n; i++) {
		vec p1 = pairs[i].p1 - centroid1;
		vec p2 = pairs[i].p2 - centroid2;
		const vec &n1 = pairs[i].n1;
		const vec &n2 = pairs[i].n2;
		vec nn = n1 + n2;
		vec c = scale * (p1 CROSS n2 + p2 CROSS n1);
		vec p12 = p1 - p2;
		float d = scale * (p12 DOT nn);

		float x[6] = { c[0], c[1], c[2], nn[0], nn[1], nn[2] };
		for (int j = 0; j < 6; j++) {
			b[j] += d * x[j];
			for (int k = j; k < 6; k++)
				evec[j][k] += x[j] * x[k];
		}

		// Regularization for rotational component - point to point
		float reg = REGULARIZATION * sqr(scale);
		evec[0][0] += reg * (sqr(p2[1]) + sqr(p2[2]));
		evec[0][1] -= reg * p2[0] * p2[1];
		evec[0][2] -= reg * p2[0] * p2[2];
		evec[1][1] += reg * (sqr(p2[0]) + sqr(p2[2]));
		evec[1][2] -= reg * p2[1] * p2[2];
		evec[2][2] += reg * (sqr(p2[0]) + sqr(p2[1]));
		vec c2 = p2 CROSS p12;
		x[0] += reg * c2[0];
		x[1] += reg * c2[1];
		x[2] += reg * c2[2];
	}

	// Regularization for translational component
	evec[3][3] += REGULARIZATION * n;
	evec[4][4] += REGULARIZATION * n;
	evec[5][5] += REGULARIZATION * n;

	// Make matrix symmetric
	for (int j = 0; j < 6; j++)
		for (int k = 0; k < j; k++)
			evec[j][k] = evec[k][j];

	// Eigen-decomposition and inverse
	float eval[6];
	eigdc<float,6>(evec, eval);
	for (int i = 0; i < 6; i++)
		einv[i] = 1.0f / eval[i];

	// Solve system
	eigmult<float,6>(evec, einv, b);

	// Extract rotation and translation
	vec rot(b[0], b[1], b[2]), trans(b[3], b[4], b[5]);
	float rotangle = atan(len(rot));
	float cosangle = cos(rotangle);
	trans *= cosangle / ((1.0f + cosangle) * scale);

	alignxf = xform::trans(trans + centroid1) *
	          xform::rot(rotangle, rot) *
	          xform::trans(trans - centroid2);
}


// Do symmetric point-to-plane translation-only alignment
static void align_pt2pl_trans(const vector<PtPair> &pairs,
                              const point &centroid1, const point &centroid2,
                              xform &alignxf)
{
	size_t n = pairs.size();

	float evec[3][3] = { { 0 } }, einv[3] = { 0 };
	vec b;
	for (size_t i = 0; i < n; i++) {
		vec p1 = pairs[i].p1 - centroid1;
		vec p2 = pairs[i].p2 - centroid2;
		const vec &n1 = pairs[i].n1;
		const vec &n2 = pairs[i].n2;
		vec nn = n1 + n2;
		float d = (p1 - p2) DOT nn;

		for (int j = 0; j < 3; j++) {
			b[j] += d * nn[j];
			for (int k = 0; k < 3; k++)
				evec[j][k] += nn[j] * nn[k];
		}
	}

	// Regularization
	evec[0][0] += REGULARIZATION * n;
	evec[1][1] += REGULARIZATION * n;
	evec[2][2] += REGULARIZATION * n;

	// Eigen-decomposition and inverse
	vec eval;
	eigdc<float,3>(evec, eval);
	for (int i = 0; i < 3; i++)
		einv[i] = 1.0f / eval[i];

	// Solve system
	eigmult<float,3>(evec, einv, b);
	b += centroid1 - centroid2;
	alignxf = xform::trans(b);
}


// Compute isotropic or anisotropic scale.  Assumes alignxf already contains
// a rigid-body transformation to be applied to pairs[i].p2
static void align_scale(const vector<PtPair> &pairs, xform &alignxf,
                        const point &centroid1, const point &centroid2,
                        bool do_affine)
{
	size_t n = pairs.size();

	point centroid = 0.5f * (centroid1 + alignxf * centroid2);

	// Compute covariance matrices
	float cov1[3][3] = { { 0 } };
	float cov2[3][3] = { { 0 } };
	for (size_t i = 0; i < n; i++) {
		point p1 = pairs[i].p1 - centroid;
		point p2 = alignxf * pairs[i].p2 - centroid;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				cov1[j][k] += p1[j] * p1[k];
				cov2[j][k] += p2[j] * p2[k];
			}
		}
	}

	// Compute eigenstuff of cov
	vec eval1, eval2;
	eigdc<float,3>(cov1, eval1);
	eigdc<float,3>(cov2, eval2);

	if (!do_affine) {
		// Just uniform scale
		alignxf = xform::trans(centroid) *
		          xform::scale(sqrt(eval1.sum() / eval2.sum())) *
		          xform::trans(-centroid) *
		          alignxf;
		return;
	}

	// Compute sqrt of covariance
	float csqrt1[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
	float icsqrt2[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
	for (int i = 0; i < 3; i++) {
		eigmult<float,3>(cov1, sqrt(eval1), csqrt1[i]);
		eigmult<float,3>(cov2, sqrt(1.0f / eval2), icsqrt2[i]);
	}

	alignxf = xform::trans(centroid) *
	          xform::fromarray(csqrt1) *
	          xform::fromarray(icsqrt2) *
	          xform::trans(-centroid) *
	          alignxf;
}


// Compute point-to-point or point-to-plane squared distances
static void compute_dist2(const vector<PtPair> &pairs,
	vector<float> &distances2, ICP_iter_type iter_type)
{
	size_t np = pairs.size();
	distances2.resize(np);
	if (iter_type == ICP_POINT_TO_POINT) {
		for (size_t i = 0; i < np; i++)
			distances2[i] = dist2(pairs[i].p1, pairs[i].p2);
	} else /* ICP_POINT_TO_PLANE */ {
		float p2pl_mult = 0.5f / (1.0f + REGULARIZATION);
		float p2pt_mult = REGULARIZATION / (1.0f + REGULARIZATION);
		for (size_t i = 0; i < np; i++) {
			vec v = pairs[i].p1 - pairs[i].p2;
			distances2[i] = p2pl_mult * (sqr(v DOT pairs[i].n1) +
				sqr(v DOT pairs[i].n2));
			distances2[i] += p2pt_mult * len2(v);
		}
	}
}


// Find the median of a list of numbers - makes a copy of vals (since it's
// passed by value, not by reference) so that we can modify it
static float median(vector<float> vals)
{
	size_t n = vals.size();
	if (!n)
		return 0.0f;

	size_t mid = n / 2;
	nth_element(vals.begin(), vals.begin() + mid, vals.end());
	return vals[mid];
}


// Find the mean of a list of numbers
static float mean(const vector<float> &vals)
{
	size_t n = vals.size();
	if (!n)
		return 0.0f;

	float sum = 0.0f;
	for (size_t i = 0; i < n; i++)
		sum += vals[i];

	return sum / n;
}


// Do one iteration of ICP
static float ICP_iter(TriMesh *mesh1, TriMesh *mesh2,
                      const xform &xf1, xform &xf2,
                      const KDtree *kd1, const KDtree *kd2,
                      const vector<float> &weights1, const vector<float> &weights2,
                      vector<float> &sampcdf1, vector<float> &sampcdf2,
                      float cdfincr, bool update_cdfs,
                      float dist_thresh_mult, float normdot_thresh,
                      float &maxdist, int verbose,
                      ICP_iter_type iter_type, ICP_xform_type xform_type)
{
	const size_t nv1 = mesh1->vertices.size(), nv2 = mesh2->vertices.size();

	// Compute pairs
	timestamp t1 = now();
	if (verbose > 1) {
		dprintf("Selecting points using cdfincr = %f\n", cdfincr);
		dprintf("Matching with maxdist = %g\n", maxdist);
	}

	vector<PtPair> pairs;
	pairs.reserve(DESIRED_PAIRS);
	select_and_match(mesh1, mesh2, xf1, xf2, kd2, sampcdf1, cdfincr,
		maxdist, false, pairs);
	select_and_match(mesh2, mesh1, xf2, xf1, kd1, sampcdf2, cdfincr,
		maxdist, true, pairs);

	timestamp t2 = now();
	size_t npairs = pairs.size();
	if (verbose > 1) {
		dprintf("Generated %lu pairs in %.3f msec.\n",
			(unsigned long) npairs, (t2 - t1) * 1000.0f);
	}

	// Compute median point-to-point or point-to-plane distance.
	vector<float> distances2;
	compute_dist2(pairs, distances2, iter_type);
	float median_dist = sqrt(median(distances2));

	// Now compute threshold for rejection, as a multiple of sigma
	// (estimated robustly based on the median).
	float sigma = MEDIAN_TO_SIGMA * median_dist;
	float dist_thresh = dist_thresh_mult * sigma;
	float dist_thresh2 = sqr(dist_thresh);
	
	// We also compute the new maxdist, but that always has to be
	// based on point-to-point distances (since that's how the KDtree
	// interprets maxdist).
	if (!update_cdfs) {
		if (iter_type == ICP_POINT_TO_POINT) {
			maxdist = dist_thresh_mult * sigma;
		} else {
			vector<float> distances2_pt2pt;
			compute_dist2(pairs, distances2_pt2pt, ICP_POINT_TO_POINT);
			float sigma_pt2pt = MEDIAN_TO_SIGMA * sqrt(median(distances2_pt2pt));
			maxdist = dist_thresh_mult * sigma_pt2pt;
		}
	}


	// Reject
	if (verbose > 1)
		dprintf("Rejecting pairs with dist > %g or angle > %.1f\n",
			dist_thresh, degrees(acos(normdot_thresh)));
	float err = 0.0f;
	size_t next = 0;
	for (size_t i = 0; i < npairs; i++) {
		if (distances2[i] > dist_thresh2)
			continue;
		if ((pairs[i].n1 DOT pairs[i].n2) < normdot_thresh)
			continue;
		pairs[next++] = pairs[i];
		err += distances2[i];
	}
	pairs.erase(pairs.begin() + next, pairs.end());

	timestamp t3 = now();
	if (verbose > 1) {
		dprintf("Rejected %lu pairs in %.3f msec.\n",
			(unsigned long) (npairs - pairs.size()),
			(t3 - t2) * 1000.0f);
	}

	npairs = pairs.size();
	if (npairs < MIN_PAIRS) {
		if (verbose)
			dprintf("Too few point pairs.\n");
		return -1.0f;
	}

	err = sqrt(err / npairs);
	const char *dist_type = (iter_type == ICP_POINT_TO_POINT) ?
		"point-to-point" : "point-to-plane";
	if (verbose > 1)
		dprintf("RMS %s error before alignment = %g\n", dist_type, err);

	// Compute centroids and scale
	point centroid1, centroid2;
	for (size_t i = 0; i < npairs; i++) {
		centroid1 += pairs[i].p1;
		centroid2 += pairs[i].p2;
	}
	centroid1 /= npairs;
	centroid2 /= npairs;

	float scale = 0.0f;
	for (size_t i = 0; i < npairs; i++) {
		scale += dist2(pairs[i].p1, centroid1);
		scale += dist2(pairs[i].p2, centroid2);
	}
	scale = sqrt(scale / (2 * npairs));
	scale = 1.0f / scale;

	// Do the minimization
	float evec[6][6] = { { 0 } }, einv[6] = { 0 };
	xform alignxf;
	if (xform_type == ICP_TRANSLATION) {
		if (iter_type == ICP_POINT_TO_POINT) {
			alignxf = xform::trans(centroid1 - centroid2);
		} else {
			align_pt2pl_trans(pairs, centroid1, centroid2, alignxf);
		}
		update_cdfs = false;
	} else {
		// First do rigid-body alignment
		if (iter_type == ICP_POINT_TO_POINT) {
			align_pt2pt(pairs, centroid1, centroid2, alignxf);
			update_cdfs = false;
		} else {
			align_pt2pl(pairs, scale, centroid1, centroid2,
			            evec, einv, alignxf);
		}
		// ... and then estimate the scale on top of that, if required
		if (xform_type == ICP_SIMILARITY)
			align_scale(pairs, alignxf, centroid1, centroid2, false);
		else if (xform_type == ICP_AFFINE)
			align_scale(pairs, alignxf, centroid1, centroid2, true);
	}

	// Apply transform, and find distance after alignment
	xf2 = alignxf * xf2;
	if (xform_type == ICP_RIGID)
		orthogonalize(xf2);
	xform nalignxf = norm_xf(alignxf);
	for (size_t i = 0; i < npairs; i++) {
		pairs[i].p2 = alignxf * pairs[i].p2;
		pairs[i].n2 = nalignxf * pairs[i].n2;
	}
	compute_dist2(pairs, distances2, iter_type);
	err = sqrt(mean(distances2));

	timestamp t4 = now();
	if (verbose > 1) {
		dprintf("Computed xform in %.3f msec.\n",
			(t4 - t3) * 1000.0f);
		dprintf("RMS %s error after alignment = %g\n\n", dist_type, err);
	}

	// See if we need to update CDFs
	if (!update_cdfs)
		return err;

	centroid2 = alignxf * centroid2;
	xform nxf1 = norm_xf(xf1), nxf2 = norm_xf(xf2);
	double sum_sampcdf1 = 0, sum_sampcdf2 = 0;
#pragma omp parallel
	{
#pragma omp for nowait reduction(+ : sum_sampcdf1)
		for (size_t i = 0; i < nv1; i++) {
			if (!weights1[i]) {
				sampcdf1[i] = 0.0;
				continue;
			}
			point p = scale * (xf1 * mesh1->vertices[i] - centroid1);
			vec n = nxf1 * mesh1->normals[i];
			vec c = p CROSS n;
			float s = 0.0f;
			for (int j = 0; j < 6; j++) {
				float tmp = evec[0][j] * c[0] + evec[1][j] * c[1] +
				            evec[2][j] * c[2] + evec[3][j] * n[0] +
				            evec[4][j] * n[1] + evec[5][j] * n[2];
				// Compromise between uniform and inverse
				// eigenvalue-directed sampling - sqrt
				s += sqrt(einv[j]) * sqr(tmp);
			}
			sum_sampcdf1 += (sampcdf1[i] = s * weights1[i]);
		}

#pragma omp for reduction(+ : sum_sampcdf2)
		for (size_t i = 0; i < nv2; i++) {
			if (!weights2[i]) {
				sampcdf2[i] = 0.0;
				continue;
			}
			point p = scale * (xf2 * mesh2->vertices[i] - centroid2);
			vec n = nxf2 * mesh2->normals[i];
			vec c = p CROSS n;
			float s = 0.0f;
			for (int j = 0; j < 6; j++) {
				float tmp = evec[0][j] * c[0] + evec[1][j] * c[1] +
				            evec[2][j] * c[2] + evec[3][j] * n[0] +
				            evec[4][j] * n[1] + evec[5][j] * n[2];
				s += sqrt(einv[j]) * sqr(tmp);
			}
			sum_sampcdf2 += (sampcdf2[i] = s * weights2[i]);
		}
	} // omp parallel

	if (!sum_sampcdf1 || !sum_sampcdf2) {
		if (verbose)
			dprintf("No overlap.\n");
		return -1.0f;
	}

	float cdf_scale1 = 1 / sum_sampcdf1;
	sampcdf1[0] *= cdf_scale1;
	for (size_t i = 1; i < nv1 - 1; i++)
		sampcdf1[i] = cdf_scale1 * sampcdf1[i] + sampcdf1[i-1];
	sampcdf1[nv1-1] = 1.0f;

	float cdf_scale2 = 1 / sum_sampcdf2;
	sampcdf2[0] *= cdf_scale2;
	for (size_t i = 1; i < nv2 - 1; i++)
		sampcdf2[i] = cdf_scale2 * sampcdf2[i] + sampcdf2[i-1];
	sampcdf2[nv2-1] = 1.0f;

	timestamp t5 = now();
	if (verbose > 1) {
		dprintf("Updated CDFs in %.3f msec.\n\n",
			(t5 - t4) * 1000.0f);
	}

	return err;
}


// Create a CDF for simple weighted sampling
static void make_uniform_cdfs(
	const vector<float> &weights1, vector<float> &sampcdf1,
	const vector<float> &weights2, vector<float> &sampcdf2)
{
	const size_t nv1 = weights1.size(), nv2 = weights2.size();
	sampcdf1.resize(nv1);
	sampcdf2.resize(nv2);

	double sum_sampcdf1 = 0, sum_sampcdf2 = 0;
#pragma omp parallel
	{
#pragma omp for nowait reduction(+ : sum_sampcdf1)
		for (size_t i = 0; i < nv1; i++)
			sum_sampcdf1 += (sampcdf1[i] = weights1[i]);
#pragma omp for reduction(+ : sum_sampcdf2)
		for (size_t i = 0; i < nv2; i++)
			sum_sampcdf2 += (sampcdf2[i] = weights2[i]);
	}

	float cdf_scale1 = 1 / sum_sampcdf1;
	sampcdf1[0] *= cdf_scale1;
	for (size_t i = 1; i < nv1 - 1; i++)
		sampcdf1[i] = cdf_scale1 * sampcdf1[i] + sampcdf1[i-1];
	sampcdf1[nv1-1] = 1.0f;

	float cdf_scale2 = 1 / sum_sampcdf2;
	sampcdf2[0] *= cdf_scale2;
	for (size_t i = 1; i < nv2 - 1; i++)
		sampcdf2[i] = cdf_scale2 * sampcdf2[i] + sampcdf2[i-1];
	sampcdf2[nv2-1] = 1.0f;
}


// Do ICP.  Aligns mesh mesh2 to mesh1, updating xf2 with the new transform.
// Returns alignment error, or -1 on failure
float ICP(TriMesh *mesh1, TriMesh *mesh2,
          const xform &xf1, xform &xf2,
          const KDtree *kd1, const KDtree *kd2,
          vector<float> &weights1, vector<float> &weights2,
          float maxdist /* = 0.0f */, int verbose /* = 0 */,
          ICP_xform_type xform_type /* = ICP_RIGID */ )
{
	timestamp t_start = now();

	// Precompute normals, connectivity (used to determine boundaries),
	// and grids (used for fast overlap computation)
	Grid *g1, *g2;
#pragma omp parallel for
	for (int i = 0; i < 2; i++) {
		TriMesh *mesh = (i == 0) ? mesh1 : mesh2;
		Grid * &g = (i == 0) ? g1 : g2;

		mesh->need_normals();
		if (REJECT_BDY) {
			mesh->need_faces();
			mesh->need_neighbors();
			mesh->need_adjacentfaces();
		}
		g = new Grid(mesh->vertices);
	}

	// Initial maxdist, thresholds
	if (maxdist <= 0.0f)
		maxdist = min(g1->bbox_size(), g2->bbox_size());
	float dist_thresh_mult = DIST_THRESH_MULT_START;
	float normdot_thresh = NORMDOT_THRESH_START;

	// Weights and initial (uniform) CDFs
	const size_t nv1 = mesh1->vertices.size(), nv2 = mesh2->vertices.size();
	bool had_weights = (weights1.size() == nv1 && weights2.size() == nv2);
	weights1.resize(nv1, 1.0f);
	weights2.resize(nv2, 1.0f);

	vector<float> sampcdf1, sampcdf2;
	make_uniform_cdfs(weights1, sampcdf1, weights2, sampcdf2);
	float cdfincr = 2.0f / DESIRED_PAIRS;

	timestamp t_initial_iters = now();
	if (verbose > 1) {
		dprintf("\nTime for preprocessing: %.3f msec.\n\n",
			(t_initial_iters - t_start) * 1000.0f);
	}

	// First, do a few point-to-point iterations for stability in case the
	// initial misalignment is big.
	ICP_xform_type iter_xform_type =
		(xform_type == ICP_TRANSLATION) ? ICP_TRANSLATION : ICP_RIGID;
	float err;
	for (int iter = 0; iter < INITIAL_ITERS; iter++) {
		err = ICP_iter(mesh1, mesh2, xf1, xf2, kd1, kd2,
		               weights1, weights2, sampcdf1, sampcdf2,
		               cdfincr, false, dist_thresh_mult,
		               normdot_thresh, maxdist, verbose,
		               ICP_POINT_TO_POINT, iter_xform_type);
		if (err < 0) {
			delete g1;
			delete g2;
			if (!had_weights) {
				weights1.clear();
				weights2.clear();
			}
			return err;
		}
	}

	timestamp t_main_iters = now();
	if (verbose > 1) {
		dprintf("Time for initial iterations: %.3f msec.\n\n",
		       (t_main_iters - t_initial_iters) * 1000.0f);
	}

	// Now the real (point-to-plane) iterations
	float min_err = 0.0f;
	int iter_of_min_err = -1;
	int iter;
	for (iter = 0; iter < MAX_ITERS; iter++) {
		// Update thresholds
		dist_thresh_mult = mix(dist_thresh_mult, DIST_THRESH_MULT_FINAL,
			THRESH_RATE_OF_CHANGE);
		normdot_thresh = mix(normdot_thresh, NORMDOT_THRESH_FINAL,
			THRESH_RATE_OF_CHANGE);

		// Should we recompute overlaps and CDFs?
		bool recompute = ((iter % CDF_UPDATE_INTERVAL) == 0);
		if (!had_weights && recompute && iter != 0)
			compute_overlaps(mesh1, mesh2, xf1, xf2,
			                 NULL, NULL, g1, g2,
			                 weights1, weights2,
			                 maxdist, verbose);

		// If we're recomputing CDFs, use uniform sampling on this
		// iteration to make sure that covariance is unbiased
		if (recompute && iter != 0)
			make_uniform_cdfs(weights1, sampcdf1,
			                  weights2, sampcdf2);

		// Do an iteration
		err = ICP_iter(mesh1, mesh2, xf1, xf2, kd1, kd2,
		               weights1, weights2, sampcdf1, sampcdf2,
		               cdfincr, recompute, dist_thresh_mult,
		               normdot_thresh, maxdist, verbose,
		               ICP_POINT_TO_PLANE, iter_xform_type);
		if (err < 0) {
			delete g1;
			delete g2;
			if (!had_weights) {
				weights1.clear();
				weights2.clear();
			}
			return err;
		}

		if ((err < min_err || iter_of_min_err < 0) && !recompute) {
			min_err = err;
			iter_of_min_err = iter;
		}

		// Stop if we've gone at least TERMINATION_ITER_THRESH
		// iterations without seeing a new minimum error
		if (iter - iter_of_min_err >= TERMINATION_ITER_THRESH &&
		    iter_of_min_err >= 0 && !recompute &&
		    xform_type != ICP_SIMILARITY && xform_type != ICP_AFFINE) {
			iter++; // Get #-of-iters printf correct
			break;
		}

		// If we're optimizing for similarity or affine, switch on
		// those transformations after MAX_ITERS/2 iterations
		if (iter == MAX_ITERS / 2 &&
		    (xform_type == ICP_SIMILARITY || xform_type == ICP_AFFINE))
			iter_xform_type = xform_type;
	}

	// Some final iterations at a higher sampling rate...
	if (verbose > 1) {
		dprintf("Time for %d iterations: %.3f msec.\n\n",
			iter, (now() - t_main_iters) * 1000.0f);
	}

	for (iter = 0; iter < FINAL_ITERS; iter++) {
		if (iter == 0) {
			cdfincr *= (float) DESIRED_PAIRS / DESIRED_PAIRS_FINAL;
			// Use uniform sampling so that the final error
			// we return is unbiased
			make_uniform_cdfs(weights1, sampcdf1, weights2, sampcdf2);
		}
		err = ICP_iter(mesh1, mesh2, xf1, xf2, kd1, kd2,
		               weights1, weights2, sampcdf1, sampcdf2,
		               cdfincr, false, DIST_THRESH_MULT_FINAL,
		               NORMDOT_THRESH_FINAL, maxdist, verbose,
		               ICP_POINT_TO_PLANE, iter_xform_type);
		if (err < 0) {
			delete g1;
			delete g2;
			if (!had_weights) {
				weights1.clear();
				weights2.clear();
			}
			return err;
		}
	}
	if (verbose == 1) {
		dprintf("ICP error = %g\n", err);
	} else if (verbose > 1) {
		// err already printed out in ICP_iter
		dprintf("Time for ICP: %.3f msec.\n\n",
		       (now() - t_start) * 1000.0f);
	}

	delete g1;
	delete g2;
	if (!had_weights) {
		weights1.clear();
		weights2.clear();
	}
	return err;
}


// Easier-to-use interfaces to ICP
float ICP(TriMesh *mesh1, TriMesh *mesh2,
          const xform &xf1, xform &xf2,
          const KDtree *kd1, const KDtree *kd2,
          int verbose /* = 0 */, ICP_xform_type xform_type /* = ICP_RIGID */)
{
	vector<float> weights1, weights2;
	return ICP(mesh1, mesh2, xf1, xf2, kd1, kd2,
	           weights1, weights2, 0.0f, verbose, xform_type);
}

float ICP(TriMesh *mesh1, TriMesh *mesh2,
          const xform &xf1, xform &xf2,
          int verbose /* = 0 */, ICP_xform_type xform_type /* = ICP_RIGID */)
{
	if (verbose > 1)
		dprintf("\nBuilding KDtrees... ");
	timestamp t = now();

	KDtree *kd1, *kd2;
#pragma omp parallel for
	for (int i = 0; i < 2; i++) {
		if (i == 0)
			kd1 = new KDtree(mesh1->vertices);
		else
			kd2 = new KDtree(mesh2->vertices);
	}

	if (verbose > 1)
		dprintf("Done.  %.3f msec.\n", (now() - t) * 1000.0f);

	vector<float> weights1, weights2;
	float icperr = ICP(mesh1, mesh2, xf1, xf2, kd1, kd2,
	                   weights1, weights2, 0.0f, verbose, xform_type);
	delete kd2;
	delete kd1;
	return icperr;
}


// Compatibility interfaces to ICP from before we had ICP_xform_type
float ICP(TriMesh *mesh1, TriMesh *mesh2,
          const xform &xf1, xform &xf2,
          const KDtree *kd1, const KDtree *kd2,
          ::std::vector<float> &weights1, ::std::vector<float> &weights2,
          float maxdist, int verbose, bool do_scale, bool do_affine)
{
	ICP_xform_type xform_type = do_affine ? ICP_AFFINE :
		do_scale ? ICP_SIMILARITY: ICP_RIGID;
	return ICP(mesh1, mesh2, xf1, xf2, kd1, kd2,
	           weights1, weights2, maxdist, verbose, xform_type);
}

float ICP(TriMesh *mesh1, TriMesh *mesh2, const xform &xf1, xform &xf2,
          const KDtree *kd1, const KDtree *kd2,
          int verbose, bool do_scale, bool do_affine)
{
	ICP_xform_type xform_type = do_affine ? ICP_AFFINE :
		do_scale ? ICP_SIMILARITY: ICP_RIGID;
	return ICP(mesh1, mesh2, xf1, xf2, kd1, kd2, verbose, xform_type);
}

float ICP(TriMesh *mesh1, TriMesh *mesh2, const xform &xf1, xform &xf2,
          int verbose, bool do_scale, bool do_affine)
{
	ICP_xform_type xform_type = do_affine ? ICP_AFFINE :
		do_scale ? ICP_SIMILARITY: ICP_RIGID;
	return ICP(mesh1, mesh2, xf1, xf2, verbose, xform_type);
}


} // namespace trimesh
