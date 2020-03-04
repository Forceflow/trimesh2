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


#define MAX_ITERS 150
#define TERMINATION_ITER_THRESH 25
#define FINAL_ITERS 2
#define MIN_PAIRS 10
#define DESIRED_PAIRS 200
#define DESIRED_PAIRS_FINAL 5000
#define CDF_UPDATE_INTERVAL 20
#define APPROX_EPS 0.05f
#define REJECT_BDY false
#define USE_NORMCOMPAT true
#define HUBER_THRESH_MULT 2.0f
#define REGULARIZATION 0.0002f
#define DIST_THRESH_MULT 12.0f
#define DIST_THRESH_MULT_FINAL 4.0f
#define ANGLE_THRESH_MULT 1.5f
#define ANGLE_THRESH_MIN 0.1f
#define ANGLE_THRESH_MAX 1.0f
#define dprintf TriMesh::dprintf


namespace trimesh {


// A pair of points, with associated normals
struct PtPair {
	vec p1, n1, p2, n2;
	PtPair(const point &p1_, const vec &n1_,
	       const point &p2_, const vec &n2_) :
			p1(p1_), n1(n1_), p2(p2_), n2(n2_)
		{}
};


// A class for evaluating compatibility of normals during KDtree searches
class NormCompat : public KDtree::CompatFunc {
private:
	TriMesh *m;
	const vec &n;
	float thresh;

public:
	NormCompat(const vec &n_, TriMesh *m_, float thresh_):
		m(m_), n(n_), thresh(thresh_)
		{}
	virtual bool operator () (const float *p) const
	{
		int idx = (const point *) p - &(m->vertices[0]);
		return (n DOT m->normals[idx]) > thresh;
	}
};


// Determine which points on mesh1 and mesh2 overlap the other,
// filling in o1 and o2
void compute_overlaps(TriMesh *mesh1, TriMesh *mesh2,
                      const xform &xf1, const xform &xf2,
                      const KDtree *kd1, const KDtree *kd2,
                      vector<float> &o1, vector<float> &o2,
                      float maxdist, int verbose)
{
	timestamp t = now();

	const int nv1 = mesh1->vertices.size(), nv2 = mesh2->vertices.size();
	o1.clear(); o1.resize(nv1);
	o2.clear(); o2.resize(nv2);

	xform xf12 = inv(xf2) * xf1;
	xform xf21 = inv(xf1) * xf2;

#pragma omp parallel
	{
#pragma omp for nowait
		for (int i = 0; i < nv1; i++) {
			point p = xf12 * mesh1->vertices[i];
			if (kd2->exists_pt_within(p, maxdist))
				o1[i] = 1;
		}

#pragma omp for
		for (int i = 0; i < nv2; i++) {
			point p = xf21 * mesh2->vertices[i];
			if (kd1->exists_pt_within(p, maxdist))
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
                             float maxdist, float angle_thresh, bool do_flip,
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
	float normdot_thresh = cos(angle_thresh);

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
			NormCompat nc(n12, mesh2, normdot_thresh);
			match = kd2->closest_to_pt(p12, maxdist2, &nc, APPROX_EPS);
		} else {
			match = kd2->closest_to_pt(p12, maxdist2, APPROX_EPS);
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


// Do symmetric point-to-plane alignment, returning alignxf
// as well as eigenvectors and inverse eigenvalues
static void align_symm(const vector<PtPair> &pairs, float scale,
                       const point &centroid1, const point &centroid2,
                       float median_dist, xform &alignxf)
{
	float huber_thresh = HUBER_THRESH_MULT * scale * median_dist;
	size_t npairs = pairs.size();
	float A[6][6] = { { 0 } }, b[6] = { 0 };
	for (size_t i = 0; i < npairs; i++) {
		vec p1 = scale * (pairs[i].p1 - centroid1);
		vec p2 = scale * (pairs[i].p2 - centroid2);
		vec n = pairs[i].n1 + pairs[i].n2;
		vec p = p1 + p2;
		vec c = p CROSS n;
		vec d = p1 - p2;

		float x[6] = { c[0], c[1], c[2], n[0], n[1], n[2] };
		float dn = d DOT n;

		// Huber weights, used for IRLS
		float wt = huber_thresh / max(fabs(dn), huber_thresh);

		for (int j = 0; j < 6; j++) {
			b[j] += wt * dn * x[j];
			for (int k = j; k < 6; k++)
				A[j][k] += wt * x[j] * x[k];
		}
	}

	// Make matrix symmetric
	for (int j = 1; j < 6; j++)
		for (int k = 0; k < j; k++)
			A[j][k] = A[k][j];

	// Eigen-decomposition and inverse
	float eval[6], einv[6];
	eigdc<float,6>(A, eval);
	for (int i = 0; i < 6; i++)
		einv[i] = 1.0f / (eval[i] + REGULARIZATION * eval[5]);

	// Solve system
	eigmult<float,6>(A, einv, b);

	// Extract rotation and translation
	vec rot(b[0], b[1], b[2]), trans(b[3], b[4], b[5]);
	float rotangle = atan(len(rot));
	trans *= cos(rotangle);
	trans *= 1.0f / scale;

	xform R = xform::rot(rotangle, rot);
	alignxf = xform::trans(centroid1) *
	          R * xform::trans(trans) * R *
	          xform::trans(-centroid2);
}


// Do symmetric point-to-plane translation-only alignment
static void align_pt2pl_trans(const vector<PtPair> &pairs,
                              const point &centroid1, const point &centroid2,
                              xform &alignxf)
{
	size_t npairs = pairs.size();

	float evec[3][3] = { { 0 } }, einv[3] = { 0 };
	vec b;
	for (size_t i = 0; i < npairs; i++) {
		vec p1 = pairs[i].p1 - centroid1;
		vec p2 = pairs[i].p2 - centroid2;
		vec n = 0.5f * (pairs[i].n1 + pairs[i].n2);
		float d = (p1 - p2) DOT n;

		for (int j = 0; j < 3; j++) {
			b[j] += d * n[j];
			for (int k = 0; k < 3; k++)
				evec[j][k] += n[j] * n[k];
		}
	}

	// Eigen-decomposition and inverse
	vec eval;
	eigdc<float,3>(evec, eval);
	for (int i = 0; i < 3; i++)
		einv[i] = 1.0f / (eval[i] + REGULARIZATION * eval[2]);

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
	size_t npairs = pairs.size();

	point centroid = 0.5f * (centroid1 + alignxf * centroid2);

	// Compute covariance matrices
	float cov1[3][3] = { { 0 } };
	float cov2[3][3] = { { 0 } };
	for (size_t i = 0; i < npairs; i++) {
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


// Do one iteration of ICP
static float ICP_iter(TriMesh *mesh1, TriMesh *mesh2,
                      const xform &xf1, xform &xf2,
                      const KDtree *kd1, const KDtree *kd2,
                      const vector<float> &weights1, const vector<float> &weights2,
                      vector<float> &sampcdf1, vector<float> &sampcdf2,
                      int desired_pairs, float &cdfincr, bool update_cdfs,
                      float &maxdist, float &angle_thresh,
		      int verbose, ICP_xform_type xform_type)
{
	const int nv1 = mesh1->vertices.size(), nv2 = mesh2->vertices.size();

	// Compute pairs
	timestamp t1 = now();
	if (verbose > 1) {
		dprintf("Selecting points using cdfincr = %f\n", cdfincr);
		dprintf("Matching with maxdist = %g\n", maxdist);
	}

	vector<PtPair> pairs;
	pairs.reserve(desired_pairs);
	select_and_match(mesh1, mesh2, xf1, xf2, kd2, sampcdf1, cdfincr,
		maxdist, angle_thresh, false, pairs);
	select_and_match(mesh2, mesh1, xf2, xf1, kd1, sampcdf2, cdfincr,
		maxdist, angle_thresh, true, pairs);

	timestamp t2 = now();
	size_t npairs = pairs.size();
	if (verbose > 1) {
		dprintf("Generated %lu pairs in %.3f msec.\n",
			(unsigned long) npairs, (t2 - t1) * 1000.0f);
	}

	// Compute median point-to-point distance and angle
	vector<float> distances2(npairs), normdots(npairs);
	for (size_t i = 0; i < npairs; i++) {
		distances2[i] = dist2(pairs[i].p1, pairs[i].p2);
		normdots[i] = pairs[i].n1 DOT pairs[i].n2;
	}
	float median_dist = sqrt(median(distances2));
	float median_angle = acos(median(normdots));

	// Compute rejection thresholds, which will also serve as
	// thresholds for the next iteration
	if (median_dist)
		maxdist = DIST_THRESH_MULT * median_dist;
	float dist2_thresh = sqr(maxdist);
	angle_thresh = clamp(ANGLE_THRESH_MULT * median_angle,
		ANGLE_THRESH_MIN, ANGLE_THRESH_MAX);
	float normdot_thresh = cos(angle_thresh);

	// Reject
	if (verbose > 1)
		dprintf("Rejecting pairs with dist > %g or angle > %.1f\n",
			maxdist, degrees(angle_thresh));
	size_t next = 0;
	for (size_t i = 0; i < npairs; i++) {
		if (distances2[i] > dist2_thresh ||
		    normdots[i] < normdot_thresh)
			continue;
		pairs[next++] = pairs[i];
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
	cdfincr *= (float) npairs / desired_pairs;

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
	xform alignxf;
	if (xform_type == ICP_TRANSLATION) {
		align_pt2pl_trans(pairs, centroid1, centroid2, alignxf);
	} else {
		// First do rigid-body alignment
		align_symm(pairs, scale, centroid1, centroid2, median_dist, alignxf);
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

	float err = 0.0f;
	for (size_t i = 0; i < npairs; i++) {
		vec n = normalized(pairs[i].n1 + nalignxf * pairs[i].n2);
		err += sqr((pairs[i].p1 - pairs[i].p2) DOT n);
	}
	err = sqrt(err / npairs);

	timestamp t4 = now();
	if (verbose > 1) {
		dprintf("Computed xform in %.3f msec.\n",
			(t4 - t3) * 1000.0f);
		dprintf("RMS error after alignment = %g\n\n", err);
	}

	// See if we need to update CDFs
	if (!update_cdfs)
		return err;

	// Construct CDFs based on inverse covariance of normals
	float A1[3][3] = { { 0 } }, A2[3][3] = { { 0 } };
	for (int i = 0; i < nv1; i++) {
		if (!weights1[i])
			continue;
		for (int j = 0; j < 3; j++)
			for (int k = j; k < 3; k++)
				A1[j][k] += mesh1->normals[i][j] *
					    mesh1->normals[i][k];
	}
	for (int i = 0; i < nv2; i++) {
		if (!weights2[i])
			continue;
		for (int j = 0; j < 3; j++)
			for (int k = j; k < 3; k++)
				A2[j][k] += mesh2->normals[i][j] *
					    mesh2->normals[i][k];
	}

	for (int j = 1; j < 3; j++) {
		for (int k = 0; k < j; k++) {
			A1[j][k] = A1[k][j];
			A2[j][k] = A2[k][j];
		}
	}

	float eval1[3], eval2[3];
	eigdc<float,3>(A1, eval1);
	eigdc<float,3>(A2, eval2);
	for (int i = 0; i < 3; i++) {
		eval1[i] = 1.0f / (eval1[i] + REGULARIZATION * eval1[2]);
		eval2[i] = 1.0f / (eval2[i] + REGULARIZATION * eval2[2]);
	}

	double sum_sampcdf1 = 0, sum_sampcdf2 = 0;
	for (int i = 0; i < nv1; i++) {
		if (!weights1[i]) {
			sampcdf1[i] = 0.0;
			continue;
		}
		float s = 0.0f;
		for (int j = 0; j < 3; j++)
			s += eval1[j] * sqr(
				A1[0][j] * mesh1->normals[i][0] +
				A1[1][j] * mesh1->normals[i][1] +
				A1[2][j] * mesh1->normals[i][2]);
		sum_sampcdf1 += (sampcdf1[i] = s * weights1[i]);
	}

	for (int i = 0; i < nv2; i++) {
		if (!weights2[i]) {
			sampcdf2[i] = 0.0;
			continue;
		}
		float s = 0.0f;
		for (int j = 0; j < 3; j++)
			s += eval2[j] * sqr(
				A2[0][j] * mesh2->normals[i][0] +
				A2[1][j] * mesh2->normals[i][1] +
				A2[2][j] * mesh2->normals[i][2]);
		sum_sampcdf2 += (sampcdf2[i] = s * weights2[i]);
	}

	if (!sum_sampcdf1 || !sum_sampcdf2) {
		if (verbose)
			dprintf("No overlap.\n");
		return -1.0f;
	}

	float cdf_scale1 = 1 / sum_sampcdf1;
	sampcdf1[0] *= cdf_scale1;
	for (int i = 1; i < nv1 - 1; i++)
		sampcdf1[i] = cdf_scale1 * sampcdf1[i] + sampcdf1[i-1];
	sampcdf1[nv1-1] = 1.0f;

	float cdf_scale2 = 1 / sum_sampcdf2;
	sampcdf2[0] *= cdf_scale2;
	for (int i = 1; i < nv2 - 1; i++)
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
	const int nv1 = weights1.size(), nv2 = weights2.size();
	sampcdf1.resize(nv1);
	sampcdf2.resize(nv2);

	double sum_sampcdf1 = 0, sum_sampcdf2 = 0;
	for (int i = 0; i < nv1; i++)
		sum_sampcdf1 += (sampcdf1[i] = weights1[i]);
	for (int i = 0; i < nv2; i++)
		sum_sampcdf2 += (sampcdf2[i] = weights2[i]);

	float cdf_scale1 = 1 / sum_sampcdf1;
	sampcdf1[0] *= cdf_scale1;
	for (int i = 1; i < nv1 - 1; i++)
		sampcdf1[i] = cdf_scale1 * sampcdf1[i] + sampcdf1[i-1];
	sampcdf1[nv1-1] = 1.0f;

	float cdf_scale2 = 1 / sum_sampcdf2;
	sampcdf2[0] *= cdf_scale2;
	for (int i = 1; i < nv2 - 1; i++)
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

	// Precompute normals and connectivity (used to determine boundaries)
	mesh1->need_normals();
	mesh2->need_normals();
	if (REJECT_BDY) {
		mesh1->need_faces();
		mesh1->need_neighbors();
		mesh1->need_adjacentfaces();
		mesh2->need_faces();
		mesh2->need_neighbors();
		mesh2->need_adjacentfaces();
	}

	// Initial distance and angle thresholds
	if (maxdist <= 0.0f) {
		mesh1->need_bsphere();
		mesh2->need_bsphere();
		maxdist = dist(xf1 * mesh1->bsphere.center,
		               xf2 * mesh2->bsphere.center);
		maxdist += mesh1->bsphere.r + mesh2->bsphere.r;
	}
	float angle_thresh = ANGLE_THRESH_MAX;

	// Weights and initial (uniform) CDFs
	const size_t nv1 = mesh1->vertices.size(), nv2 = mesh2->vertices.size();
	bool had_weights = (weights1.size() == nv1 && weights2.size() == nv2);
	if (!had_weights) {
		weights1.resize(nv1, 1.0f);
		weights2.resize(nv2, 1.0f);
	}

	vector<float> sampcdf1, sampcdf2;
	make_uniform_cdfs(weights1, sampcdf1, weights2, sampcdf2);
	float cdfincr = 2.0f / DESIRED_PAIRS;

	timestamp t_iters = now();
	if (verbose > 1) {
		dprintf("\nTime for preprocessing: %.3f msec.\n\n",
			(t_iters - t_start) * 1000.0f);
	}

	// Now the main ICP iterations
	ICP_xform_type iter_xform_type =
		(xform_type == ICP_TRANSLATION) ? ICP_TRANSLATION : ICP_RIGID;
	float err, min_err = 0.0f;
	int iter_of_min_err = -1;
	int iter;
	for (iter = 0; iter < MAX_ITERS; iter++) {
		// Should we recompute overlaps and CDFs?
		// The below uses the constant "2" so that we don't compute
		// overlaps until we have a few iterations under our belt,
		// in case things started out really far away.
		bool recompute = ((iter % CDF_UPDATE_INTERVAL) == 2);

		if (recompute) {
			if (!had_weights) {
				compute_overlaps(mesh1, mesh2, xf1, xf2,
						 kd1, kd2,
						 weights1, weights2,
						 maxdist, verbose);
			}

			// If we're recomputing CDFs, use uniform sampling
			// on this iteration to make sure that covariance
			// is unbiased.
			make_uniform_cdfs(weights1, sampcdf1,
			                  weights2, sampcdf2);
		}

		// Do an iteration
		err = ICP_iter(mesh1, mesh2, xf1, xf2, kd1, kd2,
		               weights1, weights2, sampcdf1, sampcdf2,
			       DESIRED_PAIRS, cdfincr, recompute,
			       maxdist, angle_thresh,
			       verbose, iter_xform_type);

		// Check resulting error
		if (err < 0) {
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
			iter, (now() - t_iters) * 1000.0f);
	}

	for (iter = 0; iter < FINAL_ITERS; iter++) {
		if (iter == 0) {
			cdfincr *= (float) DESIRED_PAIRS / DESIRED_PAIRS_FINAL;
			// Use uniform sampling so that the final error
			// we return is unbiased
			make_uniform_cdfs(weights1, sampcdf1, weights2, sampcdf2);
		}
		maxdist *= DIST_THRESH_MULT_FINAL / DIST_THRESH_MULT;
		err = ICP_iter(mesh1, mesh2, xf1, xf2, kd1, kd2,
		               weights1, weights2, sampcdf1, sampcdf2,
			       DESIRED_PAIRS_FINAL, cdfincr, false,
			       maxdist, angle_thresh,
			       verbose, iter_xform_type);
		if (err < 0) {
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
