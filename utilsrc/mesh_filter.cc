/*
Szymon Rusinkiewicz
Princeton University

mesh_filter.cc
Apply a variety of tranformations to a mesh
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
using namespace std;
using namespace trimesh;

#define ATOF(x) ((float) atof(x))


// Is this argument a floating-point number?
static bool isanumber(const char *c)
{
	if (!c || !*c)
		return false;
	char *endptr;
	strtod(c, &endptr);
	return (endptr && *endptr == '\0');
}


// Is this argument an integer?
static bool isanint(const char *c)
{
	if (!c || !*c)
		return false;
	char *endptr;
	strtol(c, &endptr, 10);
	return (endptr && *endptr == '\0');
}


// Transform the mesh by a matrix read from a file
void apply_xform(TriMesh *mesh, const char *xffilename)
{
	xform xf;
	if (!xf.read(xffilename))
		fprintf(stderr, "Couldn't open %s\n", xffilename);
	else
		apply_xform(mesh, xf);
}


// Transform the mesh by inverse of a matrix read from a file
void apply_ixform(TriMesh *mesh, const char *xffilename)
{
	xform xf;
	if (!xf.read(xffilename)) {
		fprintf(stderr, "Couldn't open %s\n", xffilename);
	} else {
		invert(xf);
		apply_xform(mesh, xf);
	}
}


// Clip mesh to the given bounding box file
bool clip(TriMesh *mesh, const char *bboxfilename)
{
	box b;
	if (!b.read(bboxfilename)) {
		fprintf(stderr, "Couldn't read bounding box %s\n", bboxfilename);
		return false;
	}

	clip(mesh, b);
	return true;
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s infile [options] [outfile]\n", myname);
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "	-color		Add per-vertex color\n");
	fprintf(stderr, "	-nocolor	Remove per-vertex color\n");
	fprintf(stderr, "	-conf		Add per-vertex confidence\n");
	fprintf(stderr, "	-noconf		Remove per-vertex confidence\n");
	fprintf(stderr, "	-tstrip		Convert to use triangle strips\n");
	fprintf(stderr, "	-notstrip	Unpack triangle strips to faces\n");
	fprintf(stderr, "	-nogrid		Unpack range grid to faces\n");
	fprintf(stderr, "	-nofaces	Delete all tris/tstrips/grid\n");
	fprintf(stderr, "	-reorder	Optimize order of vertices\n");
	fprintf(stderr, "	-orient		Auto-orient faces within the mesh\n");
	fprintf(stderr, "	-faceflip	Flip the order of vertices within each face\n");
	fprintf(stderr, "	-edgeflip	Optimize triangle connectivity by flipping edges\n");
	fprintf(stderr, "	-subdiv		Subdivide faces (planar)\n");
	fprintf(stderr, "	-loop		Perform Loop subdivision\n");
	fprintf(stderr, "	-fly		Perform butterfly subdivision\n");
	fprintf(stderr, "	-smooth s	Smooth surface with sigma=s*edgelength\n");
	fprintf(stderr, "	-bilat sd sr	Bilateral surface smoothing with domain, range sigmas\n");
	fprintf(stderr, "	-sharpen s	Sharpen surface with sigma=s*edgelength\n");
	fprintf(stderr, "	-smoothnorm s	Diffuse normals with sigma=s*edgelength\n");
	fprintf(stderr, "	-usmooth n	Perform n iterations of simple umbrella smoothing\n");
	fprintf(stderr, "	-tsmooth n	Perform n iterations of tangent-plane umbrella smoothing\n");
	fprintf(stderr, "	-lmsmooth n	Perform n iterations of Taubin's lambda-mu smoothing\n");
	fprintf(stderr, "	-nsmooth n	Perform n iterations of umbrella smoothing on the normals\n");
	fprintf(stderr, "	-inflate s	Create offset surface s*edgelength away\n");
	fprintf(stderr, "	-noisify s	Add O(s*edgelength) noise to each vertex\n");
	fprintf(stderr, "	-share tol	Merge (\"share\") vertices within tol*edgelength\n");
	fprintf(stderr, "	-clip bbox	Clip to the given bbox (file has 6 numbers)\n");
	fprintf(stderr, "	-xform file.xf	Transform by the given matrix\n");
	fprintf(stderr, "	-ixform file.xf	Transform by inverse of matrix\n");
	fprintf(stderr, "	-rot r x y z	Rotate r degrees around axis (x,y,z)\n");
	fprintf(stderr, "	-trans x y z	Translate by (x,y,z)\n");
	fprintf(stderr, "	-scale s	Uniform scale by s\n");
	fprintf(stderr, "	-scale x y z	Scale by (x,y,z)\n");
	fprintf(stderr, "	-scale s x y z	Scale by s in direction (x,y,z)\n");
	fprintf(stderr, "	-center		Translate so center of mass is at (0,0,0)\n");
	fprintf(stderr, "	-bbcenter	Translate so center of bbox is at (0,0,0)\n");
	fprintf(stderr, "	-varnorm	Scale so variance (RMS distance) from center is 1\n");
	fprintf(stderr, "	-bbnorm		Scale so bbox has maximum extent 1\n");
	fprintf(stderr, "	-pcarot		Rotate so that principal axes lie along X, Y, Z\n");
	fprintf(stderr, "	-pcasnap	As above, but only rotate by 90/180 degrees\n");
	fprintf(stderr, "	-rmunused	Remove unreferenced vertices\n");
	fprintf(stderr, "	-rmslivers	Remove long, skinny faces\n");
	fprintf(stderr, "	-erode		Enlarge boundaries by removing boundary vertices\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 3)
		usage(argv[0]);
	const char *filename = argv[1];

	TriMesh *themesh = TriMesh::read(filename);
	if (!themesh)
		usage(argv[0]);

	bool have_tstrips = !themesh->tstrips.empty();
	for (int i = 2; i < argc; i++) {
		if (!strcmp(argv[i], "-color") ||
		    !strcmp(argv[i], "-colors")) {
			if (themesh->colors.empty()) {
				int nv = themesh->vertices.size();
				themesh->colors.resize(nv, Color::white());
			}
		} else if (!strcmp(argv[i], "-nocolor") ||
		           !strcmp(argv[i], "-nocolors")) {
			themesh->colors.clear();
		} else if (!strcmp(argv[i], "-conf")) {
			if (isanumber(argv[i+1])) {
				float desired_conf = ATOF(argv[i+1]);
				int nv = themesh->vertices.size();
				themesh->confidences.clear();
				themesh->confidences.resize(nv, desired_conf);
				i++;
			} else if (themesh->confidences.empty()) {
				int nv = themesh->vertices.size();
				themesh->confidences.resize(nv, 1);
			}
		} else if (!strcmp(argv[i], "-noconf")) {
			themesh->confidences.clear();
		} else if (!strcmp(argv[i], "-tstrip") ||
		           !strcmp(argv[i], "-tstrips") ||
		           !strcmp(argv[i], "-strip") ||
		           !strcmp(argv[i], "-strips")) {
			themesh->need_tstrips();
			have_tstrips = true;
		} else if (!strcmp(argv[i], "-notstrip") ||
		           !strcmp(argv[i], "-notstrips") ||
		           !strcmp(argv[i], "-nostrip") ||
		           !strcmp(argv[i], "-nostrips") ||
		           !strcmp(argv[i], "-unstrip")) {
			themesh->need_faces();
			themesh->tstrips.clear();
			have_tstrips = false;
		} else if (!strcmp(argv[i], "-nogrid")) {
			themesh->need_faces();
			themesh->grid.clear();
		} else if (!strcmp(argv[i], "-nofaces")) {
			themesh->faces.clear();
			themesh->tstrips.clear();
			themesh->grid.clear();
		} else if (!strcmp(argv[i], "-reorder")) {
			reorder_verts(themesh);
		} else if (!strcmp(argv[i], "-orient")) {
			orient(themesh);
		} else if (!strcmp(argv[i], "-faceflip")) {
			faceflip(themesh);
		} else if (!strcmp(argv[i], "-edgeflip")) {
			edgeflip(themesh);
		} else if (!strcmp(argv[i], "-subdiv")) {
			subdiv(themesh, SUBDIV_PLANAR);
		} else if (!strcmp(argv[i], "-loop")) {
			subdiv(themesh, SUBDIV_LOOP);
		} else if (!strcmp(argv[i], "-loop2")) {
			subdiv(themesh, SUBDIV_LOOP_ORIG);
		} else if (!strcmp(argv[i], "-loop3")) {
			subdiv(themesh, SUBDIV_LOOP_NEW);
		} else if (!strcmp(argv[i], "-fly")) {
			subdiv(themesh, SUBDIV_BUTTERFLY_MODIFIED);
		} else if (!strcmp(argv[i], "-fly2")) {
			subdiv(themesh, SUBDIV_BUTTERFLY);
		} else if (!strcmp(argv[i], "-smooth")) {
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-smooth requires one float parameter: s\n\n");
				usage(argv[0]);
			}
			float amount = ATOF(argv[i]) * themesh->feature_size();
			smooth_mesh(themesh, amount);
			themesh->pointareas.clear();
			themesh->normals.clear();
		} else if (!strcmp(argv[i], "-bilat")) {
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-smooth requires two float parameters: sd, sr\n\n");
				usage(argv[0]);
			}
			float fs = themesh->feature_size();
			float sd = ATOF(argv[i]) * fs;
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-smooth requires two float parameters: sd, sr\n\n");
				usage(argv[0]);
			}
			float sr = ATOF(argv[i]) * fs;
			bilateral_smooth_mesh(themesh, sd, sr);
			themesh->normals.clear();
		} else if (!strcmp(argv[i], "-sharpen")) {
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-sharpen requires one float parameter: s\n\n");
				usage(argv[0]);
			}
			float amount = ATOF(argv[i]) * themesh->feature_size();
			vector<point> origverts = themesh->vertices;
			smooth_mesh(themesh, amount);
			for (size_t v = 0; v < themesh->vertices.size(); v++)
				themesh->vertices[v] += 2.0f *
					(origverts[v] - themesh->vertices[v]);
			themesh->pointareas.clear();
			themesh->normals.clear();
		} else if (!strcmp(argv[i], "-smoothnorm")) {
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-smoothnorm requires one float parameter: s\n\n");
				usage(argv[0]);
			}
			float amount = ATOF(argv[i]) * themesh->feature_size();
			diffuse_normals(themesh, amount);
		} else if (!strcmp(argv[i], "-usmooth")) {
			i++;
			if (!(i < argc && isanint(argv[i]))) {
				fprintf(stderr, "\n-usmooth requires one int parameter: n\n\n");
				usage(argv[0]);
			}
			int niters = atoi(argv[i]);
			for (int iter = 0; iter < niters; iter++)
				umbrella(themesh, 0.5f);
		} else if (!strcmp(argv[i], "-tsmooth")) {
			i++;
			if (!(i < argc && isanint(argv[i]))) {
				fprintf(stderr, "\n-tsmooth requires one int parameter: n\n\n");
				usage(argv[0]);
			}
			int niters = atoi(argv[i]);
			for (int iter = 0; iter < niters; iter++)
				umbrella(themesh, 0.5f, true);
		} else if (!strcmp(argv[i], "-lmsmooth")) {
			i++;
			if (!(i < argc && isanint(argv[i]))) {
				fprintf(stderr, "\n-lmsmooth requires one int parameter: n\n\n");
				usage(argv[0]);
			}
			int niters = atoi(argv[i]);
			lmsmooth(themesh, niters);
		} else if (!strcmp(argv[i], "-nsmooth")) {
			i++;
			if (!(i < argc && isanint(argv[i]))) {
				fprintf(stderr, "\n-nsmooth requires one int parameter: n\n\n");
				usage(argv[0]);
			}
			int niters = atoi(argv[i]);
			for (int iter = 0; iter < niters; iter++)
				numbrella(themesh, 0.5f);
		} else if (!strcmp(argv[i], "-inflate")) {
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-inflate requires one float parameter: s\n\n");
				usage(argv[0]);
			}
			float amount = ATOF(argv[i]) * themesh->feature_size();
			inflate(themesh, amount);
		} else if (!strcmp(argv[i], "-noisify")) {
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-noisify requires one float parameter: s\n\n");
				usage(argv[0]);
			}
			float amount = ATOF(argv[i]) * themesh->feature_size();
			noisify(themesh, amount);
		} else if (!strcmp(argv[i], "-share")) {
			i++;
			if (!(i < argc && isanumber(argv[i]))) {
				fprintf(stderr, "\n-share requires one float parameter: tol\n\n");
				usage(argv[0]);
			}
			float tol = ATOF(argv[i]) * themesh->feature_size();
			shared(themesh, tol);
		} else if (!strcmp(argv[i], "-clip")) {
			i++;
			if (!(i < argc)) {
				fprintf(stderr, "\n-clip requires one argument\n\n");
				usage(argv[0]);
			}
			if (!clip(themesh, argv[i]))
				usage(argv[0]);
		} else if (!strcmp(argv[i], "-xf") ||
		           !strcmp(argv[i], "-xform")) {
			i++;
			if (!(i < argc)) {
				fprintf(stderr, "\n-xform requires one argument\n\n");
				usage(argv[0]);
			}
			apply_xform(themesh, argv[i]);
		} else if (!strcmp(argv[i], "-ixf") ||
		           !strcmp(argv[i], "-ixform")) {
			i++;
			if (!(i < argc)) {
				fprintf(stderr, "\n-ixform requires one argument\n\n");
				usage(argv[0]);
			}
			apply_ixform(themesh, argv[i]);
		} else if (!strcmp(argv[i], "-rot") ||
		           !strcmp(argv[i], "-rotate")) {
			i += 4;
			if (!(i < argc &&
			      isanumber(argv[i]) && isanumber(argv[i-1]) &&
			      isanumber(argv[i-2]) && isanumber(argv[i-3]))) {
				fprintf(stderr, "\n-rot requires four arguments\n\n");
				usage(argv[0]);
			}
			vec ax(ATOF(argv[i-2]), ATOF(argv[i-1]), ATOF(argv[i]));
			float ang = radians(ATOF(argv[i-3]));
			rot(themesh, ang, ax);
		} else if (!strcmp(argv[i], "-trans") ||
		           !strcmp(argv[i], "-translate")) {
			i += 3;
			if (!(i < argc && isanumber(argv[i]) &&
			      isanumber(argv[i-1]) && isanumber(argv[i-2]))) {
				fprintf(stderr, "\n-trans requires three arguments\n\n");
				usage(argv[0]);
			}
			vec t(ATOF(argv[i-2]), ATOF(argv[i-1]), ATOF(argv[i]));
			trans(themesh, t);
		} else if (!strcmp(argv[i], "-scale")) {
			int nargs = 0;
			float args[4];
			while (nargs < 4) {
				if (++i >= argc)
					break;
				if (!isanumber(argv[i]) ||
				    !sscanf(argv[i], "%f", &(args[nargs]))) {
					--i;
					break;
				}
				nargs++;
			}
			if (!(i < argc) || nargs == 0 || nargs == 2) {
				fprintf(stderr, "\n-scale requires 1, 3, or 4 arguments\n\n");
				usage(argv[0]);
			}
			xform s = xform::scale(args[0]);
			if (nargs == 3)
				s = xform::scale(args[0], args[1], args[2]);
			else if (nargs == 4)
				s = xform::scale(args[0], args[1], args[2], args[3]);
			apply_xform(themesh, s);
		} else if (!strcmp(argv[i], "-center")) {
			trans(themesh, -mesh_center_of_mass(themesh));
		} else if (!strcmp(argv[i], "-bbcenter")) {
			themesh->need_bbox();
			trans(themesh, -themesh->bbox.center());
		} else if (!strcmp(argv[i], "-varnorm")) {
			normalize_variance(themesh);
		} else if (!strcmp(argv[i], "-bbnorm")) {
			themesh->need_bbox();
			vec l = themesh->bbox.size();
			float ll = max(max(l[0], l[1]), l[2]);
			trans(themesh, -themesh->bbox.center());
			float s = 1.0f / ll;
			scale(themesh, s);
			trans(themesh, themesh->bbox.center());
		} else if (!strcmp(argv[i], "-pcarot")) {
			pca_rotate(themesh);
		} else if (!strcmp(argv[i], "-pcasnap")) {
			pca_snap(themesh);
		} else if (!strcmp(argv[i], "-rmunused")) {
			remove_unused_vertices(themesh);
		} else if (!strcmp(argv[i], "-rmslivers")) {
			remove_sliver_faces(themesh);
		} else if (!strcmp(argv[i], "-erode")) {
			erode(themesh);
		} else if (i == argc - 1 &&
		           (argv[i][0] != '-' || argv[i][1] == '\0')) {
			if (have_tstrips && themesh->tstrips.empty())
				themesh->need_tstrips();
			themesh->write(argv[i]);
		} else {
			fprintf(stderr, "\nUnrecognized option [%s]\n\n", argv[i]);
			usage(argv[0]);
		}
	}
}
