#ifndef TRIMESH_H
#define TRIMESH_H
/*
Szymon Rusinkiewicz
Princeton University

TriMesh.h
Class for triangle meshes.
*/

#include "Vec.h"
#include "Box.h"
#include "Color.h"
#include "strutil.h"
#include <vector>


namespace trimesh {

template <class T>
static inline void clear_and_release(::std::vector<T> &v)
{
	// Standard trick to release a vector's storage, since clear() doesn't
	::std::vector<T>().swap(v);
}

class TriMesh {
public:
	//
	// Types
	//
	typedef Vec<3,int> Face;
	typedef Box<3,float> BBox;

	struct BSphere {
		point center;
		float r;
		bool valid;
		BSphere() : valid(false)
			{}
	};

	//
	// Enums
	//
	enum TstripRep { TSTRIP_LENGTH, TSTRIP_TERM };
	enum { GRID_INVALID = -1 };
	enum StatOp {
		STAT_MIN, STAT_MINABS, STAT_MAX, STAT_MAXABS,
		STAT_SUM, STAT_SUMABS, STAT_SUMSQR,
		STAT_MEAN, STAT_MEANABS, STAT_RMS,
		STAT_MEDIAN, STAT_STDEV };
	enum StatVal { STAT_VALENCE, STAT_FACEAREA, STAT_ANGLE,
		STAT_DIHEDRAL, STAT_EDGELEN, STAT_X, STAT_Y, STAT_Z };

	//
	// Constructor
	//
	TriMesh() : grid_width(-1), grid_height(-1), flag_curr(0)
		{}

	//
	// Members
	//

	// The basics: vertices and faces
	::std::vector<point> vertices;
	::std::vector<Face> faces;

	// Triangle strips
	::std::vector<int> tstrips;

	// Grid, if present
	::std::vector<int> grid;
	int grid_width, grid_height;

	// Other per-vertex properties
	::std::vector<Color> colors;
	::std::vector<float> confidences;
	::std::vector<unsigned> flags;
	unsigned flag_curr;

	// Computed per-vertex properties
	::std::vector<vec> normals;
	::std::vector<vec> pdir1, pdir2;
	::std::vector<float> curv1, curv2;
	::std::vector< Vec<4,float> > dcurv;
	::std::vector<vec> cornerareas;
	::std::vector<float> pointareas;

	// Bounding structures
	BBox bbox;
	BSphere bsphere;

	// Connectivity structures:
	//  For each vertex, all neighboring vertices
	::std::vector< ::std::vector<int> > neighbors;
	//  For each vertex, all neighboring faces
	::std::vector< ::std::vector<int> > adjacentfaces;
	//  For each face, the three faces attached to its edges
	//  (for example, across_edge[3][2] is the number of the face
	//   that's touching the edge opposite vertex 2 of face 3)
	::std::vector<Face> across_edge;

	//
	// Compute all this stuff...
	//
	void need_faces()
	{
		if (!faces.empty())
			return;
		if (!tstrips.empty())
			unpack_tstrips();
		else if (!grid.empty())
			triangulate_grid();
	}
	void need_tstrips(TstripRep rep = TSTRIP_LENGTH);
	void convert_strips(TstripRep rep);
	void unpack_tstrips();
	void resize_grid(int width, int height)
	{
		grid_width = width;
		grid_height = height;
		grid.clear();
		grid.resize(grid_width * grid_height, GRID_INVALID);
	}
	void triangulate_grid(bool remove_slivers = true);
	void need_normals(bool simple_area_weighted = false);
	void need_curvatures();
	void need_dcurv();
	void need_pointareas();
	void need_bbox();
	void need_bsphere();
	void need_neighbors();
	void need_adjacentfaces();
	void need_across_edge();

	//
	// Delete everything and release storage
	//
	void clear_vertices()      { clear_and_release(vertices); }
	void clear_faces()         { clear_and_release(faces); }
	void clear_tstrips()       { clear_and_release(tstrips); }
	void clear_grid()          { clear_and_release(grid);
	                             grid_width = grid_height = -1;}
	void clear_colors()        { clear_and_release(colors); }
	void clear_confidences()   { clear_and_release(confidences); }
	void clear_flags()         { clear_and_release(flags); flag_curr = 0; }
	void clear_normals()       { clear_and_release(normals); }
	void clear_curvatures()    { clear_and_release(pdir1);
	                             clear_and_release(pdir2);
	                             clear_and_release(curv1);
	                             clear_and_release(curv2); }
	void clear_dcurv()         { clear_and_release(dcurv); }
	void clear_pointareas()    { clear_and_release(pointareas);
	                             clear_and_release(cornerareas); }
	void clear_bbox()          { bbox.clear(); }
	void clear_bsphere()       { bsphere.valid = false; }
	void clear_neighbors()     { clear_and_release(neighbors); }
	void clear_adjacentfaces() { clear_and_release(adjacentfaces); }
	void clear_across_edge()   { clear_and_release(across_edge); }
	void clear()
	{
		clear_vertices(); clear_faces(); clear_tstrips(); clear_grid();
		clear_colors(); clear_confidences(); clear_flags();
		clear_normals(); clear_curvatures(); clear_dcurv();
		clear_pointareas(); clear_bbox(); clear_bsphere();
		clear_neighbors(); clear_adjacentfaces(); clear_across_edge();
	}

	//
	// Input and output
	//
protected:
	static bool read_helper(const char *filename, TriMesh *mesh);
public:
	static TriMesh *read(const char *filename);
	static TriMesh *read(const ::std::string &filename);
	bool write(const char *filename);
	bool write(const ::std::string &filename);


	//
	// Useful queries
	//

	// Is vertex v on the mesh boundary?
	inline bool is_bdy(int v)
	{
		if (unlikely(neighbors.empty())) need_neighbors();
		if (unlikely(adjacentfaces.empty())) need_adjacentfaces();
		return neighbors[v].size() != adjacentfaces[v].size();
	}

	// Centroid of face f
	inline vec centroid(int f)
	{
		if (unlikely(faces.empty())) need_faces();
		return (1.0f / 3.0f) *
			(vertices[faces[f][0]] +
			 vertices[faces[f][1]] +
			 vertices[faces[f][2]]);
	}

	// Normal of face f
	inline vec trinorm(int f)
	{
		if (unlikely(faces.empty())) need_faces();
		return trimesh::trinorm(vertices[faces[f][0]],
			vertices[faces[f][1]], vertices[faces[f][2]]);
	}

	// Angle of corner j in triangle i
	inline float cornerangle(int i, int j)
	{
		using namespace ::std;

		if (unlikely(faces.empty())) need_faces();
		const point &p0 = vertices[faces[i][j]];
		const point &p1 = vertices[faces[i][NEXT_MOD3(j)]];
		const point &p2 = vertices[faces[i][PREV_MOD3(j)]];
		return acos((p1 - p0) DOT (p2 - p0));
	}

	// Dihedral angle between face i and face across_edge[i][j]
	inline float dihedral(int i, int j)
	{
		if (unlikely(across_edge.empty())) need_across_edge();
		if (unlikely(across_edge[i][j] < 0)) return 0.0f;
		vec mynorm = trinorm(i);
		vec othernorm = trinorm(across_edge[i][j]);
		float ang = angle(mynorm, othernorm);
		vec towards = 0.5f * (vertices[faces[i][NEXT_MOD3(j)]] +
		                      vertices[faces[i][PREV_MOD3(j)]]) -
		              vertices[faces[i][j]];
		if ((towards DOT othernorm) < 0.0f)
			return M_PIf + ang;
		else
			return M_PIf - ang;
	}

	// Statistics
	float stat(StatOp op, StatVal val);
	float feature_size();

	//
	// Debugging
	//

	// Debugging printout, controllable by a "verbose"ness parameter
	static int verbose;
	static void set_verbose(int);
	static void (*dprintf_hook)(const char *);
	static void set_dprintf_hook(void (*hook)(const char *));
	static void dprintf(const char *format, ...);

	// Same as above, but fatal-error printout
	static void (*eprintf_hook)(const char *);
	static void set_eprintf_hook(void (*hook)(const char *));
	static void eprintf(const char *format, ...);

};

} // namespace trimesh

#endif
