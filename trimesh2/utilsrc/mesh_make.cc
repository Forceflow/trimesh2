/*
Szymon Rusinkiewicz
Princeton University

mesh_make.cc
Create various kinds of meshes for testing...
*/

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "noise3d.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
using namespace std;
using namespace trimesh;

#define ATOF(x) ((float) atof(x))


static inline void mkpoint(TriMesh *mesh, float x, float y, float z)
{
	mesh->vertices.push_back(point(x,y,z));
}

static inline void mkface(TriMesh *mesh, int v1, int v2, int v3)
{
	mesh->faces.push_back(TriMesh::Face(v1, v2, v3));
}

static inline void mkquad(TriMesh *mesh, int ll, int lr, int ul, int ur)
{
	mkface(mesh, ll, lr, ur);
	mkface(mesh, ll, ur, ul);
}

static inline void tess4(TriMesh *mesh, int v1, int v2, int v3, int v4)
{
	point c = 0.25f * (mesh->vertices[v1] + mesh->vertices[v2] +
			   mesh->vertices[v3] + mesh->vertices[v4]);
	mkpoint(mesh, c[0], c[1], c[2]);
	int ci = mesh->vertices.size() - 1;
	mkface(mesh, ci, v1, v2);
	mkface(mesh, ci, v2, v3);
	mkface(mesh, ci, v3, v4);
	mkface(mesh, ci, v4, v1);
}

static inline void tess5(TriMesh *mesh, int v1, int v2, int v3, int v4, int v5)
{
	point c = 0.2f * (mesh->vertices[v1] + mesh->vertices[v2] +
			  mesh->vertices[v3] + mesh->vertices[v4] +
			  mesh->vertices[v5]);
	mkpoint(mesh, c[0], c[1], c[2]);
	int ci = mesh->vertices.size() - 1;
	mkface(mesh, ci, v1, v2);
	mkface(mesh, ci, v2, v3);
	mkface(mesh, ci, v3, v4);
	mkface(mesh, ci, v4, v5);
	mkface(mesh, ci, v5, v1);
}


TriMesh *make_plane(int tess_x, int tess_y)
{
	if (tess_x < 1)
		tess_x = 1;
	if (tess_y < 1)
		tess_y = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve((tess_x+1)*(tess_y+1));
	for (int j = 0; j < tess_y+1; j++) {
		float y = -1.0f + 2.0f * j / tess_y;
		for (int i = 0; i < tess_x+1; i++) {
			float x = -1.0f + 2.0f * i / tess_x;
			mkpoint(mesh, x, y, 0);
		}
	}

	mesh->faces.reserve(2*tess_x*tess_y);
	for (int j = 0; j < tess_y; j++) {
		for (int i = 0; i < tess_x; i++) {
			int ind = i + j * (tess_x+1);
			mkquad(mesh, ind, ind+1, ind+tess_x+1, ind+tess_x+2);
		}
	}

	return mesh;
}


TriMesh *make_bump(int tess, float sigma)
{
	if (tess < 2)
		tess = 2;

	TriMesh *mesh = make_plane(tess, tess);
	int nv = mesh->vertices.size();
	float escale = -0.5f * sqr(1.0f / sigma);
	for (int i = 0; i < nv; i++) {
		point &p = mesh->vertices[i];
		p[2] = exp(escale * (sqr(p[0]) + sqr(p[1])));
	}
	return mesh;
}


TriMesh *make_wave(int tess, float omega)
{
	if (tess < 2)
		tess = 2;

	TriMesh *mesh = make_plane(tess, tess);
	int nv = mesh->vertices.size();
	float scale = 1.0f / omega;
	for (int i = 0; i < nv; i++) {
		point &p = mesh->vertices[i];
		p[2] = scale * sin(omega * p[0]) * sin(omega * p[1]);
	}
	return mesh;
}


TriMesh *make_frac(int tess)
{
	if (tess < 2)
		tess = 2;

	TriMesh *mesh = make_plane(tess, tess);
	int nv = mesh->vertices.size();
	PerlinNoise3D noise(tess, tess, tess);
	for (int i = 0; i < nv; i++) {
		point &p = mesh->vertices[i];
		p[2] = noise.lookup(p[0], p[1], 0.5f);
	}
	return mesh;
}

TriMesh *make_cube(int tess)
{
	if (tess < 1)
		tess = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(6*sqr(tess)+2);
	for (int j = 0; j < tess+1; j++) {
		float y = 1.0f - 2.0f * j / tess;
		for (int i = 0; i < tess+1; i++) {
			float x = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, x, y, -1);
		}
	}
	for (int j = 1; j < tess; j++) {
		float z = -1.0f + 2.0f * j / tess;
		for (int i = 0; i < tess; i++) {
			float x = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, x, -1, z);
		}
		for (int i = 0; i < tess; i++) {
			float y = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, 1, y, z);
		}
		for (int i = 0; i < tess; i++) {
			float x = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, x, 1, z);
		}
		for (int i = 0; i < tess; i++) {
			float y = 1.0f - 2.0f * i / tess;
			mkpoint(mesh, -1, y, z);
		}
	}
	for (int j = 0; j < tess+1; j++) {
		float y = -1.0f + 2.0f * j / tess;
		for (int i = 0; i < tess+1; i++) {
			float x = -1.0f + 2.0f * i / tess;
			mkpoint(mesh, x, y, 1);
		}
	}

	mesh->faces.reserve(12*sqr(tess));
	for (int j = 0; j < tess; j++) {
		for (int i = 0; i < tess; i++) {
			int ind = i + j * (tess+1);
			mkquad(mesh, ind, ind+tess+1, ind+1, ind+tess+2);
		}
	}

	int topstart = sqr(tess+1) + 4*tess*(tess-1);
	for (int j = 0; j < tess; j++) {
		int next = sqr(tess+1) + 4*tess*(j-1);
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = sqr(tess+1)-1 - i;
				lr = ll - 1;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = tess*(tess+1) - i*(tess+1);
				lr = ll - (tess+1);
			}
			if (j == tess-1) {
				ul = topstart + tess + i*(tess+1);
				ur = ul + (tess+1);
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = i;
				lr = i + 1;
			}
			if (j == tess-1) {
				ul = topstart + sqr(tess+1)-1 - i;
				ur = ul - 1;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
		for (int i = 0; i < tess; i++) {
			int ll = next++;
			int lr = ll + 1;
			int ul = ll + 4*tess;
			int ur = ul + 1;
			if (j == 0) {
				ll = tess + i*(tess+1);
				lr = ll + (tess+1);
			}
			if (j == tess-1) {
				ul = topstart + tess*(tess+1) - i*(tess+1);
				ur = ul - (tess+1);
			}
			if (i == tess-1) {
				if (j != 0)
					lr -= 4*tess;
				if (j != tess-1)
					ur -= 4*tess;
			}
			mkquad(mesh, ll, lr, ul, ur);
		}
	}
	for (int j = 0; j < tess; j++) {
		for (int i = 0; i < tess; i++) {
			int ind = topstart + i + j * (tess+1);
			mkquad(mesh, ind, ind+1, ind+tess+1, ind+tess+2);
		}
	}

	return mesh;
}


TriMesh *make_disc(int tess_th, int tess_r)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_r < 1)
		tess_r = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(1+tess_th*tess_r);
	mkpoint(mesh, 0, 0, 0);
	for (int j = 1; j <= tess_r; j++) {
		float r = (float) j / tess_r;
		for (int i = 0; i < tess_th; i++) {
			float th = M_TWOPIf * i / tess_th;
			mkpoint(mesh, r*cos(th), r*sin(th), 0);
		}
	}

	mesh->faces.reserve(2*tess_th*tess_r-tess_th);
	for (int i = 0; i < tess_th; i++)
		mkface(mesh, 0, i+1, ((i+1)%tess_th)+1);
	for (int j = 1; j < tess_r; j++) {
		int base = 1 + (j-1) * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i, base+tess_th+i1,
				base+i, base+i1);
		}
	}

	return mesh;
}


TriMesh *make_cyl(int tess_th, int tess_h, float r = 1.0f)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_h < 1)
		tess_h = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(tess_th * (tess_h+1));
	for (int j = 0; j <= tess_h; j++) {
		float z = -1.0f + 2.0f * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_TWOPIf * i / tess_th;
			mkpoint(mesh, r*cos(th), r*sin(th), z);
		}
	}

	mesh->faces.reserve(2*tess_th*tess_h);
	for (int j = 0; j < tess_h; j++) {
		int base = j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base + i, base + i1,
				base+tess_th+i, base+tess_th+i1);
		}
	}

	return mesh;
}


TriMesh *make_ccyl(int tess_th, int tess_h, float r = 1.0f)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_h < 1)
		tess_h = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(2+3*tess_th*tess_h-tess_th);

	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j <= tess_h; j++) {
		float rr = r * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_TWOPIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), -1);
		}
	}
	int side_start = mesh->vertices.size();
	for (int j = 1; j < tess_h; j++) {
		float z = -1.0f + 2.0f * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_TWOPIf * i / tess_th;
			mkpoint(mesh, r*cos(th), r*sin(th), z);
		}
	}
	int top_start = mesh->vertices.size();
	for (int j = tess_h; j > 0; j--) {
		float rr = r * j / tess_h;
		for (int i = 0; i < tess_th; i++) {
			float th = M_TWOPIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), 1);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(6*tess_th*tess_h - 2*tess_th);

	for (int i = 0; i < tess_th; i++)
		mkface(mesh, 0, ((i+1)%tess_th)+1, i+1);
	for (int j = 1; j < tess_h; j++) {
		int base = 1 + (j-1) * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i1, base+tess_th+i,
				base+i1, base+i);
		}
	}

	for (int j = 0; j < tess_h; j++) {
		int base = side_start - tess_th + j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base + i, base + i1,
				base+tess_th+i, base+tess_th+i1);
		}
	}

	for (int j = 0; j < tess_h-1; j++) {
		int base = top_start + j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i1, base+tess_th+i,
				base+i1, base+i);
		}
	}
	int base = top_start + (tess_h-1)*tess_th;
	for (int i = 0; i < tess_th; i++)
		mkface(mesh, base+i, base+((i+1)%tess_th), base+tess_th);

	return mesh;
}


TriMesh *make_cone(int tess_th, int tess_r, float r = 1.0f)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_r < 1)
		tess_r = 1;

	TriMesh *mesh = make_disc(tess_th, tess_r);
	for (int j = 0; j <= tess_r; j++) {
		float z = 1.0f - 2.0f * j / tess_r;
		if (j == 0) {
			mesh->vertices[0][2] = z;
			continue;
		}
		for (int i = 0; i < tess_th; i++) {
			int ind = 1 + i + (j-1)*tess_th;
			mesh->vertices[ind][0] *= r;
			mesh->vertices[ind][1] *= r;
			mesh->vertices[ind][2] = z;
		}
	}
	return mesh;
}


TriMesh *make_ccone(int tess_th, int tess_r, float r = 1.0f)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_r < 1)
		tess_r = 1;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(2+2*tess_th*tess_r-tess_th);

	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j <= tess_r; j++) {
		float rr = r * j / tess_r;
		for (int i = 0; i < tess_th; i++) {
			float th = M_TWOPIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), -1);
		}
	}
	int side_start = mesh->vertices.size();
	for (int j = 1; j < tess_r; j++) {
		float z = -1.0f + 2.0f * j / tess_r;
		float rr = r * (tess_r - j) / tess_r;
		for (int i = 0; i < tess_th; i++) {
			float th = M_TWOPIf * i / tess_th;
			mkpoint(mesh, rr*cos(th), rr*sin(th), z);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(4*tess_th*tess_r - 2*tess_th);

	for (int i = 0; i < tess_th; i++)
		mkface(mesh, 0, ((i+1)%tess_th)+1, i+1);
	for (int j = 1; j < tess_r; j++) {
		int base = 1 + (j-1) * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base+tess_th+i1, base+tess_th+i,
				base+i1, base+i);
		}
	}

	for (int j = 0; j < tess_r-1; j++) {
		int base = side_start - tess_th + j * tess_th;
		for (int i = 0; i < tess_th; i++) {
			int i1 = (i+1)%tess_th;
			mkquad(mesh, base + i, base + i1,
				base+tess_th+i, base+tess_th+i1);
		}
	}

	int base = side_start + (tess_r-2)*tess_th;
	for (int i = 0; i < tess_th; i++)
		mkface(mesh, base+i, base+((i+1)%tess_th), base+tess_th);

	return mesh;
}


TriMesh *make_torus(int tess_th, int tess_ph, float r = 0.25f)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_cyl(tess_ph, tess_th);

	for (int i = 0; i < tess_ph; i++)
		mesh->vertices.pop_back();
	for (size_t i = 0; i < mesh->faces.size(); i++) {
		mesh->faces[i][0] %= mesh->vertices.size();
		mesh->faces[i][1] %= mesh->vertices.size();
		mesh->faces[i][2] %= mesh->vertices.size();
	}

	for (int j = 0; j < tess_th; j++) {
		float th = M_TWOPIf * j / tess_th;
		vec circlepos(cos(th), sin(th), 0);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_TWOPIf * i / tess_ph;
			mesh->vertices[i+j*tess_ph] = circlepos +
						      cos(ph)*r*circlepos +
						      sin(ph)*r*vec(0,0,-1);
		}
	}
	return mesh;
}


TriMesh *make_knot(int tess_th, int tess_ph, float r = 0.2f)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_torus(tess_th, tess_ph, r);

	for (int j = 0; j < tess_th; j++) {
		float th = M_TWOPIf * j / tess_th;
		vec pos(2.0f*sin(2.0f*th) + cos(th),
			2.0f*cos(2.0f*th) + sin(th),
			cos(3.0f*th));
		pos /= 3.0f;
		vec vel( 4.0f*cos(2.0f*th) - sin(th),
			-4.0f*sin(2.0f*th) + cos(th),
			-3.0f*sin(3.0f*th));
		normalize(vel);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_TWOPIf * i / tess_ph;
			vec u = vel CROSS vec(0,0,1);
			normalize(u);
			vec v = u CROSS vel;
			mesh->vertices[i+j*tess_ph] =
				pos + cos(ph)*r*v + sin(ph)*r*u;
		}
	}
	return mesh;
}

TriMesh *make_klein(int tess_th, int tess_ph)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_torus(tess_th, tess_ph, 0.2f);

	for (int j = 0; j < tess_th; j++) {
		float v = 8.0f * M_PIf * j / tess_th;
		int part = 0;
		for (int i = 0; i < 7; i++)
			if (v > M_PIf)
				v -= M_PIf, part++;
		float flipx = 1.0f;
		if (part > 3)
			part -= 4, flipx = -1.0f;
		if (part == 0 || part == 2)
			v = M_PIf - v;
		for (int i = 0; i < tess_ph; i++) {
			float u = M_TWOPIf * i / tess_ph;
			// Based on formulae from
			// http://emsh.calarts.edu/~mathart/sw/klein/Klein.html
			if (part == 0)
				mesh->vertices[i+j*tess_ph] = point(
					(2.5f + 1.5f * cos(v)) * flipx * cos(u),
					(2.5f + 1.5f * cos(v)) * sin(u),
					-2.5f * sin(v));
			else if (part == 1)
				mesh->vertices[i+j*tess_ph] = point(
					(2.5f + 1.5f * cos(v)) * flipx * cos(u),
					(2.5f + 1.5f * cos(v)) * sin(u),
					3.0f * v);
			else if (part == 2)
				mesh->vertices[i+j*tess_ph] = point(
					2.0f + (2.0f - flipx * cos(u)) * cos(v),
					sin(u),
					3.0f * M_PIf + (2.0f - flipx * cos(u)) * sin(v));
			else
				mesh->vertices[i+j*tess_ph] = point(
					2.0f + 2.0f * cos(v) - flipx * cos(u),
					sin(u),
					3.0f * (M_PIf - v));
		}
	}
	return mesh;
}


TriMesh *make_helix(int tess_th, int tess_ph, float turns, float r = 0.25f)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = make_cyl(tess_ph, tess_th);

	const float angle = atan(1.0f / M_PIf);
	for (int j = 0; j <= tess_th; j++) {
		float z = (-1.0f + 2.0f * j / tess_th) * turns;
		float th = M_TWOPIf * j / tess_th * turns;
		vec helixpos(cos(th), sin(th), z);
		vec xdir(cos(th), sin(th), 0);
		vec ydir = xform::rot(angle, xdir) * vec(0, 0, -1);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_TWOPIf * i / tess_ph;
			mesh->vertices[i+j*tess_ph] = helixpos +
						      r*cos(ph)*xdir +
						      r*sin(ph)*ydir;
		}
	}
	return mesh;
}


TriMesh *make_sphere_polar(int tess_ph, int tess_th)
{
	if (tess_th < 3)
		tess_th = 3;
	if (tess_ph < 3)
		tess_ph = 3;

	TriMesh *mesh = new TriMesh;
	mesh->vertices.reserve(2+tess_ph*(tess_th-1));

	mkpoint(mesh, 0, 0, -1);
	for (int j = 1; j < tess_th; j++) {
		float th = M_PIf * j / tess_th;
		float z = -cos(th);
		float r = sin(th);
		for (int i = 0; i < tess_ph; i++) {
			float ph = M_TWOPIf * i / tess_ph;
			mkpoint(mesh, r*cos(ph), r*sin(ph), z);
		}
	}
	mkpoint(mesh, 0, 0, 1);

	mesh->faces.reserve(2*tess_th*tess_ph - 2*tess_ph);

	for (int i = 0; i < tess_ph; i++)
		mkface(mesh, 0, ((i+1)%tess_ph)+1, i+1);

	for (int j = 0; j < tess_th-2; j++) {
		int base = 1 + j * tess_ph;
		for (int i = 0; i < tess_ph; i++) {
			int i1 = (i+1)%tess_ph;
			mkquad(mesh, base + i, base + i1,
				base+tess_ph+i, base+tess_ph+i1);
		}
	}

	int base = 1 + (tess_th-2)*tess_ph;
	for (int i = 0; i < tess_ph; i++)
		mkface(mesh, base+i, base+((i+1)%tess_ph), base+tess_ph);

	return mesh;
}


TriMesh *make_platonic(int nfaces)
{
	float phi = 0.5f * (1.0f + sqrt(5.0f));
	float phibar = phi - 1.0f;

	TriMesh *mesh = new TriMesh;
	switch (nfaces) {
		case 4:
			mkpoint(mesh, 0, 0, 3);
			mkpoint(mesh, sqrt(8.0f), 0.0f, -1.0f);
			mkpoint(mesh, sqrt(8.0f) - sqrt(18.0f),
				      sqrt(6.0f), -1.0f);
			mkpoint(mesh, sqrt(8.0f) - sqrt(18.0f),
				     -sqrt(6.0f), -1.0f);
			for (int i = 0; i < 4; i++)
				normalize(mesh->vertices[i]);
			mkface(mesh, 0, 1, 2);
			mkface(mesh, 0, 2, 3);
			mkface(mesh, 3, 2, 1);
			mkface(mesh, 3, 1, 0);
			break;

		case 6:
			mkpoint(mesh, -1, -1, -1);
			mkpoint(mesh, -1,  1, -1);
			mkpoint(mesh,  1,  1, -1);
			mkpoint(mesh,  1, -1, -1);
			mkpoint(mesh,  1, -1,  1);
			mkpoint(mesh,  1,  1,  1);
			mkpoint(mesh, -1,  1,  1);
			mkpoint(mesh, -1, -1,  1);
			for (int i = 0; i < 8; i++)
				normalize(mesh->vertices[i]);
			tess4(mesh, 0, 1, 2, 3);
			tess4(mesh, 5, 4, 3, 2);
			tess4(mesh, 4, 5, 6, 7);
			tess4(mesh, 1, 0, 7, 6);
			tess4(mesh, 6, 5, 2, 1);
			tess4(mesh, 0, 3, 4, 7);
			break;

		case 8:
			mkpoint(mesh, 1, 0, 0);
			mkpoint(mesh, 0, 1, 0);
			mkpoint(mesh, 0, 0, 1);
			mkpoint(mesh, -1, 0, 0);
			mkpoint(mesh, 0, -1, 0);
			mkpoint(mesh, 0, 0, -1);
			mkface(mesh, 0, 1, 2);
			mkface(mesh, 0, 2, 4);
			mkface(mesh, 0, 4, 5);
			mkface(mesh, 0, 5, 1);
			mkface(mesh, 3, 2, 1);
			mkface(mesh, 3, 1, 5);
			mkface(mesh, 3, 5, 4);
			mkface(mesh, 3, 4, 2);
			break;

		case 12:
			mkpoint(mesh, -1, -1, -1);
			mkpoint(mesh, -1,  1, -1);
			mkpoint(mesh,  1,  1, -1);
			mkpoint(mesh,  1, -1, -1);
			mkpoint(mesh,  1, -1,  1);
			mkpoint(mesh,  1,  1,  1);
			mkpoint(mesh, -1,  1,  1);
			mkpoint(mesh, -1, -1,  1);
			mkpoint(mesh,  phi,  phibar, 0);
			mkpoint(mesh,  phi, -phibar, 0);
			mkpoint(mesh, -phi, -phibar, 0);
			mkpoint(mesh, -phi,  phibar, 0);
			mkpoint(mesh,  phibar, 0,  phi);
			mkpoint(mesh, -phibar, 0,  phi);
			mkpoint(mesh, -phibar, 0, -phi);
			mkpoint(mesh,  phibar, 0, -phi);
			mkpoint(mesh, 0,  phi,  phibar);
			mkpoint(mesh, 0,  phi, -phibar);
			mkpoint(mesh, 0, -phi, -phibar);
			mkpoint(mesh, 0, -phi,  phibar);
			for (int i = 0; i < 20; i++)
				normalize(mesh->vertices[i]);
			tess5(mesh, 13, 12, 5, 16, 6);
			tess5(mesh, 12, 13, 7, 19, 4);
			tess5(mesh, 12, 4, 9, 8, 5);
			tess5(mesh, 9, 4, 19, 18, 3);
			tess5(mesh, 2, 17, 16, 5, 8);
			tess5(mesh, 0, 14, 15, 3, 18);
			tess5(mesh, 0, 10, 11, 1, 14);
			tess5(mesh, 0, 18, 19, 7, 10);
			tess5(mesh, 7, 13, 6, 11, 10);
			tess5(mesh, 3, 15, 2, 8, 9);
			tess5(mesh, 1, 17, 2, 15, 14);
			tess5(mesh, 1, 11, 6, 16, 17);
			break;

		case 20:
			mkpoint(mesh, 0, -1,  phi);
			mkpoint(mesh, 0,  1,  phi);
			mkpoint(mesh, 0, -1, -phi);
			mkpoint(mesh, 0,  1, -phi);
			mkpoint(mesh,  phi, 0,  1);
			mkpoint(mesh,  phi, 0, -1);
			mkpoint(mesh, -phi, 0,  1);
			mkpoint(mesh, -phi, 0, -1);
			mkpoint(mesh,  1,  phi, 0);
			mkpoint(mesh,  1, -phi, 0);
			mkpoint(mesh, -1,  phi, 0);
			mkpoint(mesh, -1, -phi, 0);
			for (int i = 0; i < 12; i++)
				normalize(mesh->vertices[i]);
			mkface(mesh, 0, 4, 1);
			mkface(mesh, 0, 9, 4);
			mkface(mesh, 9, 5, 4);
			mkface(mesh, 4, 5, 8);
			mkface(mesh, 4, 8, 1);
			mkface(mesh, 8, 10, 1);
			mkface(mesh, 8, 3, 10);
			mkface(mesh, 5, 3, 8);
			mkface(mesh, 5, 2, 3);
			mkface(mesh, 2, 7, 3);
			mkface(mesh, 7, 10, 3);
			mkface(mesh, 7, 6, 10);
			mkface(mesh, 7, 11, 6);
			mkface(mesh, 11, 0, 6);
			mkface(mesh, 0, 1, 6);
			mkface(mesh, 6, 1, 10);
			mkface(mesh, 9, 0, 11);
			mkface(mesh, 9, 11, 2);
			mkface(mesh, 9, 2, 5);
			mkface(mesh, 7, 2, 11);
			break;

		default:
			TriMesh::eprintf("\nNumber of faces must be 4, 6, 8, 12, or 20\n\n"); 
			delete mesh;
			return 0;
	}
	return mesh;
}


TriMesh *make_sphere_subdiv(int nfaces, int nsubdiv)
{
	TriMesh *mesh = make_platonic(nfaces);
	if (!mesh)
		return 0;
	for (size_t i = 0; i < mesh->vertices.size(); i++)
		normalize(mesh->vertices[i]);
	for (int i = 0; i < nsubdiv; i++)
		subdiv(mesh);
	for (size_t i = 0; i < mesh->vertices.size(); i++)
		normalize(mesh->vertices[i]);
	return mesh;
}


TriMesh *make_rd()
{
	TriMesh *mesh = make_platonic(6);
	for (int i = 8; i < 14; i++)
		mesh->vertices[i] *= 2.0f;
	return mesh;
}

TriMesh *make_rt()
{
	const float phi3 = 0.5f * (5.0f - sqrt(5.0f));
	TriMesh *mesh = make_platonic(12);
	for (int i = 20; i < 32; i++)
		mesh->vertices[i] *= phi3;
	return mesh;
}

TriMesh *make_sor(int tess_th, const char *curve_pts)
{
	FILE *f = fopen(curve_pts, "r");
	if (!f) {
		TriMesh::eprintf("Couldn't open %s\n", curve_pts);
		exit(1);
	}

	vector<point> pts;
	float x, y, z;
	while (1) {
		if (fscanf(f, "%f%f%f", &x, &y, &z) != 3)
			break;
		pts.push_back(point(x,y,z));
	}
	fclose(f);

	int n = pts.size();
	if (n < 2) {
		TriMesh::eprintf("Not enough points read from %s\n", curve_pts);
		exit(1);
	}

	TriMesh *mesh = make_cyl(tess_th, n - 1);

	for (int i = 0; i < tess_th; i++) {
		xform xf = xform::rot(2.0 * M_PI * i / tess_th, 0, 0, 1);
		for (int j = 0; j < n; j++)
			mesh->vertices[i + j * tess_th] = xf * pts[j];
	}

	return mesh;
}

void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s shape outfile\n", myname);
	fprintf(stderr, "Shapes:\n");
	fprintf(stderr, "	plane m [n]	m x n tessellated square (default n = m)\n");
	fprintf(stderr, "	bump n sigma	n x n tessellated Gaussian bump of width sigma\n");
	fprintf(stderr, "	wave n omega	n x n tessellated sine wave of frequency omega\n");
	fprintf(stderr, "	frac n		n x n fractal landscape\n");
	fprintf(stderr, "	cube n		n x n tessellated cube\n");
	fprintf(stderr, "	disc n m	Circular disc, tessellated with m rings of n points\n");
	fprintf(stderr, "	cyl n m [r]	Cylinder of radius r (default 1)\n");
	fprintf(stderr, "	ccyl n m [r]	Capped cylinder\n");
	fprintf(stderr, "	cone n m [r]	Cone\n");
	fprintf(stderr, "	ccone n m [r]	Capped cone\n");
	fprintf(stderr, "	torus n m [r]	Torus of minor radius r (default 0.25)\n");
	fprintf(stderr, "	knot n m [r]	Trefoil knot of minor radius r (default 0.2)\n");
	fprintf(stderr, "	klein n m	Klein bottle\n");
	fprintf(stderr, "	helix n m t [r]	Helix of minor radius r, with t turns\n");
	fprintf(stderr, "	sphere n m	Sphere, tessellated in polar coordinates\n");
	fprintf(stderr, "	platonic n	Platonic solid with n sides\n");
	fprintf(stderr, "	ssphere n m	Sphere, subdivided m times from a Platonic of n sides\n");
	fprintf(stderr, "	rd		Rhombic dodecahedron\n");
	fprintf(stderr, "	rt		Rhombic triacontahedron\n");
	fprintf(stderr, "	sor n curve.txt	Surface of revolution: n copies of curve, rot z axis\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	if (argc < 3)
		usage(argv[0]);

	TriMesh *mesh = NULL;
	const char *outfilename = argv[2];

	if (!strcmp(argv[1], "plane") && argc > 4) {
		mesh = make_plane(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "plane") && argc > 3) {
		mesh = make_plane(atoi(argv[2]), atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "bump") && argc > 4) {
		mesh = make_bump(atoi(argv[2]), ATOF(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "wave") && argc > 4) {
		mesh = make_wave(atoi(argv[2]), ATOF(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "frac") && argc > 3) {
		mesh = make_frac(atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "cube") && argc > 3) {
		mesh = make_cube(atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "disc") && argc > 4) {
		mesh = make_disc(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "disk") && argc > 4) {
		mesh = make_disc(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "cyl") && argc > 5) {
		mesh = make_cyl(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "cyl") && argc > 4) {
		mesh = make_cyl(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "ccyl") && argc > 5) {
		mesh = make_ccyl(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "ccyl") && argc > 4) {
		mesh = make_ccyl(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "cone") && argc > 5) {
		mesh = make_cone(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "cone") && argc > 4) {
		mesh = make_cone(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "ccone") && argc > 5) {
		mesh = make_ccone(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "ccone") && argc > 4) {
		mesh = make_ccone(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "torus") && argc > 5) {
		mesh = make_torus(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "torus") && argc > 4) {
		mesh = make_torus(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "knot") && argc > 5) {
		mesh = make_knot(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "knot") && argc > 4) {
		mesh = make_knot(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "klein") && argc > 4) {
		mesh = make_klein(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "helix") && argc > 6) {
		mesh = make_helix(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]), ATOF(argv[5]));
		outfilename = argv[6];
	} else if (!strcmp(argv[1], "helix") && argc > 5) {
		mesh = make_helix(atoi(argv[2]), atoi(argv[3]), ATOF(argv[4]));
		outfilename = argv[5];
	} else if (!strcmp(argv[1], "sphere") && argc > 4) {
		mesh = make_sphere_polar(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "platonic") && argc > 3) {
		mesh = make_platonic(atoi(argv[2]));
		outfilename = argv[3];
	} else if (!strcmp(argv[1], "ssphere") && argc > 4) {
		mesh = make_sphere_subdiv(atoi(argv[2]), atoi(argv[3]));
		outfilename = argv[4];
	} else if (!strcmp(argv[1], "rd") && argc > 2) {
		mesh = make_rd();
		outfilename = argv[2];
	} else if (!strcmp(argv[1], "rt") && argc > 2) {
		mesh = make_rt();
		outfilename = argv[2];
	} else if (!strcmp(argv[1], "sor") && argc > 4) {
		mesh = make_sor(atoi(argv[2]), argv[3]);
		outfilename = argv[4];
	}

	if (mesh) {
		mesh->need_tstrips();
		mesh->write(outfilename);
	} else {
		usage(argv[0]);
	}

}

