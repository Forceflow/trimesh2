/*
Szymon Rusinkiewicz
Princeton University

mesh_view.cc
Simple viewer
*/

#include "TriMesh.h"
#include "XForm.h"
#include "GLCamera.h"
#include "ICP.h"
#include "strutil.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <GL/glut.h>
using namespace std;
using namespace trimesh;


// Globals
vector<TriMesh *> meshes;
vector<xform> xforms;
vector<bool> visible;
vector<string> filenames;

TriMesh::BSphere global_bsph;
xform global_xf;
GLCamera camera;

int current_mesh = -1;

bool draw_edges = false;
bool draw_points = false;
bool draw_2side = false;
bool draw_shiny = true;
bool draw_lit = true;
bool draw_falsecolor = false;
bool draw_index = false;
bool white_bg = false;
bool grab_only = false;
int point_size = 1, line_width = 1;


// Signal a redraw
void need_redraw()
{
	glutPostRedisplay();
}


// Clear the screen
void cls()
{
	glDisable(GL_DITHER);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	if (white_bg)
		glClearColor(1, 1, 1, 0);
	else
		glClearColor(0.08f, 0.08f, 0.08f, 0);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}


// Set up lights and materials
void setup_lighting(int id)
{
	Color c(1.0f);
	if (draw_falsecolor)
		c = Color::hsv(-3.88f * id, 0.6f + 0.2f * sin(0.42f * id), 1);
	glColor3fv(c);

	if (!draw_lit || meshes[id]->normals.empty()) {
		glDisable(GL_LIGHTING);
		return;
	}

	GLfloat mat_specular[4] = { 0.18f, 0.18f, 0.18f, 0.18f };
	if (!draw_shiny) {
		mat_specular[0] = mat_specular[1] =
		mat_specular[2] = mat_specular[3] = 0.0f;
	}
	GLfloat mat_shininess[] = { 64 };
	GLfloat global_ambient[] = { 0.02f, 0.02f, 0.05f, 0.05f };
	GLfloat light0_ambient[] = { 0, 0, 0, 0 };
	GLfloat light0_diffuse[] = { 0.85f, 0.85f, 0.8f, 0.85f };
	if (current_mesh >= 0 && id != current_mesh) {
		light0_diffuse[0] *= 0.5f;
		light0_diffuse[1] *= 0.5f;
		light0_diffuse[2] *= 0.5f;
	}
	GLfloat light1_diffuse[] = { -0.01f, -0.01f, -0.03f, -0.03f };
	GLfloat light0_specular[] = { 0.85f, 0.85f, 0.85f, 0.85f };
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, draw_2side);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
}


// Draw triangle strips.  They are stored as length followed by values.
void draw_tstrips(const TriMesh *themesh)
{
	static bool use_glArrayElement = false;
	static bool tested_renderer = false;
	if (!tested_renderer) {
		use_glArrayElement = !!strstr(
			(const char *) glGetString(GL_RENDERER), "Intel");
		tested_renderer = true;
	}

	const int *t = &themesh->tstrips[0];
	const int *end = t + themesh->tstrips.size();
	if (use_glArrayElement) {
		while (likely(t < end)) {
			glBegin(GL_TRIANGLE_STRIP);
			int striplen = *t++;
			for (int i = 0; i < striplen; i++)
				glArrayElement(*t++);
			glEnd();
		}
	} else {
		while (likely(t < end)) {
			int striplen = *t++;
			glDrawElements(GL_TRIANGLE_STRIP, striplen, GL_UNSIGNED_INT, t);
			t += striplen;
		}
	}
}


// Draw the mesh
void draw_mesh(int i)
{
	const TriMesh *themesh = meshes[i];

	glPushMatrix();
	glMultMatrixd(xforms[i]);

	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);

	if (draw_2side) {
		glDisable(GL_CULL_FACE);
	} else {
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
	}

	// Vertices
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT,
			sizeof(themesh->vertices[0]),
			&themesh->vertices[0][0]);

	// Normals
	if (!themesh->normals.empty() && !draw_index) {
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT,
				sizeof(themesh->normals[0]),
				&themesh->normals[0][0]);
	} else {
		glDisableClientState(GL_NORMAL_ARRAY);
	}

	// Colors
	if (!themesh->colors.empty() && !draw_falsecolor && !draw_index) {
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_FLOAT,
			       sizeof(themesh->colors[0]),
			       &themesh->colors[0][0]);
	} else {
		glDisableClientState(GL_COLOR_ARRAY);
	}

	// Main drawing pass
	if (draw_points || themesh->tstrips.empty()) {
		// No triangles - draw as points
		glPointSize(float(point_size));
		glDrawArrays(GL_POINTS, 0, themesh->vertices.size());
		glPopMatrix();
		return;
	}

	if (draw_edges) {
		glPolygonOffset(10.0f, 10.0f);
		glEnable(GL_POLYGON_OFFSET_FILL);
	}

	draw_tstrips(themesh);
	glDisable(GL_POLYGON_OFFSET_FILL);

	// Edge drawing pass
	if (draw_edges) {
		glPolygonMode(GL_FRONT, GL_LINE);
		glLineWidth(float(line_width));
		glDisableClientState(GL_COLOR_ARRAY);
		glDisable(GL_COLOR_MATERIAL);
		GLfloat global_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
		GLfloat light0_diffuse[] = { 0.8f, 0.8f, 0.8f, 0.0f };
		GLfloat light1_diffuse[] = { -0.2f, -0.2f, -0.2f, 0.0f };
		GLfloat light0_specular[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
		GLfloat mat_diffuse[4] = { 0.0f, 0.0f, 1.0f, 1.0f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		glColor3f(0, 0, 1); // Used iff unlit
		draw_tstrips(themesh);
		glPolygonMode(GL_FRONT, GL_FILL);
	}

	glPopMatrix();
}


// Draw the scene
void redraw()
{
	timestamp t = now();
	camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
	glPushMatrix();
	glMultMatrixd(global_xf);
	cls();
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		setup_lighting(i);
		draw_mesh(i);
	}

	glPopMatrix();
	glutSwapBuffers();
	if (grab_only) {
		void dump_image();
		dump_image();
		exit(0);
	}
	printf("\r                        \r%.1f msec.", 1000.0f * (now() - t));
	fflush(stdout);
}


// Update global bounding sphere.
void update_bsph()
{
	point boxmin(1e38f, 1e38f, 1e38f);
	point boxmax(-1e38f, -1e38f, -1e38f);
	bool some_vis = false;
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])	
			continue;
		some_vis = true;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		for (int j = 0; j < 3; j++) {
			boxmin[j] = min(boxmin[j], c[j]-r);
			boxmax[j] = max(boxmax[j], c[j]+r);
		}
	}
	if (!some_vis)
		return;
	point &gc = global_bsph.center;
	float &gr = global_bsph.r;
	gc = 0.5f * (boxmin + boxmax);
	gr = 0.0f;
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])	
			continue;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		gr = max(gr, dist(c, gc) + r);
	}
}


// Set the view...
void resetview()
{
	camera.stopspin();

	// Reload mesh xforms
	for (size_t i = 0; i < meshes.size(); i++)
		if (!xforms[i].read(xfname(filenames[i])))
			xforms[i] = xform();

	update_bsph();

	// Set camera to first ".camxf" if we have it...
	for (size_t i = 0; i < filenames.size(); i++) {
		if (global_xf.read(replace_ext(filenames[i], "camxf")))
			return;
	}

	// else default view
	global_xf = xform::trans(0, 0, -5.0f * global_bsph.r) *
		    xform::trans(-global_bsph.center);
}


// Make some mesh current
void set_current(int i)
{
	camera.stopspin();
	if (i >= 0 && i < (int)meshes.size() && visible[i])
		current_mesh = i;
	else
		current_mesh = -1;
	need_redraw();
}


// Change visiblility of a mesh
void toggle_vis(int i)
{
	if (i >= 0 && i < (int)meshes.size())
		visible[i] = !visible[i];
	if (current_mesh == i && !visible[i])
		set_current(-1);
	update_bsph();
	need_redraw();
}


// Save the current image to a PPM file.
// Uses the next available filename matching filenamepattern
void dump_image()
{
	// Find first non-used filename
	const char filenamepattern[] = "img%d.ppm";
	int imgnum = 0;
	FILE *f;
	while (1) {
		char filename[1024];
		sprintf(filename, filenamepattern, imgnum++);
		f = fopen(filename, "rb");
		if (!f) {
			f = fopen(filename, "wb");
			printf("\n\nSaving image %s... ", filename);
			fflush(stdout);
			break;
		}
		fclose(f);
	}

	// Read pixels
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	GLint width = V[2], height = V[3];
	char *buf = new char[width*height*3];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n#\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;

	printf("Done.\n\n");
}


// Save scan transforms
void save_xforms()
{
	for (size_t i = 0; i < xforms.size(); i++) {
		string xffile = xfname(filenames[i]);
		printf("Writing %s\n", xffile.c_str());
		xforms[i].write(xffile);
	}
}


// Save camera xform
void save_cam_xform()
{
	std::string camfile = replace_ext(filenames[0], "camxf");
	printf("Writing %s\n", camfile.c_str());
	global_xf.write(camfile);
}


// ICP
void do_icp(int n)
{
	camera.stopspin();

	if (current_mesh < 0 || current_mesh >= (int)meshes.size())
		return;
	if (n < 0 || n >= (int)meshes.size())
		return;
	if (!visible[n] || !visible[current_mesh] || n == current_mesh)
		return;
	ICP(meshes[n], meshes[current_mesh], xforms[n], xforms[current_mesh], 2);
	update_bsph();
	need_redraw();
}


// Handle mouse button and motion events
static unsigned buttonstate = 0;

void doubleclick(int button, int x, int y)
{
	// Render and read back ID reference image
	camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
	glDisable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	draw_index = true;
	glPushMatrix();
	glMultMatrixd(global_xf);
	for (size_t i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		glColor3ub((i >> 16) & 0xff,
			   (i >> 8)  & 0xff,
			    i        & 0xff);
		draw_mesh(i);
	}
	glPopMatrix();
	draw_index = false;
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	y = int(V[1] + V[3]) - 1 - y;
	unsigned char pix[3];
	glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, pix);
	int n = (pix[0] << 16) + (pix[1] << 8) + pix[2];

	if (button == 0 || buttonstate == (1 << 0)) {
		// Double left click - select a mesh
		set_current(n);
	} else if (button == 2 || buttonstate == (1 << 2)) {
		// Double right click - ICP current to clicked-on
		do_icp(n);
	}
}

void mousehelperfunc(int x, int y)
{
	static const Mouse::button physical_to_logical_map[] = {
		Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
		Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
	};

	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else if (buttonstate & (1 << 30))
		b = Mouse::LIGHT;
	else
		b = physical_to_logical_map[buttonstate & 7];

	if (current_mesh < 0) {
		camera.mouse(x, y, b,
			     global_xf * global_bsph.center, global_bsph.r,
			     global_xf);
	} else {
		xform tmp_xf = global_xf * xforms[current_mesh];
		camera.mouse(x, y, b,
			     tmp_xf * meshes[current_mesh]->bsphere.center,
			     meshes[current_mesh]->bsphere.r,
			     tmp_xf);
		xforms[current_mesh] = inv(global_xf) * tmp_xf;
		update_bsph();
	}
}

void mousemotionfunc(int x, int y)
{
	mousehelperfunc(x,y);
	if (buttonstate)
		need_redraw();
}

void mousebuttonfunc(int button, int state, int x, int y)
{
	static timestamp last_click_time;
	static unsigned last_click_buttonstate = 0;
	static float doubleclick_threshold = 0.4f;

	if (glutGetModifiers() & GLUT_ACTIVE_CTRL)
		buttonstate |= (1 << 30);
	else
		buttonstate &= ~(1 << 30);

	if (state == GLUT_DOWN) {
		buttonstate |= (1 << button);
		if (buttonstate == last_click_buttonstate &&
		    now() - last_click_time < doubleclick_threshold) {
			doubleclick(button, x, y);
			last_click_buttonstate = 0;
		} else {
			last_click_time = now();
			last_click_buttonstate = buttonstate;
		}
	} else {
		buttonstate &= ~(1 << button);
	}

	mousehelperfunc(x, y);
	if (buttonstate & ((1 << 3) | (1 << 4))) // Wheel
		need_redraw();
	if ((buttonstate & 7) && (buttonstate & (1 << 30))) // Light
		need_redraw();
}


// Idle callback
void idle()
{
	xform tmp_xf = global_xf;
	if (current_mesh >= 0)
		tmp_xf = global_xf * xforms[current_mesh];

	if (camera.autospin(tmp_xf))
		need_redraw();
	else
		usleep(10000);

	if (current_mesh >= 0) {
		xforms[current_mesh] = inv(global_xf) * tmp_xf;
		update_bsph();
	} else {
		global_xf = tmp_xf;
	}
}


// Keyboard
#define Ctrl (1-'a')
void keyboardfunc(unsigned char key, int, int)
{
	switch (key) {
		case ' ':
			if (current_mesh < 0)
				resetview();
			else
				set_current(-1);
			break;
		case '@': // Shift-2
			draw_2side = !draw_2side; break;
		case 'e':
			draw_edges = !draw_edges; break;
		case 'f':
			draw_falsecolor = !draw_falsecolor; break;
		case 'l':
			draw_lit = !draw_lit; break;
		case 'p':
			point_size++; break;
		case 'P':
			point_size = max(1, point_size - 1); break;
		case Ctrl+'p':
			draw_points = !draw_points; break;
		case 's':
			draw_shiny = !draw_shiny; break;
		case 't':
			line_width++; break;
		case 'T':
			line_width = max(1, line_width - 1); break;
		case 'w':
			white_bg = !white_bg; break;
		case 'I':
			dump_image(); break;
		case Ctrl+'x':
			save_xforms();
			break;
		case Ctrl+'c':
			save_cam_xform();
			break;
		case '\033': // Esc
		case Ctrl+'q':
		case 'Q':
		case 'q':
			exit(0);
		case '0':
			toggle_vis(9); break;
		case '-':
			toggle_vis(10); break;
		case '=':
			toggle_vis(11); break;
		default:
			if (key >= '1' && key <= '9') {
				int m = key - '1';
				toggle_vis(m);
			}
	}
	need_redraw();
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-grab] infile...\n", myname);
	exit(1);
}

int main(int argc, char *argv[])
{
	glutInitWindowSize(512, 512);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInit(&argc, argv);

	if (argc < 2)
		usage(argv[0]);

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-grab")) {
			grab_only = true;
			continue;
		}
		const char *filename = argv[i];
		TriMesh *themesh = TriMesh::read(filename);
		if (!themesh)
			usage(argv[0]);
		themesh->need_normals();
		themesh->need_tstrips();
		themesh->need_bsphere();
		meshes.push_back(themesh);
		xforms.push_back(xform());
		visible.push_back(true);
		filenames.push_back(filename);
	}

	glutCreateWindow(argv[1]);
	glutDisplayFunc(redraw);
	glutMouseFunc(mousebuttonfunc);
	glutMotionFunc(mousemotionfunc);
	glutKeyboardFunc(keyboardfunc);
	glutIdleFunc(idle);

	resetview();
	glutMainLoop();
}


