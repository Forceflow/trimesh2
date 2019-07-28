
This is trimesh2, a C++ library and set of utilities for input, output,
and basic manipulation of 3D triangle meshes.  The goals of the code are
ease of use and efficiency, possibly at the expense of some generality.


The library includes the following:

 - Support for reading PLY, OFF, 3DS, VVD, and Wavefront OBJ files,
   together with a few other formats of mostly local interest (SM, RAY).

 - Support for writing PLY, OFF, and OBJ files.

 - A templated C++ class for constant-length vectors, with support for the
   usual arithmetic operations (add, subtract, componentwise multiply and
   divide, dot product, cross product, etc.)

 - A class for rigid-body transformations.

 - An OpenGL trackball/arcball implementation, with automatic selection of
   rotation center.

 - Algorithms for subdivision, smoothing, curvature estimation, triangle
   stripping, and various other simple mesh manipulations.


Bundled together with the library are:

 - Bernd Gartner's "miniball" code for efficient exact bounding sphere
   computation:  http://www.inf.ethz.ch/personal/gaertner/minibal.html

 - An ever-so-slightly tweaked version of the freeglut library, an
   open source GLUT replacement:  http://freeglut.sourceforge.net/

 - GLUI, a user interface library that sits on top of OpenGL, hence is
   portable to many different systems:  http://www.nigels.com/glt/glui/


In addition, the following utility programs are included:

 - mesh_view: a simple 3D mesh viewer

 - mesh_make: create arbitrarily-tessellated meshes of various simple shapes

 - mesh_filter: applies a variety of simple transformations to a mesh,
   such as converting formats, flipping faces, subdivision, smoothing,
   rigid-body transformations, etc.

 - mesh_cc: list and/or extract connected components from a mesh

 - mesh_cat: combine several meshes into a single file

 - mesh_align: align 2 meshes using ICP


The author of trimesh2 is Szymon Rusinkiewicz, and the library is
distributed under the GNU General Public License (GPL).  The various
libraries distributed with trimesh2 are open source, and are believed
to be covered by GPL-compatible licenses.  Please see the COPYING file.


Running the utilities:
----------------------

All of the programs must be run from the command line.  Each displays
a list of options when run without arguments.

In the mesh_view program, the view is changed by dragging with particular
mouse buttons as follows:

        Left:        Rotate (release while moving to spin the object)
        Middle:      Translate left, right, up, down
        Left+right:  Translate left, right, up, down
        Right:       Translate forward and back
        Mouse wheel: Translate forward and back
        Ctrl+Left:   Move the light

In addition, the following keys invoke various options:

        Space bar:   Reset to initial view
        e:           Toggle display of mesh edges
        l:           Toggle lighting
        f:           Toggle false-color rendering
        q:           Exit
        s:           Toggle specular component in reflection
        Shift-2:     Toggle 2-sided lighting
        x:           Write current viewpoint as a ".xf" file
        1, 2, etc.:  Toggle the n-th mesh off and on


Using the library:
------------------

There is no complete documentation for the library - poke around in
include/*.h to see what is available, and look in utilsrc/*.cc to
see how it's used.  Here's a trivial program to help you get started:


#include "TriMesh.h"
#include <cstdlib>
#include <iostream>
using namespace trimesh;
using namespace std;

int main(int argc, char *argv[])
{
        const char *filename = "foo.ply";
        TriMesh *m = TriMesh::read(filename);
        if (!m)
                exit(1);

        cout << "There are " << m->vertices.size() << " vertices" << endl;
        cout << "Vertex 0 is at " << m->vertices[0] << endl;

        // Convert triangle strips to faces, if necessary
        m->need_faces();
        cout << "Face 0 has vertices " << m->faces[0][0] << ", "
             << m->faces[0][1] << ", and " << m->faces[0][2] << endl;

        m->need_normals();
        cout << "Vertex 0 has normal " << m->normals[0] << endl;
}


Compiling:
----------

The code is written in C++, and is known to compile using a recent version
of g++ (4.0.1 or later, but avoid 4.1.2).  Compiling under Windows is
possible using Cygwin or MinGW32.  Other compilers have not been tested -
let me know of any successes or failures.

After unpacking the code, edit the appropriate "Makedefs.*" file to
specify the compiler and any necessary options.  If there is no Makedefs
for your OS, run "uname" to determine what it should be named and create
one (probably starting with one of the existing files as a template).
At this point, just running "make" in the top-level directory should make
the library (placed in lib.`uname`) and utilities (placed in bin.`uname`).

If you have any problems, please double check the following:
 - You are using a recent version of the compiler
 - You are using a recent version of GNU make
 - The correct -I option is included to find all the headers you need,
   including OpenGL headers
If you still have problems, let me (smr@princeton.edu) know and I'll do
my best to help resolve them.


