# TriMesh2 (v2.16)
[![C/C++ CI](https://github.com/Forceflow/trimesh2/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/Forceflow/trimesh2/actions/workflows/c-cpp.yml) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A C++ library and set of utilities for input, output, and basic manipulation of 3D triangle meshes.

![trimesh2 logo](https://raw.githubusercontent.com/Forceflow/trimesh2/main/html/trimesh_logo.jpg)

This is a fork of the [TriMesh2 library](http://gfx.cs.princeton.edu/proj/trimesh2/) (originally by [Szymon Rusinkiewicz](https://www.cs.princeton.edu/~smr/)), which I use a lot in my [other](https://github.com/Forceflow/ooc_svo_builder) [graphics](https://github.com/Forceflow/cuda_voxelizer) [projects](https://github.com/Forceflow/gpu_suggestive_contours). I like TriMesh2 because of the low setup costs required to do model loading, as well as the robust and powerful implementation of various model manipulation techniques.

The original TriMesh2 project is quite Linux/GCC-oriented, and only has limited Win32 support (through MinGW compilation targets). The primary aim of this fork is to add a stable Visual Studio solution, for both x86 and x64 targets, whilst staying as close as possible to the original codebase (and subsequent updates).
 
## Getting started
 * Download a [prebuilt release](https://github.com/Forceflow/trimesh2/releases) of trimesh2 or build the library yourself (see further below)
   * The static library will be called `trimesh.lib`, the debug version is `trimeshd.lib`.
 * Include the header `include/TriMesh.h`, and make sure the static library is in your build path. All Trimesh2 functions will be in the `TriMesh` namespace.
 * Typical ways to get started:
   * Loading a model : `TriMesh* themesh = TriMesh::read(filename);`.
   * This mesh class contains a data member `vertices` which will be filled with all the vertices of your model, and a data member `faces`, which will tell you which vertices make up a face.
   * If your model contains vertex normals, they will be in `normals`. You can (re)compute them by calling `need_normals` on the mesh. There's also `need_bbox` for a bounding box, `need_dcurv` for curvature, etc.
   * For inspiration on how to use the library and its various features, check out `include/TriMesh.h` and the utilities in the `utilsrc` folder.
## Building
### Dependencies
 * The library itself has no dependencies other than the standard C++ STL
 * The (optional) ``mesh_view`` utility has dependencies on OpenGL and Freeglut, because it needs to display a window with a textured model

### Build steps
  * For **Windows**, build solutions for VS2017 and VS2019 are provided in the `mscv`folder, verified working with the [free Community Editions](https://visualstudio.microsoft.com/vs/community/) of Visual Studio. The solutions contain both Debug and Release profiles for 32-bit and 64-bit builds.
    * The built libraries will be placed in a folder named `lib.(architecture).(visual studio version)` in the trimesh2 root folder. For example, for a 64-bit Visual Studio 2017 build, it will be `lib.win64.vs141`. The utilities will be placed in `util.(architecture).(visual studio version)`. This naming scheme is in place to avoid clashing trimesh2 versions.
    * For **Linux**, a makefile is provided. You might need additional packages before you can build the utilities on your system. On Ubuntu these are: `mesa-common-dev libglu1-mesa-dev libxi-dev`.
   * For **OSX**, I'm being told it builds using the makefile, but I have no way to check. If you encounter problems, please, file an issue report :)

## Info
For the original TriMesh2 project, see [the Trimesh2 homepage](http://gfx.cs.princeton.edu/proj/trimesh2/).
 
Features: 

 * Support for reading/writing PLY, OFF, OBJ files. Read-only: 3DS, SM, RAY.
 * Vec: a templated C++ class for constant-length vectors, with support for the usual arithmetic operations and XForm: a class for rigid-body transformations.
 * An OpenGL trackball/arcball implementation, with automatic selection of rotation center.
 * Algorithms for subdivision, smoothing, curvature estimation, triangle stripping, and various other simple mesh manipulations.

The following utility programs are included:

 * **mesh_view**: A simple 3D mesh viewer
 * **mesh_make**: Create arbitrarily-tessellated meshes of various simple shapes
 * **mesh_filter**: Applies a variety of simple transformations to a mesh, such as converting formats, flipping faces, subdivision, smoothing, rigid-body transformations, etc.
 * **mesh_cc**: List and/or extract connected components from a mesh
 * **mesh_cat**: Combine several meshes into a single file
 * **mesh_align**: Align 2 meshes using ICP
 * **mesh_shade**: A few procedural shaders for adding per-vertex color
 * **mesh_check**: Check for some kinds of topological oddities (e.g., more than 2 faces at an edge) in a mesh file.
 * **mesh_crunch**: Quick-n-dirty mesh decimation using the Rossignac-Borrel method of vertex collapse
 * **mesh_info**: Print out some information about a mesh
 * **xf**: Create or compose transformations in .xf files

## Fork Details

This fork stays as close as possible to the original [trimesh2](http://gfx.cs.princeton.edu/proj/trimesh2/) code, only changing the actual source files when a solution for compilation errors cannot be reached through VS pre-build steps or preprocessor magic.

Notable changes compared to vanilla trimesh2
 * MSVC project for Visual Studio Community Edition 2022, 2019 and 2017
 * Fixes for FreeGlut / Gluit compilation
 * Fixes for wingetopt replacement in MSVC
 * Added 64-bit MSVC compilation support

## See Also

Other software for importing and manipulating 3D models:
 * [Tiny OBJ Loader](https://github.com/syoyo/tinyobjloader) (by @syoyo)
 * [Tiny PLY](https://github.com/ddiakopoulos/tinyply) (by @ddiakopoulos)
 * [Open Asset Import Library](http://www.assimp.org/) (ASSIMP)
