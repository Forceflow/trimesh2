# TriMesh2 (2.13)
[![Build Status](https://travis-ci.org/Forceflow/trimesh2.svg?branch=master)](https://travis-ci.org/Forceflow/trimesh2) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A C++ library and set of utilities for input, output, and basic manipulation of 3D triangle meshes.

![trimesh2 logo](https://raw.githubusercontent.com/Forceflow/trimesh2/master/html/trimesh_logo.jpg)

This is a fork of the [TriMesh2 library](http://gfx.cs.princeton.edu/proj/trimesh2/) (originally by [Szymon Rusinkiewicz](https://www.cs.princeton.edu/~smr/)), which I use a lot in my other graphics projects. I like TriMesh2 because of the low setup costs required to do model loading, as well as the robust and powerful implementation of model manipulation.

The original TriMesh2 project is quite Linux/GCC-oriented, and only has limited Win32 support (through MinGW compilation targets). The primary aim of this fork is to add a stable Visual Studio solution, for both x86 and x64 targets.
 
## Getting started
 * Download a [prebuilt release](https://github.com/Forceflow/trimesh2/releases) of trimesh2 or build the library yourself. The static library will be called `trimesh.lib`, the debug version is `trimeshd.lib`.
  * For **Windows**, build solutions for VS2017 and VS2019 are provided in the `mscv`folder. All solutions work in the free Community Editions of VS.
    * The built library will be placed in a folder named `lib.(architecture).(visual studio version)`. For example, for a 64-bit Visual Studio 2017 build, it will be `lib.win64.vs141`.
    * The builds of the utilities will be placed in a folder named `util.(architecture).(visual studio version)`. For example, for a 32-bit Visual Studio 2019, they will be in `util.win32.vs142`.
   * For **Linux**, a makefile is provided. You might need additional packages before you can build on your system. On Ubuntu these are: `mesa-common-dev libglu1-mesa-dev libxi-dev`.
   * I'm being told it builds on **OSX** using the provided makefile too, but I have no way to check. If you encounter problems, please, file an issue report :)
 * In your own project, make sure you include the header `include/TriMesh.h`, and make sure the static library is in your build path. All Trimesh2 functions will be in the `TriMesh` namespace.
 * For inspiration on how to use the library, check out the utilities in the `utilsrc` folder, or just start by loading a model : `TriMesh* themesh = TriMesh::read(filename);`

## Info
Legacy:

 * For the original TriMesh2 project, see [the Trimesh2 homepage](http://gfx.cs.princeton.edu/proj/trimesh2/).
 * Original MSVC 2012 project by Bengt Rosenberger.
 
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

Notable changes/additions to vanilla trimesh2
 * Fixes for MSVC 14.0 compilation
 * MSVC project for Visual Studio Community Edition 2017
 * Fixes for FreeGlut / Gluit compilation in VS 2017
 * Added (experimental) 64-bit MSVC compilation support

For todo/planned features, see todo.txt.
