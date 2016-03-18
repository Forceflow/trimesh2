# TriMesh2

**TriMesh2: ** C++ library and set of utilities for input, output, and basic manipulation of 3D triangle meshes

This is a fork of the trimesh2 library (originally by [Szymon Rusinkiewicz](https://www.cs.princeton.edu/~smr/)), which I use a lot in my other graphics projects. I'm staying close to the original, and adding my fixes for MSVC compilation.

 * For the Trimesh2 project, see [the original Trimesh2 homepage](http://gfx.cs.princeton.edu/proj/trimesh2/).
 * Original MSVC 2012 project by Bengt Rosenberger.

Notable changes/additions to vanilla trimesh2
 * Fixes for MSVC 12.0 / 14.0 compilation
 * MSVC project for Visual Studio 2013 professional
 * MSVC project for Visual Studio Community Edition 2015
 * Fixes for OpenMP compilation in VS 2015
 * Fixes for FreeGlut / Gluit compilation in VS 2015
 * Added (experimental) 64-bit MSVC compilation support

Todo
 * Build trimesh2 tools (mesh_view, etc.) using MSVC project
