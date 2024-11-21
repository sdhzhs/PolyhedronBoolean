# PolyhedronBoolean
Background

The C++ codes in this repository are extracted and extended from one pre-process library of a Computer Aided Engineering (CAE) project. The main goal is to construct three-dimensional (3D) polyhedra from two-dimensional (2D) polygons according to some technology rules, and make boolean operations among constructed three-dimensional polyhedra or between polyhedra and elementary geometries (cube, block...).  The codes are deeply based on one package in Computational Geometry Algorithm Library (CGAL), which is called Nef_Polyhedron. Thus this repository actually belongs to one kind of middle-level library which calls data structures and algorithms in CGAL as low-level interfaces. The main feature in this repository is that it includes supplementary data structures and algorithms which constitute complete steps of 3D boolean operations. The codes in this repository can be used in addtion to the CGAL engine which makes the feasible implemetation of geometric pre-process for CAE developers, and build a bridge between commonly used geometric objects and CGAL objects which have deep mathematical background.

Data structure

The concepts of Nef_Polyhedron are very suitable for the target of doing three-dimensional boolean operations between two polyhedra. The Nef_Polyhedron can be considered as mathematical extension of commonly used 3D geometric concepts - polyhedron. It can treat many corner and degenerated scenarios when doing 3D boolean operations. The data structure of Nef_Polyhedron includes sets of half-faces, half-edges and vertices, please refer to the documentation of CGAL for detail of organization of these half-kind geometric elements. CGAL provides sufficient iterators to traverse these geometric elements without the knowledge of internal data structure in Nef_Polyhedron.

There are supplementary data structures in this repository which represent surface mesh to build the polyhedra (inputs) or extracted faces/vertices lists (outputs). Only the C++ STL containers instantiated by 3D double-precision points or CGAL internal points are used to represent the additional data structures so that they are easy to understand for a C++ developer.
 
Algorithm

The operations of predication and construction are important to do 3D boolean operations between two polyhedra. Note that there are some geometic kernels in CGAL which can be selected with the Nef_Polyhedron package. The most widely used one is called exact predicates and exact constructions kernel, which is suitable for input of Cartesian coordinates using double-precision floating number. The "exact" method in constructions of points among linear geometric objects is implemented using multi-precision library which can do exact arithmetic for rational number without overflow and accuracy loss. It is the way CGAL uses to overcome issues caused by geometric tolerance. Refer to documentations of CGAL for other geometic kernels which can be used in Nef_Polyhedron package and the detail of predication and construction when doing 3D boolean operations between two Nef_Polyhedra, even though the documentation is written by the authors very mathematically.

The typical steps of a 3D boolean operation include
1. Construction of original Nef_Polyhedra
2. 3D boolean operations on Nef_Polyhedra
3. Extraction of B-rep faces in result Nef_Polyhedron to supplementary data structure

The ideas of original Nef_Polyhedra construction include the way from a complete surface mesh and the way using concepts of Constructive Solid Geometry (CSG, built from multiple isolated faces or planes). Please refer to comments in the test codes of this repository.

Usage

There is a hierarchical Makefile system in this repository, just type

make

in a Unix-like OS environment in which g++ has been deployed, then an exec named "CGALNefOp" will be generated. Moreover, type

make clean

can clean all executable and object files.

The codes in source file Main.cpp represent some test examples which call middle-level interfaces in this repository. The examples include
1.  Test division method of rational numbers from different multi-precision libraries
2.  Test OFF_to_nef_3 interface in Nef_Polyhedron CGAL Package
3.  Test a complete boolean operation process of two simple cubes. Boolean operation between a Nef_Polyhedron and a plane is also tested
4.  Test union boolean operations between Nef_Polyhedra which are not manifold (single or multiple faces which are not closed)
5.  Test union boolean operations between Nef_Polyhedra which are not manifold (single or multiple faces which are not closed), this is a degenerated case compared to the above example
6.  Test the consistency between self-defined surfaces file format and built-in Nef_Polyhedron file format, some mesh conversion functions are also tested
7.  Test the possibility of generating a 3D complex from 2D layout polygons
8.  Test intersection boolean operations between Nef_Polyhedra which are planes which represent infinite half spaces, this is another way to construct 3D closed polyhedra


