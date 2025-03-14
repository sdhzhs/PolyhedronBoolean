#!/bin/sh
compare() {
  echo "compare $1:"
  diff $1 golden/$1
}

if [ -z "${MACRO+x}" ]; then
  compare "cgal_surf.out"
  compare "cgal_mesh.out"
  compare "cgal_mesh.off"
  compare "cgal_plane.out"
  compare "cgal_surf_union.out"
  compare "cgal_surf_union_taper.out"
  compare "cgal_vert_union_taper.out"
  compare "cgal_surf_div.out"
  compare "cgal_mesh_cube.out"
  compare "cgal_mesh_layout.out"
  compare "cgal_surf_layout.out"
  compare "cgal_nef_layout.out"
  compare "cgal_mesh_layout.off"
fi
if [ -n "${MACRO+x}" ]; then
  compare "cgal_mesh_interplane.off"
  compare "cgal_surf_interplane.out"
  compare "cgal_nef_interplane.out"
fi
