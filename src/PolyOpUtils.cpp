#include <cstring>
#include <algorithm>
#include <limits>
#include "PolyOpUtils.h"

#include <boost/foreach.hpp>

using namespace std;

namespace cgalutil {

void CG_tools::createBasePolyhedron(vector<vector<int> > &polyIndexs , const vector<point3d> &points , Polyhedron& p)
{
  CGAL_Build_PolySet<Polyhedron> builder(polyIndexs, points);
  p.delegate(builder);
  if(!(p.is_closed())){
    cerr<<" this polyhedron is not closed "<<endl;
  }
  p.normalize_border();
  return;
}

void CG_tools::createNefPolyhedron(vector<vector<int> > &polyIndexs , const vector<point3d> &points , Nef_polyhedron& poly)
{
  Polyhedron p;
  createBasePolyhedron(polyIndexs, points, p);
  poly = Nef_polyhedron(p);
  return;
}

bool CG_tools::getNefPolyhedronProp(Nef_polyhedron& poly, polyprop& nefprop)
{
  if(poly.is_empty())
    return false;

  nefprop.volnum = poly.number_of_volumes();
  nefprop.facenum = poly.number_of_facets();
  nefprop.edgenum = poly.number_of_edges();
  nefprop.vertnum = poly.number_of_vertices();

  return true;
}

bool CG_tools::getNefPolyhedronVertex(Nef_polyhedron& poly, vector<point3d>& vertices)
{
  vertices.clear();

  if(poly.is_empty())
    return false;

  vector<CGAL_point> outputs;

  for(Nef_vertex_const_iterator vertex = poly.vertices_begin(); vertex != poly.vertices_end(); ++vertex) {
    const CGAL_point& cpt = vertex->point();
    outputs.push_back(cpt);
  }

  convertCgalPointTo3d(outputs, vertices);

  if(vertices.empty())
    return false;
  else
    return true;
}

bool CG_tools::getPolyhedronBoundBox(Nef_polyhedron& poly, boundbox3d& boundcoord)
{
  double xmin, xmax, ymin, ymax, zmin, zmax;
  vector<point3d> vertPts;

  if(poly.is_empty())
    return false;

  bool div = getNefPolyhedronVertex(poly, vertPts);
  if(!div) {
    cerr<<"Divide no vertex!"<<endl;
    return false;
  }

  double DBLMAX = std::numeric_limits<double>::max();
  xmin = DBLMAX;
  xmax = -DBLMAX;
  ymin = DBLMAX;
  ymax = -DBLMAX;
  zmin = DBLMAX;
  zmax = -DBLMAX;

  for(uint j=0; j!=vertPts.size(); ++j)
  {
    if(vertPts[j]._xx<xmin) xmin=vertPts[j]._xx;
    if(vertPts[j]._xx>xmax) xmax=vertPts[j]._xx;
    if(vertPts[j]._yy<ymin) ymin=vertPts[j]._yy;
    if(vertPts[j]._yy>ymax) ymax=vertPts[j]._yy;
    if(vertPts[j]._zz<zmin) zmin=vertPts[j]._zz;
    if(vertPts[j]._zz>zmax) zmax=vertPts[j]._zz;
  }

  boundcoord = boundbox3d(xmin, xmax, ymin, ymax, zmin, zmax);

  return true;
}

bool CG_tools::getNefPolyhedronFace(Nef_polyhedron& nefPoly, vector<vector<point3d>>& boundarys, vector<vector<vector<point3d>>>& holes)
{
  boundarys.clear();
  holes.clear();

  if(nefPoly.is_empty())
    return false;

  for(Nef_half_const_iterator face = nefPoly.halffacets_begin() ; face != nefPoly.halffacets_end() ; face++){
    Volume_const_handle vlo = face->incident_volume();
    if(!(vlo->mark())) continue;
    Nef_half_facet_const_iterator fci;
    bool first(true);
    vector<vector<point3d>> holesin;
    holesin.clear();
    for(fci = face->facet_cycles_begin(); fci != face->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()){
        SHalfedge_const_handle csedge = SHalfedge_const_handle(fci);
        SHalfedge_const_handle sh_end(csedge);
        vector<CGAL_point> outputs;
        outputs.clear();
        do {
          const CGAL_point& cpt = csedge->source()->source()->point();
          outputs.push_back(cpt);
          csedge = csedge->next();
        } while(csedge != sh_end);
        vector<point3d>  pts;
        pts.clear();
        CG_tools::convertCgalPointTo3d(outputs, pts);
        if(first){//boundary
          first = false;
          boundarys.push_back(pts);
        }else{//hole
          holesin.push_back(pts);
        }
      }
    }
    holes.push_back(holesin);
  }

  if(boundarys.empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronCGALFace(Nef_polyhedron& nefPoly, vector<vector<CGAL_point>>& boundarys, vector<vector<vector<CGAL_point>>>& holes)
{
  boundarys.clear();
  holes.clear();

  if(nefPoly.is_empty())
    return false;

  for(Nef_half_const_iterator face = nefPoly.halffacets_begin() ; face != nefPoly.halffacets_end() ; face++){
    Volume_const_handle vlo = face->incident_volume();
    if(!(vlo->mark())) continue;
    Nef_half_facet_const_iterator fci;
    bool first(true);
    vector<vector<CGAL_point>> holesin;
    holesin.clear();
    for(fci = face->facet_cycles_begin(); fci != face->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()){
        SHalfedge_const_handle csedge = SHalfedge_const_handle(fci);
        SHalfedge_const_handle sh_end(csedge);
        vector<CGAL_point> outputs;
        outputs.clear();
        do {
          const CGAL_point& cpt = csedge->source()->source()->point();
          outputs.push_back(cpt);
          csedge = csedge->next();
        } while(csedge != sh_end);
        if(first){//boundary
          first = false;
          boundarys.push_back(outputs);
        }else{//hole
          holesin.push_back(outputs);
        }
      }
    }
    holes.push_back(holesin);
  }

  if(boundarys.empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronFaces(Nef_polyhedron& nefPoly, vector<vector<point3d>>& boundarys, vector<vector<vector<point3d>>>& holes)
{
  boundarys.clear();
  holes.clear();

  if(nefPoly.is_empty())
    return false;

  Nef_half_const_iterator face;
  CGAL_forall_facets(face, nefPoly) {
    Nef_half_facet_const_iterator fci;
    bool first(true);
    vector<vector<point3d>> holesin;
    holesin.clear();
    for(fci = face->facet_cycles_begin(); fci != face->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()){
        SHalfedge_const_handle csedge = SHalfedge_const_handle(fci);
        SHalfedge_const_handle sh_end(csedge);
        vector<CGAL_point> outputs;
        outputs.clear();
        do {
          const CGAL_point& cpt = csedge->source()->source()->point();
          outputs.push_back(cpt);
          csedge = csedge->next();
        } while(csedge != sh_end);
        vector<point3d> pts;
        pts.clear();
        CG_tools::convertCgalPointTo3d(outputs, pts);
        if(first){//boundary
          first = false;
          boundarys.push_back(pts);
        }else{//hole
          holesin.push_back(pts);
        }
      }
    }
    holes.push_back(holesin);
  }

  if(boundarys.empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronFaceByVolume(Nef_polyhedron& nefPoly, vector<vector<vector<point3d>>>& boundarys, vector<vector<vector<vector<point3d>>>>& holes)
{
  boundarys.clear();
  holes.clear();

  if(nefPoly.is_empty() || nefPoly.number_of_volumes() <= 1)
    return false;

  boundarys.resize(nefPoly.number_of_volumes()-1);
  holes.resize(nefPoly.number_of_volumes()-1);

  vector<Volume_const_handle> vols;

  for(Nef_half_const_iterator face=nefPoly.halffacets_begin(); face!=nefPoly.halffacets_end(); face++) {
    uint i;
    Volume_const_handle vlo = face->incident_volume();
    if(!(vlo->mark())) continue;
    for(i=0; i<vols.size(); i++) {
      if(vlo == vols[i]) {
        //cout<<"Volume: "<<i<<endl;
        break;
      }
    }
    if(i == vols.size()) {
      vols.push_back(vlo);
      //cout<<"Volume: "<<i<<endl;
    }

    Nef_half_facet_const_iterator fci;
    bool first(true);
    vector<vector<point3d>> holesin;
    holesin.clear();
    for(fci = face->facet_cycles_begin(); fci != face->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()){
        SHalfedge_const_handle csedge = SHalfedge_const_handle(fci);
        SHalfedge_const_handle sh_end(csedge);
        vector<CGAL_point> outputs;
        outputs.clear();
        do {
          const CGAL_point& cpt = csedge->source()->source()->point();
          outputs.push_back(cpt);
          csedge = csedge->next();
        } while(csedge != sh_end);
        vector<point3d>  pts;
        pts.clear();
        CG_tools::convertCgalPointTo3d(outputs, pts);
        if(first){//boundary
          first = false;
          boundarys[i].push_back(pts);
        }else{//hole
          holesin.push_back(pts);
        }
      }
    }
    holes[i].push_back(holesin);
  }

  if(boundarys[0].empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronCGALFaceByVolume(Nef_polyhedron& nefPoly, vector<vector<vector<CGAL_point>>>& boundarys, vector<vector<vector<vector<CGAL_point>>>>& holes)
{
  boundarys.clear();
  holes.clear();

  if(nefPoly.is_empty() || nefPoly.number_of_volumes() <= 1)
    return false;

  boundarys.resize(nefPoly.number_of_volumes()-1);
  holes.resize(nefPoly.number_of_volumes()-1);

  vector<Volume_const_handle> vols;

  for(Nef_half_const_iterator face=nefPoly.halffacets_begin(); face!=nefPoly.halffacets_end(); face++) {
    uint i;
    Volume_const_handle vlo = face->incident_volume();
    if(!(vlo->mark())) continue;
    for(i=0; i<vols.size(); i++) {
      if(vlo == vols[i]) break;
    }
    if(i == vols.size()) vols.push_back(vlo);

    Nef_half_facet_const_iterator fci;
    bool first(true);
    vector<vector<CGAL_point>> holesin;
    holesin.clear();
    for(fci = face->facet_cycles_begin(); fci != face->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()){
        SHalfedge_const_handle csedge = SHalfedge_const_handle(fci);
        SHalfedge_const_handle sh_end(csedge);
        vector<CGAL_point> outputs;
        outputs.clear();
        do {
          const CGAL_point& cpt = csedge->source()->source()->point();
          outputs.push_back(cpt);
          csedge = csedge->next();
        } while(csedge != sh_end);
        if(first){//boundary
          first = false;
          boundarys[i].push_back(outputs);
        }else{//hole
          holesin.push_back(outputs);
        }
      }
    }
    holes[i].push_back(holesin);
  }

  if(boundarys[0].empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronFacesByVolume(Nef_polyhedron& nefPoly, vector<vector<vector<point3d>>>& boundarys, vector<vector<vector<vector<point3d>>>>& holes)
{
  boundarys.clear();
  holes.clear();

  if(nefPoly.is_empty() || nefPoly.number_of_volumes() <= 0)
    return false;

  boundarys.resize(nefPoly.number_of_volumes());
  holes.resize(nefPoly.number_of_volumes());

  vector<Volume_const_handle> vols;

  for(Nef_half_const_iterator face=nefPoly.halffacets_begin(); face!=nefPoly.halffacets_end(); face++) {
    uint i;
    Volume_const_handle vlo = face->incident_volume();

    for(i=0; i<vols.size(); i++)
      if(vlo == vols[i]) break;
    if(i == vols.size()) vols.push_back(vlo);
    
    Nef_half_facet_const_iterator fci;
    bool first(true);
    vector<vector<point3d>> holesin;
    holesin.clear();
    for(fci = face->facet_cycles_begin(); fci != face->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()){
        SHalfedge_const_handle csedge = SHalfedge_const_handle(fci);
        SHalfedge_const_handle sh_end(csedge);
        vector<CGAL_point> outputs;
        outputs.clear();
        do {
          const CGAL_point& cpt = csedge->source()->source()->point();
          outputs.push_back(cpt);
          csedge = csedge->next();
        } while(csedge != sh_end);
        vector<point3d> pts;
        pts.clear();
        CG_tools::convertCgalPointTo3d(outputs, pts);
        if(first){//boundary
          first = false;
          boundarys[i].push_back(pts);
        }else{//hole
          holesin.push_back(pts);
        }
      }
    }
    holes[i].push_back(holesin);
  }

  if(boundarys[0].empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronCGALFacesByVolume(Nef_polyhedron& nefPoly, vector<vector<vector<CGAL_point>>>& boundarys, vector<vector<vector<vector<CGAL_point>>>>& holes)
{
  boundarys.clear();
  holes.clear();

  if(nefPoly.is_empty() || nefPoly.number_of_volumes() <= 0)
    return false;

  boundarys.resize(nefPoly.number_of_volumes());
  holes.resize(nefPoly.number_of_volumes());

  vector<Volume_const_handle> vols;

  for(Nef_half_const_iterator face=nefPoly.halffacets_begin(); face!=nefPoly.halffacets_end(); face++) {
    uint i;
    Volume_const_handle vlo = face->incident_volume();

    for(i=0; i<vols.size(); i++)
      if(vlo == vols[i]) break;
    if(i == vols.size()) vols.push_back(vlo);

    Nef_half_facet_const_iterator fci;
    bool first(true);
    vector<vector<CGAL_point>> holesin;
    holesin.clear();
    for(fci = face->facet_cycles_begin(); fci != face->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()){
        SHalfedge_const_handle csedge = SHalfedge_const_handle(fci);
        SHalfedge_const_handle sh_end(csedge);
        vector<CGAL_point> outputs;
        outputs.clear();
        do {
          const CGAL_point& cpt = csedge->source()->source()->point();
          outputs.push_back(cpt);
          csedge = csedge->next();
        } while(csedge != sh_end);
        if(first){//boundary
          first = false;
          boundarys[i].push_back(outputs);
        }else{//hole
          holesin.push_back(outputs);
        }
      }
    }
    holes[i].push_back(holesin);
  }

  if(boundarys[0].empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronSurface(Nef_polyhedron& poly , vector<vector<point3d> >& reFaces)
{
  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor    mesh_vertex_descriptor;
  typedef boost::graph_traits<Surface_mesh>::face_descriptor    mesh_face_descriptor;
  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor mesh_halfedge_descriptor;

  reFaces.clear();

  if(!isManiFold(poly))
    return false;

  Surface_mesh mesh;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(poly , mesh);
  BOOST_FOREACH(mesh_face_descriptor fb, mesh.faces()){
    mesh_halfedge_descriptor eb = halfedge(fb , mesh);
    mesh_halfedge_descriptor ee(eb);
    vector<point3d> face;
    do{
      mesh_vertex_descriptor pd = source(eb , mesh);
      CGAL_point pt = mesh.point(pd);
      const double x = CGAL::to_double(pt.x());
      const double y = CGAL::to_double(pt.y());
      const double z = CGAL::to_double(pt.z());
      face.push_back(point3d(x , y , z));
      eb = next(eb , mesh);
    } while( eb != ee);

    reFaces.push_back(face);
  }

  if(reFaces.empty())
    return false;
  else
    return true;
}

bool CG_tools::getNefPolyhedronCGALSurface(Nef_polyhedron& poly , vector<vector<CGAL_point> >& reFaces)
{
  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor    mesh_vertex_descriptor;
  typedef boost::graph_traits<Surface_mesh>::face_descriptor    mesh_face_descriptor;
  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor mesh_halfedge_descriptor;

  reFaces.clear();

  if(!isManiFold(poly))
    return false;

  Surface_mesh mesh;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(poly , mesh);
  BOOST_FOREACH(mesh_face_descriptor fb, mesh.faces()){
    mesh_halfedge_descriptor eb = halfedge(fb , mesh);
    mesh_halfedge_descriptor ee(eb);
    vector<CGAL_point> face;
    do{
      mesh_vertex_descriptor pd = source(eb , mesh);
      CGAL_point pt = mesh.point(pd);
      face.push_back(pt);
      eb = next(eb , mesh);
    } while( eb != ee);

    reFaces.push_back(face);
  }
  
  if(reFaces.empty())
    return false;
  else
    return true;
}

bool CG_tools::isCollinear(const vector<CGAL_point>& pts , const CGAL_point& ptc)
{
  int size = pts.size();
  if(size < 2 ) return false;

  const CGAL_point& pt1 = pts[size-2];
  const CGAL_point& pt2 = pts[size-1];
  return CGAL::collinear(pt1 , pt2 , ptc);
}

bool CG_tools::isClockWise(const vector<point3d>& pts, point3d& ptc)
{
  int Ibound[6];
  double bound[6];

  if(pts.size() < 3) return false;
  
  GetSurfBound(pts, bound, Ibound);
  CGAL_point pt1, pt2, pt3, pt4;
  if(bound[1]-bound[0]>0.0) {
    uint imin = Ibound[0];
    uint inext = (imin+1)%pts.size();
    uint innext = (imin+2)%pts.size();
    pt1 = CGAL_point(pts[imin]._xx, pts[imin]._yy, pts[imin]._zz);
    pt2 = CGAL_point(pts[inext]._xx, pts[inext]._yy, pts[inext]._zz);
    pt3 = CGAL_point(pts[innext]._xx, pts[innext]._yy, pts[innext]._zz);
  }
  else {
    uint imin = Ibound[2];
    uint inext = (imin+1)%pts.size();
    uint innext = (imin+2)%pts.size();
    pt1 = CGAL_point(pts[imin]._xx, pts[imin]._yy, pts[imin]._zz);
    pt2 = CGAL_point(pts[inext]._xx, pts[inext]._yy, pts[inext]._zz);
    pt3 = CGAL_point(pts[innext]._xx, pts[innext]._yy, pts[innext]._zz);
  }  
  pt4 = CGAL_point(ptc._xx, ptc._yy, ptc._zz);

  switch(CGAL::orientation(pt1, pt2, pt3, pt4)){
    case CGAL::COPLANAR:
      return false;
      break;
    case CGAL::POSITIVE:
      return true;
      break;
    case CGAL::NEGATIVE:
      return false;
      break;
    default:
      return false;
      break;
  }
  return false;
}

bool CG_tools::isClockWiseCGAL(const vector<CGAL_point>& pts, CGAL_point& ptc)
{
  if(pts.size() < 3) return false;

  CGAL_vector vec1(ptc.x()-pts[0].x() , ptc.y()-pts[0].y() , ptc.z()-pts[0].z());
  CGAL_vector vec2(pts[1].x()-pts[0].x() , pts[1].y()-pts[0].y() , pts[1].z()-pts[0].z());
  CGAL_vector vec3(pts[2].x()-pts[1].x() , pts[2].y()-pts[1].y() , pts[2].z()-pts[1].z());
  switch(CGAL::orientation(vec2, vec3, vec1)){
    case CGAL::COPLANAR:
      return false;
      break;
    case CGAL::POSITIVE:
      return true;
      break;
    case CGAL::NEGATIVE:
      return false;
      break;
    default:
      return false;
      break;
  }
  return false;
}

bool CG_tools::isClosed(vector<vector<int>> &polyIndexs, const vector<point3d> &points)
{
  Polyhedron surfs;

  createBasePolyhedron(polyIndexs, points, surfs);

  if(!(surfs.is_empty()) && surfs.is_closed())
    return true;
  else
    return false;
}

bool CG_tools::isManiFold(Nef_polyhedron& poly)
{
  if(!poly.is_empty() && poly.is_simple())
    return true;
  else
    return false;
}

bool CG_tools::convertCgalPointTo3d(const vector<CGAL_point> & originPts , vector<point3d>& pts)
{
  for(vector<CGAL_point>::const_iterator it = originPts.begin() ; it != originPts.end(); it++){
    const CGAL_point& pt = *it;
    const double x = CGAL::to_double(pt.x());
    const double y = CGAL::to_double(pt.y());
    const double z = CGAL::to_double(pt.z());
    pts.push_back(point3d(x,y,z));
  }
  return true;
}

bool CG_tools::nefPolyhedronUnion(Nef_polyhedron& nef1 , Nef_polyhedron& nef2 , Nef_polyhedron& result)
{
  result = nef1+nef2;
  if(!(result.is_empty()))
    return true;
  else
    return false;
}

bool CG_tools::nefPolyhedronDiff(Nef_polyhedron& nef1 , Nef_polyhedron& nef2 , Nef_polyhedron& result)
{
  result = nef1-nef2;
  if(!(result.is_empty()))
    return true;
  else
    return false;
}

bool CG_tools::nefPolyhedronInter(Nef_polyhedron& nef1 , Nef_polyhedron& nef2 , Nef_polyhedron& result)
{
  result = nef1*nef2;
  if(!(result.is_empty()))
    return true;
  else
    return false;
}

void CG_tools::generateNefPolyhedronByCGALPoints(const vector<vector<CGAL_point> >& faces, Nef_polyhedron& tarPoly)
{
  // the each face should be without holes
  Generate_Nef_Polyhedron_By_Faces gen;
  gen.generateByCG(faces , tarPoly);
}

void CG_tools::generateNefPolyhedronByLocalPoints(const vector<vector<point3d> >& faces, Nef_polyhedron& tarPoly)
{
  // the each face should be without holes
  Generate_Nef_Polyhedron_By_Faces gen;
  gen.generateByLocal(faces , tarPoly);
}

Surface_mesh::Vertex_index Generate_Nef_Polyhedron_By_Faces::getVertexFromMesh(const CGAL_point& pt ,  map< CGAL_point, Surface_mesh::Vertex_index ,  CGAL_Kernel3::Less_xyz_3>& maps , Surface_mesh& mesh)
{
  Surface_mesh::Vertex_index re;
  map<CGAL_point , Surface_mesh::Vertex_index,   CGAL_Kernel3::Less_xyz_3>::iterator fend , fcur;
  fcur = maps.find(pt);
  fend = maps.end();
  if(fcur == fend){
     Surface_mesh::Vertex_index vpt = mesh.add_vertex( pt );
     maps.insert(make_pair(pt , vpt));
     re = vpt;
  } else {
     re = fcur->second;
  }
  return re;
}

void Generate_Nef_Polyhedron_By_Faces::generateByCG(const vector<vector<CGAL_point> >& faces, Nef_polyhedron& tarPoly)
{
   map< CGAL_point, Surface_mesh::Vertex_index ,  CGAL_Kernel3::Less_xyz_3> dataMaps;
   Surface_mesh mesh;

   for(uint findex = 0; findex < faces.size(); findex++){
     const vector<CGAL_point>& face = faces[findex];
     vector<Surface_mesh::Vertex_index> fves;
     for(uint index = 0 ; index < face.size(); index++){
       const CGAL_point& mpt = face[index];
       Surface_mesh::Vertex_index pt = getVertexFromMesh(mpt , dataMaps , mesh);
       fves.push_back(pt);
     }
     mesh.add_face(fves);
   }

   tarPoly = Nef_polyhedron(mesh);
}

void Generate_Nef_Polyhedron_By_Faces::generateByLocal(const vector<vector<point3d> >& faces, Nef_polyhedron& tarPoly)
{
  vector<vector<CGAL_point> >  mfaces;
  for(uint findex = 0 ; findex < faces.size(); findex++){
    const vector<point3d>& face = faces[findex];
    vector<CGAL_point> mface;
    for(uint index =0; index < face.size(); index++){
      const point3d& pt = face[index];
      CGAL_point cpt(pt._xx , pt._yy , pt._zz);
      mface.push_back(cpt);
    }
    mfaces.push_back(mface);
	}
  generateByCG(mfaces , tarPoly);
}

double DotProduct(double *v1, double *v2, int ndims)
{
  int j;
  double dp;

  dp = 0.0;
  for(j=0; j<ndims; j++)
  {
    dp += v1[j]*v2[j];
  }

  return dp;
}

void CrossProduct(double v1[3], double v2[3], double v3[3])
{
  v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

double PolyArea(vector<point3d> &vertPts)
{
  int i, j;
  double A;
  double v0[3], v1[3], v2[3], v3[3]={0};

  for(i=1; i<(int) vertPts.size()-1; i++)
  {
    v0[0] = vertPts[i]._xx-vertPts[0]._xx;
    v0[1] = vertPts[i]._yy-vertPts[0]._yy;
    v0[2] = vertPts[i]._zz-vertPts[0]._zz;
    v1[0] = vertPts[i+1]._xx-vertPts[0]._xx;
    v1[1] = vertPts[i+1]._yy-vertPts[0]._yy;
    v1[2] = vertPts[i+1]._zz-vertPts[0]._zz;
    CrossProduct(v0, v1, v2);
    for(j=0; j<3; j++)
      v3[j] += v2[j];
  }

  A = 0.5*sqrt(DotProduct(v3,v3,3));

  return A;
}

void GenTriSurfMeshIds(vector<vector<int>> &vertIds, int npts)
{
  for(uint j=0; j!=vertIds.size(); ++j)
  {
    vertIds[j].clear();
  }
  vertIds.resize(2*npts+2);

  for(int j=0; j!=npts; ++j)
  {
    vertIds[0].push_back(j);
    vertIds[1].push_back(npts+j);
  }
	
	for(int j=0; j!=npts; ++j)
  {
    int cindex = j;
    int nindex = (cindex+1)%npts;
		int hnindex = npts+(npts-nindex)%npts;
    int hcindex = npts+(npts-cindex)%npts;
    vertIds[2*j+2].push_back(cindex);
    vertIds[2*j+2].push_back(hcindex);
    vertIds[2*j+2].push_back(hnindex);
    vertIds[2*j+3].push_back(hnindex);
    vertIds[2*j+3].push_back(nindex);
    vertIds[2*j+3].push_back(cindex);
  }
}

void GenTrapSurfMeshIds(vector<vector<int>> &vertIds, int npts)
{
  for(uint j=0; j!=vertIds.size(); ++j)
  {
    vertIds[j].clear();
  }
  vertIds.resize(npts+2);

  for(int j=0; j!=npts; ++j)
  {
    vertIds[0].push_back(j);
    vertIds[1].push_back(npts+j);
  }

  for(int j=2; j!=npts+2; ++j)
  {
    int cindex = j-2;
    int nindex = (cindex+1)%npts;
    int hnindex = npts+(npts-nindex)%npts;
    int hcindex = npts+(npts-cindex)%npts;
    vertIds[j].push_back(cindex);
    vertIds[j].push_back(hcindex);
		vertIds[j].push_back(hnindex);
    vertIds[j].push_back(nindex);
  }
}

void GenBlockSurfMeshPnts(vector<point3d> &vertPts, double bound[6])
{
  point3d vertpt;

  vertPts.clear();
  vertpt = point3d(bound[0], bound[2], bound[4]);
  vertPts.push_back(vertpt);
  vertpt = point3d(bound[0], bound[3], bound[4]);
  vertPts.push_back(vertpt);
  vertpt = point3d(bound[1], bound[3], bound[4]);
  vertPts.push_back(vertpt);
  vertpt = point3d(bound[1], bound[2], bound[4]);
  vertPts.push_back(vertpt);
  vertpt = point3d(bound[0], bound[2], bound[5]);
  vertPts.push_back(vertpt);
  vertpt = point3d(bound[1], bound[2], bound[5]);
  vertPts.push_back(vertpt);
  vertpt = point3d(bound[1], bound[3], bound[5]);
  vertPts.push_back(vertpt);
  vertpt = point3d(bound[0], bound[3], bound[5]);
  vertPts.push_back(vertpt);
}

double GetBodySurfArea(vector<vector<int>> &vertIds, vector<point3d> &vertPts)
{
  double Area(0);
  vector<point3d> polysurf;

  for(uint i=0; i<vertIds.size(); i++)
  {
    polysurf.clear();
    for(uint j=0; j<vertIds[i].size(); j++)
      polysurf.push_back(vertPts[vertIds[i][j]]);
    Area += PolyArea(polysurf);
  }

  return Area;
}

void GetSurfBound(const vector<point3d> &vertPts, double bound[6], int Ibound[6])
{
  double xmin, xmax, ymin, ymax, zmin, zmax;

  double DBLMAX = std::numeric_limits<double>::max();
  xmin = DBLMAX;
  xmax = -DBLMAX;
  ymin = DBLMAX;
  ymax = -DBLMAX;
  zmin = DBLMAX;
  zmax = -DBLMAX;

  for(uint j=0; j!=vertPts.size(); ++j)
  {
    if(vertPts[j]._xx<xmin) {
      xmin = vertPts[j]._xx;
      Ibound[0] = j;
    }
    if(vertPts[j]._xx>xmax) {
      xmax=vertPts[j]._xx;
      Ibound[1] = j;
    }
    if(vertPts[j]._yy<ymin) {
      ymin=vertPts[j]._yy;
      Ibound[2] = j;
    }
    if(vertPts[j]._yy>ymax) {
      ymax=vertPts[j]._yy;
      Ibound[3] = j;
    }
    if(vertPts[j]._zz<zmin) {
      zmin=vertPts[j]._zz;
      Ibound[4] = j;
    }
    if(vertPts[j]._zz>zmax) {
      zmax=vertPts[j]._zz;
      Ibound[5] = j;
    }
  }

  bound[0] = xmin;
  bound[1] = xmax;
  bound[2] = ymin;
  bound[3] = ymax;
  bound[4] = zmin;
  bound[5] = zmax;
}

void ConvertFacesToMesh(vector<vector<point3d>> &Faces, vector<point3d> &vertPts, vector<vector<int>> &vertIds)
{
  int i, j, k;
  int vertcnt, surfcnt;
  bool compt;
  double abstol(1e-8);
  point3d vertpt;
  vector<int> vertId;

  surfcnt = Faces.size();
  if(surfcnt == 0) return;

  vertPts.clear();
  vertIds.clear();
  vertId.clear();

  vertcnt = Faces[0].size();
  for(i=0; i!=vertcnt; ++i)
  {
    vertpt = Faces[0][i];
    vertPts.push_back(vertpt);
    vertId.push_back(i);
  }
  vertIds.push_back(vertId);

  for(j=1; j!=surfcnt; ++j)
  {
    vertId.clear();
    vertcnt = Faces[j].size();
    for(i=0; i!=vertcnt; ++i)
    {
      vertpt = Faces[j][i];
      compt = false;
      for(k=vertPts.size()-1; k!=-1; --k)
      {
        if(fabs(vertpt._xx-vertPts[k]._xx)<abstol && fabs(vertpt._yy-vertPts[k]._yy)<abstol && fabs(vertpt._zz-vertPts[k]._zz)<abstol)
        {
          vertId.push_back(k);
          compt = true;
          break;
        }
      }
      if(!compt)
      {
        vertPts.push_back(vertpt);
        vertId.push_back(vertPts.size()-1);
      }
    }
    vertIds.push_back(vertId);
  }
}

void MapFacesToMesh(vector<vector<point3d>> &Faces, vector<point3d> &vertPts, vector<vector<int>> &vertIds)
{
  int i, j, k;
  int vertcnt, surfcnt;
  point3d vertpt;
  vector<int> vertId;
  map< point3d, int, Less_point3d > vertsmap;
  map< point3d, int, Less_point3d >::iterator fcur;

  surfcnt = Faces.size();
  if(surfcnt == 0) return;

  vertPts.clear();
  vertIds.clear();
  vertId.clear();

  vertcnt = Faces[0].size();
  for(i=0; i!=vertcnt; ++i)
  {
    vertpt = Faces[0][i];
    vertsmap.insert(make_pair(vertpt, vertPts.size()));
    vertPts.push_back(vertpt);
    vertId.push_back(i);
  }
  vertIds.push_back(vertId);

  for(j=1; j!=surfcnt; ++j)
  {
    vertId.clear();
    vertcnt = Faces[j].size();
    for(i=0; i!=vertcnt; ++i)
    {
      vertpt = Faces[j][i];
      fcur = vertsmap.find(vertpt);
      if(fcur != vertsmap.end()) {
        vertId.push_back(fcur->second);
      } else {
        vertsmap.insert(make_pair(vertpt, vertPts.size()));
        vertPts.push_back(vertpt);
        vertId.push_back(vertPts.size()-1);
      }
    }
    vertIds.push_back(vertId);
  }
}

void SaveSurfVec(vector<vector<point3d>> &Divsurfs, vector<vector<vector<point3d>>> &Divholes, ofstream &fout)
{
  uint i, j, k;

  fout<<"surf number: "<<Divsurfs.size()<<endl;
  for(i=0; i!=Divsurfs.size(); ++i)
  {
    fout<<Divsurfs[i].size()<<endl;
    for(j=0; j!=Divsurfs[i].size(); ++j)
      fout<<scientific<<setprecision(15)<<Divsurfs[i][j]._xx<<", "<<Divsurfs[i][j]._yy<<", "<<Divsurfs[i][j]._zz<<endl;
    if(Divholes.size()>0 && Divholes[i].size()>0)
    {
      fout<<"hole number: "<<Divholes[i].size()<<endl;
      for(j=0; j!=Divholes[i].size(); ++j)
      {
        fout<<Divholes[i][j].size()<<endl;
        for(k=0; k!=Divholes[i][j].size(); ++k)
          fout<<scientific<<setprecision(15)<<Divholes[i][j][k]._xx<<", "<<Divholes[i][j][k]._yy<<", "<<Divholes[i][j][k]._zz<<endl;
      }
    }
  }

}

void SaveVertVec(vector<point3d> &vertPts, const string &filename)
{
  uint k;
  ofstream fout(filename.c_str());

  fout<<vertPts.size()<<endl;
  for(k=0; k!=vertPts.size(); ++k)
  {
    fout<<scientific<<setprecision(15)<<vertPts[k]._xx<<", "<<vertPts[k]._yy<<", "<<vertPts[k]._zz<<endl;
  }

  fout.close();
}

void SaveMeshVec(vector<point3d> &vertPts, vector<vector<int>> &vertIds, const string &filename)
{
  uint i, j;
  ofstream fout(filename.c_str());

  fout<<vertPts.size()<<", "<<vertIds.size()<<endl;
  for(j=0; j!=vertPts.size(); ++j)
  {
    fout<<scientific<<setprecision(15)<<vertPts[j]._xx<<" "<<vertPts[j]._yy<<" "<<vertPts[j]._zz<<endl;
  }
  for(j=0; j!=vertIds.size(); ++j)
  {
    for(i=0; i!=vertIds[j].size(); ++i)
    {
      fout<<vertIds[j][i]<<" ";
    }
    fout<<endl;
  }

  fout.close();
}

void CheckSmall(vector<vector<point3d>> &surfs, vector<vector<vector<point3d>>> &holes, double tol)
{
  uint j, k, l, m, surfid, vertid, holeid(9999), overlapnum(0);
  double dis, dx, dy, dz;
  bool div(false);

  for(k=0; k!=surfs.size(); ++k)
  {
    for(l=0; l!=surfs[k].size(); ++l)
    {
      j = (l+1)%surfs[k].size();
      dx = surfs[k][j]._xx-surfs[k][l]._xx;
      dy = surfs[k][j]._yy-surfs[k][l]._yy;
      dz = surfs[k][j]._zz-surfs[k][l]._zz;
      dis = sqrt(dx*dx+dy*dy+dz*dz);
      if(dis<tol)
      {
        div = true;
        surfid = k;
        vertid = l;
        overlapnum++;
        //cout<<"surfid: "<<surfid<<", "<<"vertid: "<<vertid<<", "<<"holeid: "<<holeid<<endl;
      }
    }
    if(holes.size()>0 && holes[k].size()>0)
    {
      for(l=0; l!=holes[k].size(); ++l)
      {
        for(m=0; m!=holes[k][l].size(); ++m)
        {
          j = (m+1)%holes[k][l].size();
          dx = holes[k][l][j]._xx-holes[k][l][m]._xx;
          dy = holes[k][l][j]._yy-holes[k][l][m]._yy;
          dz = holes[k][l][j]._zz-holes[k][l][m]._zz;
          dis = sqrt(dx*dx+dy*dy+dz*dz);
          if(dis<tol)
          {
            div = true;
            surfid = k;
            holeid = l;
            vertid = m;
            overlapnum++;
            //cout<<"surfid: "<<surfid<<", "<<"vertid: "<<vertid<<", "<<"holeid: "<<holeid<<endl;
          }
        }
      }
    }
  }

  if(div)
    cout<<"There are overlapped neighbor points in a facet! Number: "<<overlapnum/2<<endl;
}

}
