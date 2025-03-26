#include <limits>
#include "PolyOpUtils.h"
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/OFF_to_nef_3.h>

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
  for(mesh_face_descriptor& fb : mesh.faces()){
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
  for(mesh_face_descriptor& fb : mesh.faces()){
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

bool CG_tools::convertNefToMesh(Nef_polyhedron& nef, Surface_mesh& mesh)
{
  CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, mesh);

  return true;
}

bool CG_tools::convertOFFToNef(std::ifstream& fin, Nef_polyhedron& nef)
{
  size_t discard(0);

#ifndef Extended_Kernel
  discard = CGAL::OFF_to_nef_3(fin, nef);
#endif

  if(discard != 0)
    return false;
  else
    return true;
}

bool CG_tools::LoadNefPolyhedron(std::ifstream& fin, Nef_polyhedron& nef)
{
  fin >> nef;

  return true;
}

bool CG_tools::SaveNefPolyhedron(Nef_polyhedron& nef, std::ofstream& fout)
{
  fout << nef;

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

}
