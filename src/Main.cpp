#include "PolyOpUtils.h"

using namespace std;
using namespace cgalutil;

int main() {
  unsigned i, j, k, l, m;
  int accuracy=1, pointnum;
  bool div;
  realnum bound[6]={0}, xx, yy, zz, theta, thick, dv1, dv2, v1[2], v2[2];
  string token;
  boundbox3d bbox;
  point3d vertpt;
  CGAL_point p;
  CGAL_vector normal;
  CGAL_plane face;
  Nef_polyhedron bodyf, bodym, bodyd;
  Polyhedron polysurf;
  Surface_mesh mesh;
  vector<vector<point3d>> surfs;
  vector<vector<vector<point3d>>> holes;
  vector<vector<vector<point3d>>> volsurfs;
  vector<vector<vector<vector<point3d>>>> volholes;
  vector<vector<vector<CGAL_point>>> cvolsurfs;
  vector<vector<vector<vector<CGAL_point>>>> cvolholes;
  vector<realnum> vertshp;
  vector<point3d> vertPts;
  vector<vector<int>> vertIds;
  vector<int> vertId;
  ifstream fin;
  ofstream fout;

//Test division method of rational numbers from different multi-precision libraries
#if 1
  CGAL::Gmpq e1(1.2);
  CGAL::Gmpq e2(0.7);
  cout<<e1<<" "<<e2<<endl;

  CGAL::Gmpq e3 = e1/e2;
  cout<<e3<<endl;
  cout<<CGAL::to_double(e3)<<endl;

  CGAL::Quotient<CGAL::MP_Float> e4(1.2);
  CGAL::Quotient<CGAL::MP_Float> e5(0.7);
  cout<<e4<<" "<<e5<<endl;

  CGAL::Quotient<CGAL::MP_Float> e6 = e4/e5;
  cout<<e6<<endl;
  cout<<CGAL::to_double(e6)<<endl;

  CORE::BigFloat e7(1.2);
  CORE::BigFloat e8(0.7);
  cout<<e7<<" "<<e8<<endl;

  CORE::BigFloat e9 = sqrt(e7);
  cout<<e9<<endl;
  cout<<CGAL::to_double(e9)<<endl;

  CORE::BigFloat e10 = e7/e8;
  cout<<e10<<endl;
  cout<<CGAL::to_double(e10)<<endl;

  cout<<"--------------------------"<<endl;
#endif

//Test OFF_to_nef_3 interface in Nef_Polyhedron CGAL Package
#if 1
  fin.open("cgal_sufaces.off");

  size_t discard = CGAL::OFF_to_nef_3(fin, bodym);
  fin.close();
  cout<<"discarded: "<<discard<<endl;

  fin.open("cgal_sufaces.off");
  fin>>mesh;
  fin.close();
  bodyf = Nef_polyhedron(mesh);

  if(bodym == bodyf)
    cout<<"Same!"<<endl;
  else
    cout<<"Diff!"<<endl;

  cout<<"--------------------------"<<endl;
#endif

//Test a complete boolean operation process of two simple cubes, including construction of two original Nef_Polyhedra, difference boolean operation and extraction
//of B-rep faces in result Nef_Polyhedron into faces lists with holes lists, both divided by volumes. Boolean operation between a Nef_Polyhedron and a plane
//is also tested.
#if 1
  bound[0] = 2.0;
  bound[1] = 4.0;
  bound[2] = 2.0;
  bound[3] = 4.0;
  bound[4] = 1.343;
  bound[5] = 2.0;
  GenBlockSurfMeshPnts(vertPts, bound);
  GenTrapSurfMeshIds(vertIds, 4);

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodym);

  div = CG_tools::getNefPolyhedronSurface(bodym, surfs);
  if(!div)
    cerr<<"Divide no surface!"<<endl;

  //ConvertFacesToMesh(surfs, vertPts, vertIds);
  MapFacesToMesh(surfs, vertPts, vertIds);
  SaveMeshVec(vertPts, vertIds, "cgal_mesh.out");

  bound[0] = 2.5;
  bound[1] = 3.5;
  bound[2] = 2.5;
  bound[3] = 3.5;
  bound[4] = 0.0;
  bound[5] = 1.343;

  GenBlockSurfMeshPnts(vertPts, bound);
  GenTrapSurfMeshIds(vertIds, 4);
  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);
	
  div = CG_tools::nefPolyhedronDiff(bodym, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not diff!"<<endl;

  if(bodyf.is_empty()) cout<<"nef is empty!"<<endl;
  if(!bodyf.is_simple()) cout<<"nef is not simple!"<<endl;
  if(!bodyf.is_valid()) cout<<"nef is not valid!"<<endl;

  mesh.clear();
  CGAL::convert_nef_polyhedron_to_polygon_mesh(bodyf, mesh);
  if(!CGAL::is_closed(mesh)) cout<<"mesh is not closed"<<endl;
  if(!mesh.is_valid()) cout<<"mesh is not valid"<<endl;

  fout.open("cgal_mesh.off");
  fout<<mesh;
  fout.close();

  div = CG_tools::getNefPolyhedronFaceByVolume(bodyf, volsurfs, volholes);
  if(!div)
    cerr<<"Divide no surface!"<<endl;

  fout.open("cgal_surf.out");

  fout<<"volume number: "<<volsurfs.size()<<endl;
  for(j=0; j!=volsurfs.size(); ++j)
  {
    SaveSurfVec(volsurfs[j], volholes[j], fout);
  }

  fout.close();

  CG_tools::getPolyhedronBoundBox(bodyf, bbox);
  cout<<"xmin: "<<bbox.pmin._xx<<" xmax: "<<bbox.pmax._xx<<" ymin: "<<bbox.pmin._yy<<" ymax: "<<bbox.pmax._yy<<" zmin: "<<bbox.pmin._zz<<" zmax: "<<bbox.pmax._zz<<endl;

  p= CGAL_point(3, 3, 0);
  normal = CGAL_vector(1, 0, 0);
  face = CGAL_plane(p, normal);
  bodyf = bodym.intersection(face, Nef_polyhedron::OPEN_HALFSPACE);

  if(bodyf.is_empty()) cout<<"nef is empty!"<<endl;
  if(!bodyf.is_simple()) cout<<"nef is not simple!"<<endl;
  if(!bodyf.is_valid()) cout<<"nef is not valid!"<<endl;

  div = CG_tools::getNefPolyhedronFaceByVolume(bodyf, volsurfs, volholes);
  if(!div)
    cerr<<"Divide no surface!"<<endl;

  fout.open("cgal_plane.out");

  fout<<"volume number: "<<volsurfs.size()<<endl;
  for(j=0; j!=volsurfs.size(); ++j)
  {
    SaveSurfVec(volsurfs[j], volholes[j], fout);
  }

  fout.close();

  cout<<"--------------------------"<<endl;
#endif

//Test union boolean operations between Nef_Polyhedra which are not manifold (single or multiple faces which are not closed), the result Nef_Polyhedron can include a
//closed volume (polyhedron in geometric meaning) and the faces in this closed volume can be extracted.
#if 1
  vertPts.clear();
  vertIds.clear();

  vertPts.push_back(point3d(0.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);

  vertId.clear();
  vertId.push_back(0);
  vertPts.push_back(point3d(0.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertId.push_back(1);
  vertIds.push_back(vertId);

  vertId.clear();
  vertId.push_back(3);
  vertId.push_back(2);
  vertPts.push_back(point3d(2.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);
	
  CG_tools::createNefPolyhedron(vertIds, vertPts, bodym);

  cout<<"number of volume: "<<bodym.number_of_volumes()<<endl;
	
  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(0.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);

  div = CG_tools::nefPolyhedronUnion(bodym, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(0.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);
	
  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);

  div = CG_tools::nefPolyhedronUnion(bodyf, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;
	
  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(0.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);
	
  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);

  div = CG_tools::nefPolyhedronUnion(bodyf, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  div = CG_tools::getNefPolyhedronFaces(bodyf, surfs, holes);
  if(!div)
    cerr<<"Divide no surface!"<<endl;

  fout.open("cgal_surf_union.out");
  SaveSurfVec(surfs, holes, fout);
  fout.close();

  CG_tools::getPolyhedronBoundBox(bodyf, bbox);
  cout<<"xmin: "<<bbox.pmin._xx<<" xmax: "<<bbox.pmax._xx<<" ymin: "<<bbox.pmin._yy<<" ymax: "<<bbox.pmax._yy<<" zmin: "<<bbox.pmin._zz<<" zmax: "<<bbox.pmax._zz<<endl;
	
  cout<<"--------------------------"<<endl;
#endif

//Test union boolean operations between Nef_Polyhedra which are not manifold (single or multiple faces which are not closed), this is a degenerated case compared to
//the above example, the result Nef_Polyhedron can include more than one closed volumes and the faces can be extracted by volumes.
#if 1
  vertPts.clear();
  vertIds.clear();
  vertId.clear();

  vertPts.push_back(point3d(0.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);

  vertId.clear();
  vertId.push_back(0);
  vertPts.push_back(point3d(0.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertId.push_back(1);
  vertIds.push_back(vertId);

  vertId.clear();
  vertId.push_back(3);
  vertId.push_back(2);
  vertPts.push_back(point3d(2.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);
	
  CG_tools::createNefPolyhedron(vertIds, vertPts, bodym);

  cout<<"number of volume: "<<bodym.number_of_volumes()<<endl;
	
  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(0.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(1.5,0.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(1.5,2.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);

  div = CG_tools::nefPolyhedronUnion(bodym, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(2.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.5,2.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.5,0.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);
	
  div = CG_tools::nefPolyhedronUnion(bodyf, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;
	
  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(0.5,0.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(1.5,0.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(1.5,2.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.5,2.0,1.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);

  div = CG_tools::nefPolyhedronUnion(bodyf, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(0.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,0.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);
	
  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);
	
  div = CG_tools::nefPolyhedronUnion(bodyf, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  vertPts.clear();
  vertIds.clear();

  vertId.clear();
  vertPts.push_back(point3d(0.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(0.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,2.0));
  vertId.push_back(vertPts.size()-1);
  vertPts.push_back(point3d(2.0,2.0,0.0));
  vertId.push_back(vertPts.size()-1);
  vertIds.push_back(vertId);

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);

  div = CG_tools::nefPolyhedronUnion(bodyf, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not Inter!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  if(bodyf.interior()==bodyf) cout<<"interior success"<<endl;

  div = CG_tools::getNefPolyhedronVertex(bodyf, vertPts);
  if(!div)
    cerr<<"Divide no vertex!"<<endl;

  SaveVertVec(vertPts, "cgal_vert_union_taper.out");

  div = CG_tools::getNefPolyhedronFacesByVolume(bodyf, volsurfs, volholes);
  if(!div)
    cerr<<"Divide no surface!"<<endl; 

  fout.open("cgal_surf_union_taper.out");
  
  fout<<"volume number: "<<volsurfs.size()<<endl;
  for(j=0; j!=volsurfs.size(); ++j)
  {
    SaveSurfVec(volsurfs[j], volholes[j], fout);
  }

  fout.close();

  CheckSmall(surfs, holes, 1e-3);

  CG_tools::getPolyhedronBoundBox(bodyf, bbox);
  cout<<"xmin: "<<bbox.pmin._xx<<" xmax: "<<bbox.pmax._xx<<" ymin: "<<bbox.pmin._yy<<" ymax: "<<bbox.pmax._yy<<" zmin: "<<bbox.pmin._zz<<" zmax: "<<bbox.pmax._zz<<endl;
	
  cout<<"--------------------------"<<endl;
#endif

//Test the consistency between self-defined surfaces file format and built-in Nef_Polyhedron file format, some mesh conversion functions are also tested
#if 1
  unsigned surfcnt, vertcnt;
	
  fin.open("cgal_surf_cube.out");

  fin>>token>>token>>surfcnt;
  surfs.resize(surfcnt);
  for(k=0; k!=surfs.size(); ++k)
  {
    fin>>vertcnt;
    surfs[k].resize(vertcnt);
    for(l=0; l!=surfs[k].size(); ++l)
    {
      fin>>surfs[k][l]._xx>>token>>surfs[k][l]._yy>>token>>surfs[k][l]._zz;
    }
  }

  fin.close();
	
  ConvertFacesToMesh(surfs, vertPts, vertIds);

  SaveMeshVec(vertPts, vertIds, "cgal_mesh_cube.out");

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyf);

  fin.open("cgal_nef_cube.out");
  fin>>bodym;
  fin.close();

  if(bodym == bodyf)
    cout<<"Same!"<<endl;
  else
    cout<<"Diff!"<<endl;

  div = CG_tools::getNefPolyhedronFace(bodyf, surfs, holes);
  if(!div)
    cerr<<"Divide no surface!"<<endl;

  fout.open("cgal_surf_div.out");
  SaveSurfVec(surfs, holes, fout);
  fout.close();

  cout<<"--------------------------"<<endl;
#endif

//Test the possibility of generating a 3D complex from 2D layout polygons
#if 1
  fin.open("layout_4");
  fin>>pointnum;

  vertshp.clear();
  vertPts.clear();
  for(k=0; k!=(unsigned) pointnum; ++k)
  {
    fin>>xx>>token>>yy>>token>>zz;
    vertshp.push_back(xx);
    vertshp.push_back(yy);
    vertshp.push_back(zz);
    vertpt._xx = double(xx);
    vertpt._yy = double(yy);
    vertpt._zz = double(zz);
    vertPts.push_back(vertpt);
  }

  fin.close();
	
  reverse(vertPts.begin(), vertPts.end());
  vertshp.clear();
  for(k=0; k!=(unsigned) pointnum; ++k)
  {
    vertshp.push_back(vertPts[k]._xx);
    vertshp.push_back(vertPts[k]._yy);
    vertshp.push_back(vertPts[k]._zz);
  }

  realnum extd, cosalpha, tanalpha, orient;
	
  theta = 6*M_PI/18;
  thick = 0.2;

  GenTrapSurfMeshIds(vertIds, pointnum);
  //GenTriSurfMeshIds(vertIds, pointnum);

  extd = thick/tan(theta);
	
  for(k=0; k!=(unsigned) pointnum; ++k)
  {
    int kprev, knext;
    kprev = (pointnum+k-1)%pointnum;
    knext = (k+1)%pointnum;
    v1[0] = vertshp[3*knext]-vertshp[3*k];
    v1[1] = vertshp[3*knext+1]-vertshp[3*k+1];
    v2[0] = vertshp[3*kprev]-vertshp[3*k];
    v2[1] = vertshp[3*kprev+1]-vertshp[3*k+1];
    dv1 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
    dv2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
    v1[0] /= dv1;
    v1[1] /= dv1;
    v2[0] /= dv2;
    v2[1] /= dv2;
    cosalpha = v1[0]*v2[0]+v1[1]*v2[1];
    tanalpha = sqrt((1-cosalpha)/(1+cosalpha));
    orient = v1[0]*v2[1]-v1[1]*v2[0];
    v2[0] = v1[1];
    v2[1] = -v1[0];
    if(orient>0) {
      v1[0] = -v1[0];
      v1[1] = -v1[1];
    }
    dv2 = extd;
    dv1 = dv2/tanalpha;
    xx = vertshp[3*k]+dv1*v1[0]+dv2*v2[0];
    yy = vertshp[3*k+1]+dv1*v1[1]+dv2*v2[1];
    zz = vertshp[3*k+2]+thick;
    //zz = vertshp[3*k+2];
    vertshp.push_back(xx);
    vertshp.push_back(yy);
    vertshp.push_back(zz);
    vertpt._xx = double(xx);
    vertpt._yy = double(yy);
    vertpt._zz = double(zz);
    vertPts.push_back(vertpt);
  }

  reverse(vertPts.begin()+pointnum+1, vertPts.end());

  SaveMeshVec(vertPts, vertIds, "cgal_mesh_layout.out");

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodyd);
	
  cout<<"number of volume: "<<bodyd.number_of_volumes()<<endl;

  bound[0] = -5.0;
  bound[1] = 5.0;
  bound[2] = -5.0;
  bound[3] = 5.0;
  bound[4] = 0.0;
  bound[5] = 0.3;
  GenBlockSurfMeshPnts(vertPts, bound);

  CG_tools::createNefPolyhedron(vertIds, vertPts, bodym);
	
  cout<<"number of volume: "<<bodym.number_of_volumes()<<endl;

  div = CG_tools::nefPolyhedronDiff(bodym, bodyd, bodyf);
  if(!div)
    cerr<<"Polyhedrons not diff!"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  if(bodyf.is_empty()) cout<<"nef is empty!"<<endl;
  if(!bodyf.is_simple()) cout<<"nef is not simple!"<<endl;
  if(!bodyf.is_valid()) cout<<"nef is not valid!"<<endl;
	
  div = CG_tools::getNefPolyhedronFaceByVolume(bodyf, volsurfs, volholes);
  if(!div)
    cerr<<"Divide no surface!"<<endl;

  fout.open("cgal_surf_layout.out");

  fout<<"volume number: "<<volsurfs.size()<<endl;
  for(j=0; j!=volsurfs.size(); ++j)
  {
    SaveSurfVec(volsurfs[j], volholes[j], fout);
  }

  fout.close();

  mesh.clear();
  CGAL::convert_nef_polyhedron_to_polygon_mesh(bodyf, mesh);
  //bodyf.convert_to_polyhedron(polysurf);

  fout.open("cgal_mesh_layout.off");
  fout<<mesh;
  fout.close();

  fout.open("cgal_nef_layout.out");
  fout<<bodyf;
  fout.close();

  cout<<"--------------------------"<<endl;
#endif

//Test intersection boolean operations between Nef_Polyhedra which are planes which represent infinite half spaces, this is another way to construct 3D closed
//polyhedra and the vertices are constructed by intersections among planes, the result Nef_Polyhedron can include more than one closed volumes and the faces 
//in the closed volumes can be extracted correctly. Note that the extended geometric kernel in CGAL should be used for this test example.
#if 1
  fin.open("layout_4");
  fin>>pointnum;

  vertshp.clear();
  vertPts.clear();
  for(k=0; k!=(unsigned) pointnum; ++k)
  {
    fin>>xx>>token>>yy>>token>>zz;
    vertshp.push_back(xx);
    vertshp.push_back(yy);
    vertshp.push_back(zz);
    vertpt._xx = double(xx);
    vertpt._yy = double(yy);
    vertpt._zz = double(zz);
    vertPts.push_back(vertpt);
  }

  fin.close();
	
  reverse(vertPts.begin(), vertPts.end());
  vertshp.clear();
  for(k=0; k!=(unsigned) pointnum; ++k)
  {
    vertshp.push_back(vertPts[k]._xx);
    vertshp.push_back(vertPts[k]._yy);
    vertshp.push_back(vertPts[k]._zz);
  }

  theta = 3*M_PI/18;
  thick = 0.2;

  double px, py, pz;
  accuracy = 10;
  px = vertshp[0]*accuracy;
  py = vertshp[1]*accuracy;
  pz = vertshp[2]*accuracy;
  p = CGAL_point(px, py, pz, accuracy);
  normal = CGAL_vector(0, 0, -1);
  face = CGAL_plane(p, normal);

  bodyf = Nef_polyhedron(face);

  cout<<"success: b"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  px = vertshp[0]*accuracy;
  py = vertshp[1]*accuracy;
  pz = (vertshp[2]+thick)*accuracy;
  p = CGAL_point(px, py, pz, accuracy);
  normal = CGAL_vector(0, 0, 1);
  face = CGAL_plane(p, normal);

  bodym = Nef_polyhedron(face);
  bodyf *= bodym;

  cout<<"success: t"<<endl;

  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  for(k=0; k!=(unsigned) pointnum; ++k)
  {
    int knext;
    knext = (k+1)%pointnum;
    v1[0] = vertshp[3*knext]-vertshp[3*k];
    v1[1] = vertshp[3*knext+1]-vertshp[3*k+1];
    v2[0] = -v1[1];
    v2[1] = v1[0];
    dv2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
    zz = dv2/tan(theta);
    px = vertshp[3*k]*accuracy;
    py = vertshp[3*k+1]*accuracy;
    pz = vertshp[3*k+2]*accuracy;
    p = CGAL_point(px, py, pz, accuracy);
    px = v2[0]*accuracy;
    py = v2[1]*accuracy;
    pz = zz*accuracy;
    normal = CGAL_vector(px, py, pz, accuracy);
    face = CGAL_plane(p, normal);
    bodym = Nef_polyhedron(face);
    bodyf *= bodym;
    cout<<"success: "<<k<<endl;
  }
	
  cout<<"number of volume: "<<bodyf.number_of_volumes()<<endl;

  fout.open("cgal_nef_interplane.out");
  fout<<bodyf;
  fout.close();

  div = CG_tools::getNefPolyhedronFace(bodyf, surfs, holes);
  if(!div)
    cerr<<"Divide no surface!"<<endl;

  fout.open("cgal_surf_interplane.out");
  SaveSurfVec(surfs, holes, fout);
  fout.close();

  mesh.clear();
  CGAL::convert_nef_polyhedron_to_polygon_mesh(bodyf, mesh);

  fout.open("cgal_mesh_interplane.off");
  fout<<mesh;
  fout.close();

  cout<<"--------------------------"<<endl;
#endif

  cout<<"Task Complete!"<<endl;

  return 0;
}	
