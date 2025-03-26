#include <cmath>
#include <map>
#include <iostream>
#include <iomanip>
#include <limits>
#include "PolyBasics.h"

using namespace std;

namespace cgalutil {

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
