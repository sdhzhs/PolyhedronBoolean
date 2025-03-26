#ifndef POLY_BASICS_H
#define POLY_BASICS_H
#include <vector>
#include <string>
#include <fstream>

namespace cgalutil {

typedef unsigned int                                                          uint;
typedef double                                                                realnum;

/*
 * @brief Define Basic Data Structure for double precision Geometric Objects
 */

struct point3d{
    point3d( double xx, double yy, double zz )
      :   _xx( xx )
        , _yy( yy )
        , _zz( zz )
    {
      ;
    }
    point3d() {;}
    double _xx;
    double _yy;
    double _zz;
};

struct Less_point3d{
    bool operator()(const point3d &p1, const point3d &p2) const
    {
      if(p1._xx < p2._xx) {
        return true;
      } else if(p1._xx > p2._xx) {
        return false;
      } else {
        if(p1._yy < p2._yy) {
          return true;
        } else if(p1._yy > p2._yy) {
          return false;
        } else {
          if(p1._zz < p2._zz) {
            return true;
          } else {
            return false;
          }
        }
      }
    }
};

struct polyprop{
    polyprop(uint von, uint fn, uint en, uint ven)
      :   volnum(von)
        , facenum(fn)
        , edgenum(en)
        , vertnum(ven)
    {
      ;
    }
    polyprop():volnum(0), facenum(0), edgenum(0), vertnum(0) {;}
    uint volnum, facenum, edgenum, vertnum;
};

struct boundbox3d{
    boundbox3d( point3d& p1, point3d& p2 )
      :   pmin( p1 )
        , pmax( p2 )
    {
      ;
    }
    boundbox3d( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax )
    {
      pmin._xx = xmin;
      pmin._yy = ymin;
      pmin._zz = zmin;
      pmax._xx = xmax;
      pmax._yy = ymax;
      pmax._zz = zmax;
    }
    boundbox3d() {;}
    point3d pmin, pmax;
};

/*
 * @brief Define Geometric Ordering and Calculation Functions for Basic Operations of Data Structure
 */

double DotProduct(double *v1, double *v2, int ndims);
void CrossProduct(double v1[3], double v2[3], double v3[3]);
double PolyArea(std::vector<point3d> &vertPts);
void GenBlockSurfMeshPnts(std::vector<point3d> &vertPts, double bound[6]);
void GenTrapSurfMeshIds(std::vector<std::vector<int>> &vertIds, int npts);
void GenTriSurfMeshIds(std::vector<std::vector<int>> &vertIds, int npts);
double GetBodySurfArea(std::vector<std::vector<int>> &vertIds, std::vector<point3d> &vertPts);
void GetSurfBound(const std::vector<point3d> &vertPts, double bound[6], int Ibound[6]);
void ConvertFacesToMesh(std::vector<std::vector<point3d>> &Faces, std::vector<point3d> &vertPts, std::vector<std::vector<int>> &vertIds);
void MapFacesToMesh(std::vector<std::vector<point3d>> &Faces, std::vector<point3d> &vertPts, std::vector<std::vector<int>> &vertIds);
void SaveSurfVec(std::vector<std::vector<point3d>> &Divsurfs, std::vector<std::vector<std::vector<point3d>>> &Divholes, std::ofstream &fout);
void SaveVertVec(std::vector<point3d> &vertPts, const std::string &filename);
void SaveMeshVec(std::vector<point3d> &vertPts, std::vector<std::vector<int>> &vertIds, const std::string &filename);
void CheckSmall(std::vector<std::vector<point3d>> &Divsurfs, std::vector<std::vector<std::vector<point3d>>> &Divholes, double tol);

}

#endif
