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

/**
 * @brief 3D double precision points with constructor
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

/**
 * @brief Comparative functor between two 3D points
 */
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

/**
 * @brief Geometric properties of a polyhedron
 */
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

/**
 * @brief Boundary box of a 3D geometric object
 */
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

/**
 * @brief Dot product of two n-dimentional vectors
 */
double DotProduct(double *v1, double *v2, int ndims);
/**
 * @brief Cross product of two 3-dimentional vectors
 */
void CrossProduct(double v1[3], double v2[3], double v3[3]);
/**
 * @brief Calculate area of a 3D polygon
 */
double PolyArea(std::vector<point3d> &vertPts);
/**
 * @brief Generate all nodes in a block/cube according to boundary box
 */
void GenBlockSurfMeshPnts(std::vector<point3d> &vertPts, double bound[6]);
/**
 * @brief Generate surface mesh topology using ID for a prism with bottom face nodes number
 */
void GenTrapSurfMeshIds(std::vector<std::vector<int>> &vertIds, int npts);
/**
 * @brief Generate surface mesh topology using ID for a prism with bottom face nodes number,
 * the side faces are triangulated
 */
void GenTriSurfMeshIds(std::vector<std::vector<int>> &vertIds, int npts);
/**
 * @brief Calculate surface area of a polyhedron described by surface mesh
 */
double GetBodySurfArea(std::vector<std::vector<int>> &vertIds, std::vector<point3d> &vertPts);
/**
 * @brief Get the boundary box and index of a 3D polygon
 */
void GetSurfBound(const std::vector<point3d> &vertPts, double bound[6], int Ibound[6]);
/**
 * @brief Convert faces list to surface mesh
 */
void ConvertFacesToMesh(std::vector<std::vector<point3d>> &Faces, std::vector<point3d> &vertPts, std::vector<std::vector<int>> &vertIds);
/**
 * @brief Convert faces list to surface mesh, STL map implementation
 */
void MapFacesToMesh(std::vector<std::vector<point3d>> &Faces, std::vector<point3d> &vertPts, std::vector<std::vector<int>> &vertIds);
/**
 * @brief Save all faces and holes list in a polyhedron to a output stream
 */
void SaveSurfVec(std::vector<std::vector<point3d>> &Divsurfs, std::vector<std::vector<std::vector<point3d>>> &Divholes, std::ofstream &fout);
/**
 * @brief Save all nodes list in a polyhedron to a file
 */
void SaveVertVec(std::vector<point3d> &vertPts, const std::string &filename);
/**
 * @brief Save surface mesh in a polyhedron to a file
 */
void SaveMeshVec(std::vector<point3d> &vertPts, std::vector<std::vector<int>> &vertIds, const std::string &filename);
/**
 * @brief Check and report overlapped neighbor points in a polyhedron
 */
void CheckSmall(std::vector<std::vector<point3d>> &Divsurfs, std::vector<std::vector<std::vector<point3d>>> &Divholes, double tol);

}

#endif
