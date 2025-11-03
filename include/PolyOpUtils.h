#ifndef POLY_OP_UTILS_H
#define POLY_OP_UTILS_H
#include <vector>
#include <map>
#include "PolyBasics.h"

#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/Lazy_exact_nt.h>
//#include <CGAL/Cartesian.h>
//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Extended_homogeneous.h>
//#include <CGAL/Extended_cartesian.h>
//#include <CGAL/Filtered_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/utils.h>
//#include <CGAL/utility.h>
//#include <CGAL/minkowski_sum_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
//#include <CGAL/Nef_3/SNC_indexed_items.h>
//#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>

/*#ifndef CGAL_DONT_USE_LAZY_KERNEL
#include <CGAL/Lazy_kernel.h>
#endif*/

namespace cgalutil {
/**
 * @brief Define NumberType: FNT for Field, RNT for Ring
 */

//typedef double                                                                FNT;
//typedef CGAL::Exact_rational                                                  FNT;
//typedef CGAL::Quotient<CGAL::MP_Float>                                        FNT;
//typedef CGAL::Gmpq                                                            FNT;
//typedef CORE::BigFloat                                                        FNT;
typedef CGAL::Exact_integer                                                   RNT;

/**
 * @brief Define CGAL Kernel: Homogeneous, Cartesian and Extended Kernel
 */

//typedef CGAL::Simple_cartesian<FNT>                                           CGAL_Kernel3;
//typedef CGAL::Cartesian<FNT>                                                  CGAL_Kernel3;
//typedef CGAL::Extended_cartesian<FNT>                                         CGAL_Kernel3;
#ifdef Extended_Kernel
typedef CGAL::Extended_homogeneous<RNT>                                       CGAL_Kernel3;
#else
typedef CGAL::Exact_predicates_exact_constructions_kernel                     CGAL_Kernel3;
#endif
//typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<FNT>>                    CGAL_Kernel3;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel                   CGAL_Kernel3;
//typedef CGAL::Simple_cartesian<CGAL::Lazy_exact_nt<FNT> >                     CGAL_Kernel3;
//typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<FNT> >                       CGAL_Kernel3;

/**
 * @brief Define Basic NumberType and Geometric Objects in CGAL Kernel
 */

typedef CGAL_Kernel3::FT                                                      FT;
typedef CGAL_Kernel3::RT                                                      RT;
typedef CGAL_Kernel3::Point_3                                                 CGAL_point;
typedef CGAL_Kernel3::Vector_3                                                CGAL_vector;
typedef CGAL_Kernel3::Plane_3                                                 CGAL_plane;
typedef CGAL_Kernel3::Aff_transformation_3                                    CGAL_Afftransform_3;

/**
 * @brief Define Geometric Objects in Nef_Polyhedron CGAL Package
 */

typedef CGAL::Polyhedron_3<CGAL_Kernel3>                                      Polyhedron;          
typedef CGAL::Nef_polyhedron_3<CGAL_Kernel3>                                  Nef_polyhedron;
typedef Nef_polyhedron::Vertex_const_iterator                                 Nef_vertex_const_iterator;
typedef Nef_polyhedron::Halffacet_const_iterator                              Nef_half_const_iterator;
typedef Nef_polyhedron::Volume_const_handle                                   Volume_const_handle;
typedef Nef_polyhedron::SHalfedge_const_handle                                SHalfedge_const_handle;
typedef Nef_polyhedron::SHalfedge_around_facet_const_circulator               Nef_half_edge_around_facet_const_circulator;
typedef Nef_polyhedron::Halffacet_cycle_const_iterator                        Nef_half_facet_const_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator                          Face_edge_iterator;
typedef typename Polyhedron::Facet_iterator                                   FacetIterator;
typedef Polyhedron::Facet_handle                                              Facet_handle;
typedef CGAL::Surface_mesh<CGAL_point>                                        Surface_mesh;

/**
 * @brief Define Class to Construct Polyhedron from Surface Mesh based on Nef_Polyhedron CGAL Package
 */

template <typename Polyhedron>
class CGAL_Build_PolySet : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
{
    typedef typename Polyhedron::HalfedgeDS HDS;
    typedef CGAL::Polyhedron_incremental_builder_3<typename Polyhedron::HalfedgeDS> CGAL_Polybuilder;
  public:
    const std::vector<std::vector<int> > &_polyIndexs;
    const std::vector<point3d> &_points;
    CGAL_Build_PolySet(const std::vector<std::vector<int> > &polyIndexs , const std::vector<point3d> &points)
                      : _polyIndexs(polyIndexs)
                      , _points(points)
    {;}

    void operator()(HDS& hds) {

      CGAL_Polybuilder B(hds, true);
      B.begin_surface(_points.size(), _polyIndexs.size());

      for(unsigned i  =0 ; i < _points.size() ; i++) {
        const point3d & pt = _points[i];
        CGAL_point p( (FT)pt._xx , (FT)pt._yy , (FT)pt._zz);
        B.add_vertex(p);
      }
      for(unsigned i =0; i < _polyIndexs.size(); i++){
          std::vector<int>   pindices = _polyIndexs[i];
          if (pindices.size() >=3) {
            B.add_facet(pindices.begin(), pindices.end());
          }
       }
       B.end_surface();
    }
};

class Generate_Nef_Polyhedron_By_Faces
{
  public:
    Generate_Nef_Polyhedron_By_Faces() {;}
    ~Generate_Nef_Polyhedron_By_Faces() {;}
    void generateByCG(const std::vector<std::vector<CGAL_point> >& faces, Nef_polyhedron& tarPoly);
    void generateByLocal(const std::vector<std::vector<point3d> >& faces, Nef_polyhedron& tarPoly);
  private:
    Surface_mesh::Vertex_index getVertexFromMesh(const CGAL_point& pt ,  std::map< CGAL_point, Surface_mesh::Vertex_index ,  CGAL_Kernel3::Less_xyz_3>& maps , Surface_mesh& mesh);
};

/*
 * @brief Define static Functions to implement Construction of original Nef_Polyhedra, 3D Boolean Operations on Nef_Polyhedra and
 * Extraction of B-rep Faces in result Nef_Polyhedron to Basic Data Structure
 */

class CG_tools{
  public:
    CG_tools() {;}
    ~CG_tools() {;}

    /**
     * @brief Create Polyhedron from Surface Mesh described by Basic Data Structure
     */
    static void createBasePolyhedron(std::vector<std::vector<int> > &polyIndexs , const std::vector<point3d> &points , Polyhedron& p);

    /**
     * @brief Create Nef_Polyhedron from Surface Mesh described by Basic Data Structure
     */
    static void createNefPolyhedron(std::vector<std::vector<int> > &polyIndexs , const std::vector<point3d> &points , Nef_polyhedron& poly);

    /**
     * @brief Create Nef_Polyhedron from Faces List described by CGAL_point Lists
     */
    static void generateNefPolyhedronByCGALPoints(const std::vector<std::vector<CGAL_point> >& faces, Nef_polyhedron& tarPoly);

    /**
     * @brief Create Nef_Polyhedron from Faces List described by Basic Point Lists
     */
    static void generateNefPolyhedronByLocalPoints(const std::vector<std::vector<point3d> >& faces, Nef_polyhedron& tarPoly);

    /**
     * @brief Get Topological Properties of one Nef_Polyhedron
     */
    static bool getNefPolyhedronProp(Nef_polyhedron& poly , polyprop& nefprop);

    /**
     * @brief Get Boundary Box of one Nef_Polyhedron
     */
    static bool getPolyhedronBoundBox(Nef_polyhedron& poly, boundbox3d& boundcoord);

    /**
     * @brief Extract all Vertices in one Nef_Polyhedron into one Basic Point List
     */
    static bool getNefPolyhedronVertex(Nef_polyhedron& poly , std::vector<point3d>& vertices);

    /**
     * @brief Convert all Faces in one Nef_Polyhedron into Surface Mesh (Triangularized for Faces with Holes) and Extract all Cells into Faces List described by Basic Point
     */
    static bool getNefPolyhedronSurface(Nef_polyhedron& poly , std::vector<std::vector<point3d> >& faces);

    /**
     * @brief Convert all Faces in one Nef_Polyhedron into Surface Mesh (Triangularized for Faces with Holes) and Extract all Cells into Faces List described by CGAL_point
     */
    static bool getNefPolyhedronCGALSurface(Nef_polyhedron& poly , std::vector<std::vector<CGAL_point> >& faces);

    /**
     * @brief Extract all Marked (Constitute closed Volumes) Faces (with Holes) in one Nef_Polyhedron into one Boundary Faces List and
     * one Hole Faces List described by Basic Point, there are more than one hole faces corresponding to one boundary face
     */
    static bool getNefPolyhedronFace(Nef_polyhedron& poly , std::vector<std::vector<point3d>>& boundarys,  std::vector<std::vector<std::vector<point3d>>>& holes);

    /**
     * @brief Extract all Marked (Constitute closed Volumes) Faces (with Holes) in one Nef_Polyhedron into one Boundary Faces List and
     * one Hole Faces List described by CGAL_point, there are more than one hole faces corresponding to one boundary face
     */
    static bool getNefPolyhedronCGALFace(Nef_polyhedron& poly , std::vector<std::vector<CGAL_point>>& boundarys,  std::vector<std::vector<std::vector<CGAL_point>>>& holes);

    /**
     * @brief Extract all Half Faces (with Holes) in one Nef_Polyhedron into one Boundary Faces List and one Hole Faces List described by Basic Point
     */
    static bool getNefPolyhedronFaces(Nef_polyhedron& poly , std::vector<std::vector<point3d>>& boundarys,  std::vector<std::vector<std::vector<point3d>>>& holes);

    /**
     * @brief Extract all Marked (Constitute closed Volumes) Faces (with Holes) in one Nef_Polyhedron into one Boundary Faces List and one Hole Faces List described by Basic Point,
     * the top Level of these two Lists are Divided by closed Volumes in this Nef_Polyhedron
     */
    static bool getNefPolyhedronFaceByVolume(Nef_polyhedron& poly , std::vector<std::vector<std::vector<point3d>>>& boundarys,  std::vector<std::vector<std::vector<std::vector<point3d>>>>& holes);

    /**
     * @brief Extract all Marked (Constitute closed Volumes) Faces (with Holes) in one Nef_Polyhedron into one Boundary Faces List and one Hole Faces List described by CGAL_point,
     * the top Level of these two Lists are Divided by closed Volumes in this Nef_Polyhedron
     */
    static bool getNefPolyhedronCGALFaceByVolume(Nef_polyhedron& poly , std::vector<std::vector<std::vector<CGAL_point>>>& boundarys,  std::vector<std::vector<std::vector<std::vector<CGAL_point>>>>& holes);

    /**
     * @brief Extract all Half Faces (with Holes) in one Nef_Polyhedron into one Boundary Faces List and one Hole Faces List described by Basic Point,
     * the top Level of these two Lists are Divided by all Volumes (including Outter Space) in this Nef_Polyhedron
     */
    static bool getNefPolyhedronFacesByVolume(Nef_polyhedron& poly, std::vector<std::vector<std::vector<point3d>>>& boundarys, std::vector<std::vector<std::vector<std::vector<point3d>>>>& holes);

    /**
     * @brief Extract all Half Faces (with Holes) in one Nef_Polyhedron into one Boundary Faces List and one Hole Faces List described by CGAL_point,
     * the top Level of these two Lists are Divided by all Volumes (including Outter Space) in this Nef_Polyhedron
     */
    static bool getNefPolyhedronCGALFacesByVolume(Nef_polyhedron& poly, std::vector<std::vector<std::vector<CGAL_point>>>& boundarys, std::vector<std::vector<std::vector<std::vector<CGAL_point>>>>& holes);
    
    /**
     * @brief Predicate whether a CGAL_point is Collinear with a points List described by CGAL_point
     */
    static bool isCollinear(const std::vector<CGAL_point>& pts ,const  CGAL_point& ptc);

    /**
     * @brief Predicate whether a points List described by Basic Point is ClockWise from the perspective described by a vector
     */
    static bool isClockWise(const std::vector<point3d>& pts, point3d& ptc);

    /**
     * @brief Predicate whether a points List described by CGAL_point is ClockWise from the perspective described by a vector
     */
    static bool isClockWiseCGAL(const std::vector<CGAL_point>& pts, CGAL_point& ptc);

    /**
     * @brief Predicate whether a Nef_Polyhedron is a 2D ManiFold
     */
    static bool isManiFold(Nef_polyhedron& poly);

    /**
     * @brief Predicate whether a Nef_Polyhedron is a closed 2D ManiFold
     */
    static bool isClosed(std::vector<std::vector<int> > &polyIndexs , const std::vector<point3d> &points);

    /**
     * @brief Convert a List of CGAL_point into a List of Basic Point
     */
    static bool convertCgalPointTo3d(const std::vector<CGAL_point>& originPts , std::vector<point3d>& pts);

    /**
     * @brief Convert a Nef_Polyhedron into a Surface Mesh
     */
    static bool convertNefToMesh(Nef_polyhedron& nef, Surface_mesh& mesh);

    /**
     * @brief Read a Polyhedron from one OFF file and convert it into a Nef_Polyhedron
     */
    static bool convertOFFToNef(std::ifstream& fin, Nef_polyhedron& nef);

    /**
     * @brief Read a Nef_Polyhedron from a specially formatted file
     */
    static bool LoadNefPolyhedron(std::ifstream& fin, Nef_polyhedron& nef);

    /**
     * @brief Save a Nef_Polyhedron into a specially formatted file
     */
    static bool SaveNefPolyhedron(Nef_polyhedron& nef, std::ofstream& fout);
    
    /**
     * @brief Union Operation of two Nef_Polyhedra
     */
    static bool nefPolyhedronUnion(Nef_polyhedron& nef1 , Nef_polyhedron& nef2 , Nef_polyhedron& result);

    /**
     * @brief Difference Operation of two Nef_Polyhedra
     */
    static bool nefPolyhedronDiff(Nef_polyhedron& nef1 , Nef_polyhedron& nef2 , Nef_polyhedron& result);

    /**
     * @brief Intersection Operation of two Nef_Polyhedra
     */
    static bool nefPolyhedronInter(Nef_polyhedron& nef1 , Nef_polyhedron& nef2 , Nef_polyhedron& result);
};

}

#endif
