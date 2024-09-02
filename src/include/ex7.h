#ifndef Exercise7_HH
#define Exercise7_HH

#include "FileSystem.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "polyscope/surface_mesh.h"
#include <cmath>
#include <cstddef>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <memory>

#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "polyscope/curve_network.h"
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <iostream>
#include <memory>
#include <polyscope/surface_mesh.h>


using namespace polyscope;
using namespace geometrycentral;
using namespace geometrycentral::surface;

const double inf = std::numeric_limits<double>::infinity();

class Exercise7 {
  public:
    using Point = Vector3;
    using Vec2d = Vector2;
    using T = Eigen::Triplet<double>;

  public:
    Exercise7() : gc_mesh(nullptr), geometry(nullptr), ps_mesh(nullptr) {}
    ~Exercise7() { 
        gc_mesh = nullptr;
        geometry = nullptr;
        ps_mesh = nullptr;

        tex_mesh = nullptr;
        tex_geometry = nullptr;
        ps_tex_mesh = nullptr;

        _origin = {0., 0., 0.};

    }

    void prepare_param(PolygonInstance &polygon_instance) {
        gc_mesh = polygon_instance.mesh;
        geometry = polygon_instance.geometry;
        ps_mesh = polygon_instance.ps_mesh;

        geometry->requireVertexPositions();
        geometry->requireEdgeLengths();
        geometry->requireEdgeCotanWeights();

        // compute origin and radius for tex_geometry
        // as far as i can tell, gc doesn't have bounding box info... :/
        Point bboxMin = {inf, inf, inf};
        Point bboxMax = {-inf, -inf, -inf};
        for (Vertex v : gc_mesh->vertices()) {
            Point pos = geometry->vertexPositions[v];
            if (pos.x < bboxMin.x) {
            bboxMin.x = pos.x;
            } else if (pos.y < bboxMin.y) {
            bboxMin.y = pos.y;
            } else if (pos.z < bboxMin.z) {
            bboxMin.z = pos.z;
            }
            
            if (pos[0] > bboxMax[0]) {
            bboxMax.x = pos.x;
            } else if (pos[1] > bboxMax[1]) {
            bboxMax.y = pos.y;
            } else if (pos[2] > bboxMax[2]) {
            bboxMax.z = pos.z;
            }
        }
        _radius = norm(bboxMax - bboxMin) / 10.;
        _origin = {(bboxMin.x + bboxMax.x) / 2, (bboxMin.y + bboxMax.y) / 2, bboxMin.z};
        
        // copy to new tex_mesh & tex_geometry
        delete tex_mesh;
        delete tex_geometry;
        is_mapped = false;
        tex_mesh = new surface::ManifoldSurfaceMesh(polygon_instance.mesh->getFaceVertexList());
        tex_geometry = new surface::VertexPositionGeometry(*tex_mesh);
        _tex_coords = VertexData<Vec2d>(*gc_mesh, Vec2d{0., 0.});

        tex_geometry->requireVertexPositions();
        tex_geometry->requireEdgeLengths();
    }

    void prepare_minimal(PolygonInstance &polygon_instance) {
        gc_mesh = polygon_instance.mesh;
        geometry = polygon_instance.geometry;
        ps_mesh = polygon_instance.ps_mesh;

        geometry->requireVertexPositions();
        geometry->requireEdgeCotanWeights();
    }
    int num_iter = 10;
    int smooth_type = 0;
    bool is_mapped = false;

  public:

    void map_suface_boundary_to_circle();

    void explicitly_smooth_texture();

    void implicitly_smooth_texture();

    void minimal_surface();

  private:
    Point _origin;
    double _radius;

    VertexData<Vec2d> _tex_coords;
    std::shared_ptr<surface::SurfaceMesh> gc_mesh;
    std::shared_ptr<surface::VertexPositionGeometry> geometry;
    polyscope::SurfaceMesh *ps_mesh;

    surface::ManifoldSurfaceMesh* tex_mesh;
    surface::VertexPositionGeometry* tex_geometry;
    polyscope::SurfaceMesh* ps_tex_mesh;
};

#endif
