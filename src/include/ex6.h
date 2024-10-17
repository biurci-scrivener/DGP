#ifndef SMOOTHING_HH
#define SMOOTHING_HH

#include "FileSystem.h"
#include "polyscope/surface_mesh.h"
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <memory>
#include <polyscope/render/color_maps.h>
#include "Eigen/Sparse"

using namespace geometrycentral;
using namespace polyscope;

class Exercise6 {
  public:
    using Point = Vector3;
    using T = Eigen::Triplet<double>;

    Exercise6() : gc_mesh(nullptr), geometry(nullptr) {}
    ~Exercise6() {
        gc_mesh = nullptr;
        geometry = nullptr;
        ps_mesh = nullptr;
    }
    void prepare(PolygonInstance& polygon_instance) {
        gc_mesh = polygon_instance.mesh;
        geometry = polygon_instance.geometry;
        ps_mesh = polygon_instance.ps_mesh;

        geometry->requireVertexPositions();
        geometry->requireEdgeCotanWeights();
        geometry->requireVertexDualAreas();
        geometry->requireVertexIndices();
    }
    int smooth_type = 0;
    int num_iters = 10;
    double coefficient = 2.0;
    double timestep = 1e-5;
    void smooth() {
        switch (smooth_type) {
        case 0:
            uniform_smooth();
            break;
        case 1:
            cotan_laplacian_smooth();
            break;
        case 2:
            implicit_smooth();
            break;
        default:
            break;
        }
    }
    void enhance_feature() {
        switch (smooth_type) {
        case 0:
            uniform_laplacian_enhance_feature();
            break;
        case 1:
            cotan_laplacian_enhance_feature();
            break;
        default:
            break;
        }
    }
    void uniform_smooth();

    void cotan_laplacian_smooth();

    void implicit_smooth();

    void uniform_laplacian_enhance_feature();

    void cotan_laplacian_enhance_feature();

    void deform_surface(std::vector<int>& fixed_verts, const std::vector<int>& displaced_verts,
                         const Point& displacement_vector);

  private:
    std::shared_ptr<surface::SurfaceMesh> gc_mesh;
    std::shared_ptr<surface::VertexPositionGeometry> geometry;
    polyscope::SurfaceMesh* ps_mesh;
};
#endif
