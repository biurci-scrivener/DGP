#ifndef OPENFLIPPER_CURVESMOOTHING_HH
#define OPENFLIPPER_CURVESMOOTHING_HH

#include <geometrycentral/surface/vertex_position_geometry.h>
#include <polyscope/surface_mesh.h>

#include "FileSystem.h"
#include "geometrycentral/surface/surface_mesh.h"

using namespace geometrycentral;
using namespace polyscope;

class Exercise4 {
   public:
    using Point = Vector3;

    Exercise4()
        : gc_mesh(nullptr), geometry(nullptr), ps_mesh(nullptr) {}
    ~Exercise4() {
        gc_mesh = nullptr;
        geometry = nullptr;
        ps_mesh = nullptr;
    }

   public:
    void set_mesh(PolygonInstance& polygon_instance) {
        gc_mesh = polygon_instance.mesh;
        geometry = polygon_instance.geometry;

        ps_mesh = polygon_instance.ps_mesh;
        num_vertices = gc_mesh->nVertices();

        vertex_normals.resize(num_vertices);
        vertex_curvature_weights.resize(num_vertices);

        geometry->requireVertexDualAreas();
        geometry->requireEdgeCotanWeights();
        geometry->requireVertexNormals();
        geometry->requireVertexPositions();
        geometry->requireVertexAngleSums();
    }
    int normal_type = 0;
    int curvature_type = 0;

    void show_normal();

    void show_curvature();

    void compute_normals_with_constant_weights();

    void compute_normals_by_area_weights();

    void compute_normals_with_angle_weights();

    void calc_uniform_laplacian();

    void calc_mean_curvature();

    void calc_gauss_curvature();

    void color_coding_d(const std::vector<double> data,
                        const double _min_value = 0.0,
                        const double _max_value = 0.0, const int _bound = 20);

    double min_curvature = 0;
    double max_curvature = 0;

   private:
    size_t num_vertices;

    std::vector<std::array<double, 3>> vertex_normals;
    std::vector<double> vertex_curvature_weights;

    std::shared_ptr<surface::SurfaceMesh> gc_mesh;
    std::shared_ptr<surface::VertexPositionGeometry> geometry;
    polyscope::SurfaceMesh* ps_mesh;
};
#endif
