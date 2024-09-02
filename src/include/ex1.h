#ifndef REMB_EIGENTUTORIAL_HH
#define REMB_EIGENTUTORIAL_HH

#include <iostream>
#include "Eigen/Sparse"
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <polyscope/surface_mesh.h>

#include "FileSystem.h"
#include "geometrycentral/surface/surface_mesh.h"

using namespace geometrycentral;
using namespace polyscope;

    class Exercise1 {
    public:
        using Point = Vector3;
        using SMatd = Eigen::SparseMatrix<double>;
        using T = Eigen::Triplet<double>;
        using VecXd = Eigen::VectorXd;

        Exercise1()
            : gc_mesh(nullptr), geometry(nullptr), ps_mesh(nullptr) {}
        ~Exercise1() {
            gc_mesh = nullptr;
            geometry = nullptr;
            ps_mesh = nullptr;
        }

    public:
        void set_mesh(PolygonInstance& polygon_instance) {
            gc_mesh = polygon_instance.mesh;
            geometry = polygon_instance.geometry;
            ps_mesh = polygon_instance.ps_mesh;
        }
        void solve_sparse_linear_system();
        void generate_plot(int function_id);

    private:
        std::shared_ptr<surface::SurfaceMesh> gc_mesh;
        std::shared_ptr<surface::VertexPositionGeometry> geometry;
        polyscope::SurfaceMesh* ps_mesh;
};

#endif 
