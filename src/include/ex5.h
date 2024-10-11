#ifndef NUMERICAL_METHODS_HH
#define NUMERICAL_METHODS_HH

#include <iostream>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <filesystem>
#include "Eigen/Sparse"
#include "FileSystem.h"
#include <fast_matrix_market/app/Eigen.hpp>
#include "polyscope/point_cloud.h"
#include "polyscope/image_quantity.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/surface/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace polyscope;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using geometrycentral::pointcloud::PointPositionGeometry;
using geometrycentral::pointcloud::PointData;

class Exercise5 {
    public:
        using SMatd = Eigen::SparseMatrix<double>;
        using T = Eigen::Triplet<double>;
        using VecXd = Eigen::VectorXd;
        using Point = glm::vec3;
        using PCPointCloud = geometrycentral::pointcloud::PointCloud;

        Exercise5()
            : gc_mesh(nullptr), geometry(nullptr), ps_mesh(nullptr) {}
        ~Exercise5() {
            gc_mesh = nullptr;
            geometry = nullptr;
            ps_mesh = nullptr;
        }

        void set_mesh(PolygonInstance& polygon_instance) {
            gc_mesh = polygon_instance.mesh;
            geometry = polygon_instance.geometry;

            ps_mesh = polygon_instance.ps_mesh;

            geometry->requireVertexPositions();
        }

        void load_sparse_matrices() {

            if (std::filesystem::exists(file_path_base + "matrix1.mtx") && std::filesystem::exists(file_path_base + "matrix1.mtx")) {
                std::ifstream file_stream_1(file_path_base + "matrix1.mtx");
                std::ifstream file_stream_2(file_path_base + "matrix2.mtx");
                
                fast_matrix_market::read_matrix_market_eigen(file_stream_1, matrix_1);
                fast_matrix_market::read_matrix_market_eigen(file_stream_2, matrix_2);
                mtx_loaded = true;
                std::cout << "Loaded matrix1.mtx and matrix2.mtx." << std::endl;
            } else {
                polyscope::warning("ERROR: Couldn't load matrix1.mtx and matrix2.mtx. Run DGP from within the bin folder.");
                std::cout << "" << std::endl;
            }
            
        }

        void compare_dense_solvers();
        void analyze_sparsity_pattern();
        void compare_sparse_solvers();
        
        void find_all_eigenvectors_eigenvalues();
        void find_closest_eigenvalues();

        void gradient_descent(int func_id, double step_size);
        void newtons_method(int func_id);

        //==========================================

        void plot_function(int func_id);

        void generate_world_axes();
        void generate_ellipsoid_axes(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c);
        void p2_generate_iterate_plot(std::vector<Eigen::Vector3d> v1_it, std::vector<Eigen::Vector3d> v2_it, std::vector<Eigen::Vector3d> v3_it);

        void generate_ellipse(std::string name, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, double a_len, double b_len, double c_len);

        void generate_sparsity_plot(SMatd A);

        void p3_generate_iterate_plot(std::vector<Eigen::Vector3d> iterates, int func_id);

        void polyscope_plot_vector(std::string name, Point pos, Eigen::Vector3d dir, Point color);
        void polyscope_plot_vectors(std::string name, Point pos, std::vector<Eigen::Vector3d> dirs, std::vector<Point> colors);
        void polyscope_plot_vectors(std::string name, Point pos, std::vector<Eigen::Vector3d> dirs, Point base_color);
        void polyscope_plot_vectors(std::string name, std::vector<Vector3> pos, std::vector<Eigen::Vector3d> dirs, Point color);

        bool mtx_loaded = false;

    private:

        double PLOTTING_SCALE_f = 5e-1;
        double PLOTTING_SCALE_g = 1e-1;

        std::string file_path_base = "../matrices/";

        std::shared_ptr<surface::SurfaceMesh> gc_mesh;
        std::shared_ptr<surface::VertexPositionGeometry> geometry;
        polyscope::SurfaceMesh* ps_mesh;

        SMatd matrix_1;
        SMatd matrix_2;

};

#endif 
