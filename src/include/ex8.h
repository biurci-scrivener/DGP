#ifndef Exercise8_HH
#define Exercise8_HH

#include <geometrycentral/surface/vertex_position_geometry.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>

#include "FileSystem.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "Eigen/Sparse"

#include <queue>
#include <stdexcept>

using namespace geometrycentral;
using namespace polyscope;

class Exercise8 {
   public:
    using Point = Vector3;

    Exercise8()
        : gc_mesh(nullptr), geometry(nullptr), ps_mesh(nullptr) {}
    ~Exercise8() {
        gc_mesh = nullptr;
        geometry = nullptr;
        ps_mesh = nullptr;
    }

   public:
    void set_mesh(PolygonInstance& polygon_instance) {

        if (!polygon_instance.mesh->isManifold()) {
            polyscope::warning("ERROR: The current mesh is non-manifold");
            return;
        }

        gc_mesh = polygon_instance.mesh;
        geometry = polygon_instance.geometry;
        ps_mesh = polygon_instance.ps_mesh;
        ps_mesh->setSmoothShade(true);
        
        clear_quantities();

        geometry->requireCotanLaplacian();
        geometry->requireVertexLumpedMassMatrix();
        geometry->requireVertexDualAreas();
        geometry->requireEdgeCotanWeights();
        geometry->requireEdgeLengths();
        geometry->requireVertexNormals();
        geometry->requireFaceNormals();
        geometry->requireVertexPositions();
        geometry->requireDECOperators();

        initialized = true;
    }

    void clear_quantities() {

        polyscope::removeStructure("Tree");
        polyscope::removeStructure("Cotree");
        polyscope::removeStructure("Generators, init. edges");
        for (size_t i = 0; i < generators.size(); i++) polyscope::removeStructure("Generator " + std::to_string(i));

        ps_mesh->removeQuantity("Omega");
        ps_mesh->removeQuantity("Alpha");
        ps_mesh->removeQuantity("Beta");
        ps_mesh->removeQuantity("d-Alpha");
        ps_mesh->removeQuantity("delta-Beta");
        ps_mesh->removeQuantity("Gamma");
        for (size_t i = 0; i < generators.size(); i++) ps_mesh->removeQuantity("(Basis) Omega " + std::to_string(i));
        for (size_t i = 0; i < generators.size(); i++) ps_mesh->removeQuantity("(Basis) Gamma " + std::to_string(i));

        vertexParent.clear();
        faceParent.clear();
        generators_init.clear();
        generators.clear();
        basis_omega.clear();
        basis_harmonic.clear();
        omega_initialized = false;

    }

    void show_hodge(int component) {

        // Hodge Decomp.

        auto omega_vis = ps_mesh->getQuantity("Omega");
        if (omega_vis != nullptr) omega_vis->setEnabled(component == 0 ? true : false);

        auto alpha_vis = ps_mesh->getQuantity("Alpha");
        if (alpha_vis != nullptr) alpha_vis->setEnabled(component == 1 ? true : false);

        auto beta_vis = ps_mesh->getQuantity("Beta");
        if (beta_vis != nullptr) beta_vis->setEnabled(component == 2 ? true : false);

        auto d_alpha_vis = ps_mesh->getQuantity("d-Alpha");
        if (d_alpha_vis != nullptr) d_alpha_vis->setEnabled(component == 1 ? true : false);

        auto delta_beta_vis = ps_mesh->getQuantity("delta-Beta");
        if (delta_beta_vis != nullptr) delta_beta_vis->setEnabled(component == 2 ? true : false);

        auto gamma_vis = ps_mesh->getQuantity("Gamma");
        if (gamma_vis != nullptr) gamma_vis->setEnabled(component == 3 ? true : false);

        // Homology Generators

        if (polyscope::hasCurveNetwork("Tree")) {
            auto tree_vis = polyscope::getCurveNetwork("Tree");
            tree_vis->setEnabled(false);
        }

        if (polyscope::hasCurveNetwork("Cotree")) {
            auto cotree_vis = polyscope::getCurveNetwork("Cotree");
            cotree_vis->setEnabled(false);
        }

        if (polyscope::hasCurveNetwork("Generators, init. edges")) {
            auto gen_init_vis = polyscope::getCurveNetwork("Generators, init. edges");
            gen_init_vis->setEnabled(false);
        }
        
        for (size_t i = 0; i < generators.size(); i++) {

            if (polyscope::hasCurveNetwork("Generator " + std::to_string(i))) {
                auto gen_vis = polyscope::getCurveNetwork("Generator " + std::to_string(i));
                gen_vis->setEnabled(false);
            }

            auto omega_bas_vis = ps_mesh->getQuantity("(Basis) Omega " + std::to_string(i));
            if (omega_bas_vis != nullptr) omega_bas_vis->setEnabled(false);
            auto gamma_bas_vis = ps_mesh->getQuantity("(Basis) Gamma " + std::to_string(i));
            if (gamma_bas_vis != nullptr) gamma_bas_vis->setEnabled(false);
        }

    } 

    void show_tree_cotree() {

        // Hodge Decomp.

        auto omega_vis = ps_mesh->getQuantity("Omega");
        if (omega_vis != nullptr) omega_vis->setEnabled(false);

        auto alpha_vis = ps_mesh->getQuantity("Alpha");
        if (alpha_vis != nullptr) alpha_vis->setEnabled(false);

        auto beta_vis = ps_mesh->getQuantity("Beta");
        if (beta_vis != nullptr) beta_vis->setEnabled(false);

        auto d_alpha_vis = ps_mesh->getQuantity("d-Alpha");
        if (d_alpha_vis != nullptr) d_alpha_vis->setEnabled(false);

        auto delta_beta_vis = ps_mesh->getQuantity("delta-Beta");
        if (delta_beta_vis != nullptr) delta_beta_vis->setEnabled(false);

        auto gamma_vis = ps_mesh->getQuantity("Gamma");
        if (gamma_vis != nullptr) gamma_vis->setEnabled(false);

        // Homology Generators

        if (polyscope::hasCurveNetwork("Tree")) {
            auto tree_vis = polyscope::getCurveNetwork("Tree");
            tree_vis->setEnabled(true);
        }

        if (polyscope::hasCurveNetwork("Cotree")) {
            auto cotree_vis = polyscope::getCurveNetwork("Cotree");
            cotree_vis->setEnabled(true);
        }

        if (polyscope::hasCurveNetwork("Generators, init. edges")) {
            auto gen_init_vis = polyscope::getCurveNetwork("Generators, init. edges");
            gen_init_vis->setEnabled(false);
        }
        
        for (size_t i = 0; i < generators.size(); i++) {

            if (polyscope::hasCurveNetwork("Generator " + std::to_string(i))) {
                auto gen_vis = polyscope::getCurveNetwork("Generator " + std::to_string(i));
                gen_vis->setEnabled(false);
            }

            auto omega_bas_vis = ps_mesh->getQuantity("(Basis) Omega " + std::to_string(i));
            if (omega_bas_vis != nullptr) omega_bas_vis->setEnabled(false);
            auto gamma_bas_vis = ps_mesh->getQuantity("(Basis) Gamma " + std::to_string(i));
            if (gamma_bas_vis != nullptr) gamma_bas_vis->setEnabled(false);
        }

    } 

    void show_generators() {

        // Hodge Decomp.

        auto omega_vis = ps_mesh->getQuantity("Omega");
        if (omega_vis != nullptr) omega_vis->setEnabled(false);

        auto alpha_vis = ps_mesh->getQuantity("Alpha");
        if (alpha_vis != nullptr) alpha_vis->setEnabled(false);

        auto beta_vis = ps_mesh->getQuantity("Beta");
        if (beta_vis != nullptr) beta_vis->setEnabled(false);

        auto d_alpha_vis = ps_mesh->getQuantity("d-Alpha");
        if (d_alpha_vis != nullptr) d_alpha_vis->setEnabled(false);

        auto delta_beta_vis = ps_mesh->getQuantity("delta-Beta");
        if (delta_beta_vis != nullptr) delta_beta_vis->setEnabled(false);

        auto gamma_vis = ps_mesh->getQuantity("Gamma");
        if (gamma_vis != nullptr) gamma_vis->setEnabled(false);

        // Homology Generators

        if (polyscope::hasCurveNetwork("Tree")) {
            auto tree_vis = polyscope::getCurveNetwork("Tree");
            tree_vis->setEnabled(false);
        }

        if (polyscope::hasCurveNetwork("Cotree")) {
            auto cotree_vis = polyscope::getCurveNetwork("Cotree");
            cotree_vis->setEnabled(false);
        }

        if (polyscope::hasCurveNetwork("Generators, init. edges")) {
            auto gen_init_vis = polyscope::getCurveNetwork("Generators, init. edges");
            gen_init_vis->setEnabled(false);
        }
        
        for (size_t i = 0; i < generators.size(); i++) {

            if (polyscope::hasCurveNetwork("Generator " + std::to_string(i))) {
                auto gen_vis = polyscope::getCurveNetwork("Generator " + std::to_string(i));
                gen_vis->setEnabled(true);
            }

            auto omega_bas_vis = ps_mesh->getQuantity("(Basis) Omega " + std::to_string(i));
            if (omega_bas_vis != nullptr) omega_bas_vis->setEnabled(true);
            auto gamma_bas_vis = ps_mesh->getQuantity("(Basis) Gamma " + std::to_string(i));
            if (gamma_bas_vis != nullptr) gamma_bas_vis->setEnabled(false);
        }

    } 

    void show_generators_init() {

        // Hodge Decomp.

        auto omega_vis = ps_mesh->getQuantity("Omega");
        if (omega_vis != nullptr) omega_vis->setEnabled(false);

        auto alpha_vis = ps_mesh->getQuantity("Alpha");
        if (alpha_vis != nullptr) alpha_vis->setEnabled(false);

        auto beta_vis = ps_mesh->getQuantity("Beta");
        if (beta_vis != nullptr) beta_vis->setEnabled(false);

        auto d_alpha_vis = ps_mesh->getQuantity("d-Alpha");
        if (d_alpha_vis != nullptr) d_alpha_vis->setEnabled(false);

        auto delta_beta_vis = ps_mesh->getQuantity("delta-Beta");
        if (delta_beta_vis != nullptr) delta_beta_vis->setEnabled(false);

        auto gamma_vis = ps_mesh->getQuantity("Gamma");
        if (gamma_vis != nullptr) gamma_vis->setEnabled(false);

        // Homology Generators

        if (polyscope::hasCurveNetwork("Tree")) {
            auto tree_vis = polyscope::getCurveNetwork("Tree");
            tree_vis->setEnabled(true);
        }

        if (polyscope::hasCurveNetwork("Cotree")) {
            auto cotree_vis = polyscope::getCurveNetwork("Cotree");
            cotree_vis->setEnabled(true);
        }

        if (polyscope::hasCurveNetwork("Generators, init. edges")) {
            auto gen_init_vis = polyscope::getCurveNetwork("Generators, init. edges");
            gen_init_vis->setEnabled(true);
        }
        
        for (size_t i = 0; i < generators.size(); i++) {

            if (polyscope::hasCurveNetwork("Generator " + std::to_string(i))) {
                auto gen_vis = polyscope::getCurveNetwork("Generator " + std::to_string(i));
                gen_vis->setEnabled(false);
            }

            auto omega_bas_vis = ps_mesh->getQuantity("(Basis) Omega " + std::to_string(i));
            if (omega_bas_vis != nullptr) omega_bas_vis->setEnabled(false);
            auto gamma_bas_vis = ps_mesh->getQuantity("(Basis) Gamma " + std::to_string(i));
            if (gamma_bas_vis != nullptr) gamma_bas_vis->setEnabled(false);
        }

    } 

    void show_harmonic_basis(int idx) {

        // Hodge Decomp.

        auto omega_vis = ps_mesh->getQuantity("Omega");
        if (omega_vis != nullptr) omega_vis->setEnabled(false);

        auto alpha_vis = ps_mesh->getQuantity("Alpha");
        if (alpha_vis != nullptr) alpha_vis->setEnabled(false);

        auto beta_vis = ps_mesh->getQuantity("Beta");
        if (beta_vis != nullptr) beta_vis->setEnabled(false);

        auto d_alpha_vis = ps_mesh->getQuantity("d-Alpha");
        if (d_alpha_vis != nullptr) d_alpha_vis->setEnabled(false);

        auto delta_beta_vis = ps_mesh->getQuantity("delta-Beta");
        if (delta_beta_vis != nullptr) delta_beta_vis->setEnabled(false);

        auto gamma_vis = ps_mesh->getQuantity("Gamma");
        if (gamma_vis != nullptr) gamma_vis->setEnabled(false);

        // Homology Generators

        if (polyscope::hasCurveNetwork("Tree")) {
            auto tree_vis = polyscope::getCurveNetwork("Tree");
            tree_vis->setEnabled(false);
        }

        if (polyscope::hasCurveNetwork("Cotree")) {
            auto cotree_vis = polyscope::getCurveNetwork("Cotree");
            cotree_vis->setEnabled(false);
        }

        if (polyscope::hasCurveNetwork("Generators, init. edges")) {
            auto gen_init_vis = polyscope::getCurveNetwork("Generators, init. edges");
            gen_init_vis->setEnabled(false);
        }
        
        for (size_t i = 0; i < generators.size(); i++) {

            if (polyscope::hasCurveNetwork("Generator " + std::to_string(i))) {
                auto gen_vis = polyscope::getCurveNetwork("Generator " + std::to_string(i));
                gen_vis->setEnabled(((unsigned long)idx == i) ? true : false);
            }

            auto omega_bas_vis = ps_mesh->getQuantity("(Basis) Omega " + std::to_string(i));
            if (omega_bas_vis != nullptr) omega_bas_vis->setEnabled(false);
            auto gamma_bas_vis = ps_mesh->getQuantity("(Basis) Gamma " + std::to_string(i));
            if (gamma_bas_vis != nullptr) gamma_bas_vis->setEnabled(((idx == -1) || ((unsigned long)idx == i)) ? true : false);
        }

    }

    int computeEulerCharacteristic() {
        return gc_mesh->nVertices() - gc_mesh->nEdges() + gc_mesh->nFaces();
    }

    int nConnectedComponents() {
        return gc_mesh->nConnectedComponents();
    }

    int hasTree() {
        return (vertexParent.size() != 0);
    }

    void hodge_decomposition();
    
   
    void build_dual_spanning_cotree();
    bool in_dual_spanning_cotree(Halfedge he);
    void build_primal_spanning_tree();
    bool in_primal_spanning_tree(Halfedge he);

    std::vector<Halfedge> trace_to_cotree_root(Face f);
    void build_generators();
    void generate_harmonic_basis();
    void harmonic_decomposition();

    // ===

    Halfedge sharedHalfedge(Face f, Face g);
    void vis_hodge_decomposition();
    void buildTreeCotreeGraph();
    void buildGeneratorsGraph();
    void buildHarmonicBasisVis();
    void generate_random_one_form();
    Eigen::VectorXd solve_poisson(const Eigen::VectorXd &rho);
    Point centroid(Face f);
    bool initialized = false;
    bool omega_initialized = false;

    private:

    double VECTOR_SCALE_DEFAULT = 0.05;

    // === Tree-Cotree
    std::map<Vertex, Vertex> vertexParent;
    std::map<Face, Face> faceParent;
    std::vector<Halfedge> generators_init;
    std::vector<std::vector<Halfedge>> generators;
    std::vector<Eigen::VectorXd> basis_omega;
    std::vector<Eigen::VectorXd> basis_harmonic;

    // === Hodge Decomp.
    Eigen::VectorXd omega;
    Eigen::VectorXd alpha;
    Eigen::VectorXd beta;
    Eigen::VectorXd d_alpha;
    Eigen::VectorXd delta_beta;
    Eigen::VectorXd gamma;

    std::shared_ptr<surface::SurfaceMesh> gc_mesh;
    std::shared_ptr<surface::VertexPositionGeometry> geometry;
    polyscope::SurfaceMesh* ps_mesh;

};
#endif
