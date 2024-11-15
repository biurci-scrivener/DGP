#include "include/ex8.h"


/*
    Part 1: Hodge Decomposition

    Populate the following member variables as part of your implementation:
    - alpha (function on vertices)
    - beta (function on faces)
    - d_alpha (one-form)
    - delta_beta (one-form)
    - gamma (one-form)
*/


void Exercise8::hodge_decomposition() {
    
    // lazy way to make sure the mesh has no boundary without initializing a new ManifoldSurfaceMesh
    for (Edge e: gc_mesh->edges()) {
        if (e.isBoundary()) {
            polyscope::warning("ERROR: Please load a mesh with no boundary.");
            return;
        }
    }

    // TODO

}

/*
    Part 2: Tree-Cotree
*/

void Exercise8::build_dual_spanning_cotree() {

    /* 
    * TODO:
    * Populate the member variable <faceParent>, which is a std::map that maps each face 
    * of the input mesh to its parent in the dual spanning tree.
    */

}

bool Exercise8::in_dual_spanning_cotree(Halfedge he) {

    /*
        TODO:
        Return "True" if either he or he.twin() is in cotree.
    */

}

void Exercise8::build_primal_spanning_tree() {

    /* 
    * Populate the member variable <vertexParent>, which is a std::map that maps each vertex 
    * of the input mesh to its parent in the primal spanning tree.
    */

}

bool Exercise8::in_primal_spanning_tree(Halfedge he) {

    /*
        TODO:
        Return "True" if either he or he.twin() is in tree.
    */

}

std::vector<Halfedge> Exercise8::trace_to_cotree_root(Face f) {
    /*
        TODO:
        Returns a vector of halfedges
        corresponding to the path from a face to the root of the cotree.
    */

    std::vector<Halfedge> path;

    // TODO

    return path;
}

void Exercise8::build_generators() {

    generators_init.clear();
    generators.clear();

    /*
        TODO:
        Populate generators_init and generators.
    */ 
    

}

void Exercise8::generate_harmonic_basis() {

    basis_omega.clear();
    basis_harmonic.clear();

    /*  
        For a particular edge:
        if True: canonical orientation is lower idx to higher idx
        otherwise, higher to lower
    */ 
    EdgeData<char> e_orient = polyscopeEdgeOrientations(*gc_mesh);

    /*
        TODO:
        Populate basis_omega with a closed 1-form for each generator
        and basis_harmonic with each resulting harmonic component.
    */
    
}

void Exercise8::harmonic_decomposition() {

    if (basis_harmonic.size() != (unsigned long)(-(computeEulerCharacteristic() - 2))) {
        std::cout << "Can't run harmonic_decomposition(): basis_harmonic does not contain 2g elements" << std::endl;
        return;
    } else if ((size_t)gamma.size() != gc_mesh->nEdges()) {
        std::cout << "Can't run harmonic_decomposition() until hodge_decomposition() has been run." << std::endl;
        return;
    }

    /*
        TODO: Test your harmonic basis by finding coefficients for the
        harmonic part extracted by hodge_decomposition().

        Print the rank of the RHS, the solution to the linear system, and
        ||gamma_appox - gamma||, which should be close to 0.
    */


}

// === UTILITY FUNCTIONS BELOW ===

Halfedge Exercise8::sharedHalfedge(Face f, Face g) {

    // will always return the halfedge on "f"

    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // Should never get here!
    std::cerr << "Oops, TreeCotree::sharedHalfedge() received bad input." << std::endl;
    return f.halfedge();
}

Exercise8::Point Exercise8::centroid(Face f) {

    Halfedge he = f.halfedge();
    Point a = geometry->inputVertexPositions[he.vertex()];
    Point b = geometry->inputVertexPositions[he.next().vertex()];
    Point c = geometry->inputVertexPositions[he.next().next().vertex()];
    return (a + b + c) / 3.0;

}

void Exercise8::vis_hodge_decomposition() {

    auto vis_omega = ps_mesh->getQuantity("Omega");
    if (!(vis_omega == nullptr)) {
        vis_omega->setEnabled(false);
    }

    if ((size_t)alpha.size() == gc_mesh->nVertices()) {
        auto vis_a_ = ps_mesh->addVertexScalarQuantity("Alpha", alpha);
    } else {
        std::cout << "ERROR: beta has wrong size " << alpha.size() << ". Maybe it hasn't been set yet?" << std::endl;
    }
    
    if ((size_t)beta.size() == gc_mesh->nFaces()) {
        auto vis_b_ = ps_mesh->addFaceScalarQuantity("Beta", beta);
    } else {
        std::cout << "ERROR: beta has wrong size " << beta.size() << ". Maybe it hasn't been set yet?" << std::endl;
    }

    if ((size_t)d_alpha.size() == gc_mesh->nEdges()) {
        auto vis_a = ps_mesh->addOneFormTangentVectorQuantity("d-Alpha", d_alpha, polyscopeEdgeOrientations(*gc_mesh));
        vis_a->setEnabled(true);
        vis_a->setVectorLengthScale(VECTOR_SCALE_DEFAULT);
    } else {
        std::cout << "ERROR: d-Alpha has wrong size " << d_alpha.size() << ". Maybe it hasn't been set yet?" << std::endl;
    }
    
    if ((size_t)delta_beta.size() == gc_mesh->nEdges()) {
        auto vis_b = ps_mesh->addOneFormTangentVectorQuantity("delta-Beta", delta_beta, polyscopeEdgeOrientations(*gc_mesh));
        vis_b->setEnabled(true);
        vis_b->setVectorLengthScale(VECTOR_SCALE_DEFAULT);
    } else {
        std::cout << "ERROR: delta-Beta has wrong size " << delta_beta.size() << ". Maybe it hasn't been set yet?" << std::endl;
    }

    if ((size_t)gamma.size() == gc_mesh->nEdges()) {
        auto vis_g = ps_mesh->addOneFormTangentVectorQuantity("Gamma", gamma, polyscopeEdgeOrientations(*gc_mesh));
        vis_g->setEnabled(true);
        vis_g->setVectorLengthScale(VECTOR_SCALE_DEFAULT);
    } else {
        std::cout << "ERROR: gamma has wrong size " << gamma.size() << ". Maybe it hasn't been set yet?" << std::endl;
    }

}

void Exercise8::buildTreeCotreeGraph() {

    double meanEdgeLength = 0.;
    for (Edge e : gc_mesh->edges()) {
        meanEdgeLength += geometry->edgeLengths[e];
    }
    meanEdgeLength /= gc_mesh->nEdges();
    double lengthScale = meanEdgeLength;

    // Build primal spanning tree.
    std::vector<Vector3> vertPrimal;
    std::vector<std::array<size_t, 2>> indPrimal;
    if (vertexParent.size() > 0) {
        for (Vertex v : gc_mesh->vertices()) {
            auto vp = vertexParent.find(v);
            if (vp == vertexParent.end()) {
                throw std::runtime_error("Error: Vertex " + std::to_string(v) + " not found in tree");
            }
            vertPrimal.push_back(geometry->inputVertexPositions[v]);
            indPrimal.push_back({v.getIndex(), vp->second.getIndex()});
        }
    }

    // Build dual spanning tree.
    std::vector<Vector3> vertDual;
    std::vector<std::array<size_t, 2>> indDual;
    if (faceParent.size() > 0) {
        for (Face f : gc_mesh->faces()) {
            auto fp = faceParent.find(f);
            if (fp == faceParent.end()) {
                throw std::runtime_error("Error: Face " + std::to_string(f) + " not found in cotree");
            }
            vertDual.push_back(centroid(f));
            indDual.push_back({f.getIndex(), fp->second.getIndex()});
        }
    }

    auto treeGraph = polyscope::registerCurveNetwork("Tree", vertPrimal, indPrimal);
    auto cotreeGraph = polyscope::registerCurveNetwork("Cotree", vertDual, indDual);
    treeGraph->setColor({0, 0, 1});
    cotreeGraph->setColor({0, 1, 0});
    treeGraph->setRadius(lengthScale * 0.04);
    cotreeGraph->setRadius(lengthScale * 0.04);
    treeGraph->setEnabled(true);
    cotreeGraph->setEnabled(true);
}

void Exercise8::buildGeneratorsGraph() {

    
    for (size_t i = 0; i < generators.size(); i++) {
        std::vector<Vector3> verts;
        std::vector<std::array<size_t, 2>> inds;
        size_t offset = 0;
        size_t M = generators[i].size();
        // Add interior vertices and edges
        for (size_t j = 0; j < M; j++) {
            Halfedge he = generators[i][j];
            Vector3 a = geometry->inputVertexPositions[he.tailVertex()];
            Vector3 b = geometry->inputVertexPositions[he.tipVertex()];
            verts.push_back(centroid(he.face()));
            verts.push_back((a + b) * 0.5);
            verts.push_back(centroid(he.twin().face()));
            inds.push_back({offset + 3 * j, offset + 3 * j + 1});
            inds.push_back({offset + 3 * j + 1, offset + 3 * j + 2});
        }
        offset += 3 * M;
        auto generatorsGraph = polyscope::registerCurveNetwork("Generator " + std::to_string(i), verts, inds);
        generatorsGraph->setEnabled(true);
    }

    std::vector<Vector3> verts;
    std::vector<std::array<size_t, 2>> inds;
    for (size_t i = 0; i < generators_init.size(); i++) {
        Halfedge he = generators_init[i];
        Vector3 a = geometry->inputVertexPositions[he.tailVertex()];
        Vector3 b = geometry->inputVertexPositions[he.tipVertex()];
        verts.push_back(centroid(he.face()));
        verts.push_back((a + b) * 0.5);
        verts.push_back(centroid(he.twin().face()));
        inds.push_back({3 * i, 3 * i + 1});
        inds.push_back({3 * i + 1, 3 * i + 2});
    }
    auto generatorsInitGraph = polyscope::registerCurveNetwork("Generators, init. edges", verts, inds);
    generatorsInitGraph->setEnabled(true);
    generatorsInitGraph->setColor({1.,0.,0.});
    
}

void Exercise8::buildHarmonicBasisVis() {

    for (size_t i = 0; i < basis_omega.size(); i++) {
        auto vis_omega = ps_mesh->addOneFormTangentVectorQuantity("(Basis) Omega " + std::to_string(i), basis_omega[i], polyscopeEdgeOrientations(*gc_mesh));
        vis_omega->setEnabled(true);
        vis_omega->setVectorLengthScale(VECTOR_SCALE_DEFAULT);
    }
    
    for (size_t i = 0; i < basis_omega.size(); i++) {
        auto vis_gamma = ps_mesh->addOneFormTangentVectorQuantity("(Basis) Gamma " + std::to_string(i), basis_harmonic[i], polyscopeEdgeOrientations(*gc_mesh));
        vis_gamma->setEnabled(true);
        vis_gamma->setVectorLengthScale(VECTOR_SCALE_DEFAULT);
    }
    
}

Eigen::VectorXd Exercise8::solve_poisson(const Eigen::VectorXd &rho) {

    Eigen::SparseMatrix<double> L = geometry->cotanLaplacian;
    Eigen::SparseMatrix<double> M = geometry->vertexLumpedMassMatrix;
    Eigen::VectorXd rho_bar = Eigen::VectorXd::Constant(rho.rows(), rho.mean());

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L);
    Eigen::VectorXd res = solver.solve(-M * (rho - rho_bar));

    return res;

}

void Exercise8::generate_random_one_form() {

    // Generate random vector and scalar potentials.
    int V = gc_mesh->nVertices();
    int n = std::max(2, (int)(V / 1000));
    Eigen::VectorXd rho1 = Eigen::VectorXd::Zero(V);
    Eigen::VectorXd rho2 = Eigen::VectorXd::Zero(V);
    for (int i = 0; i < n; i++) {
        rho1[rand() % V] = rand() % 1000 - 500;
        rho2[rand() % V] = rand() % 1000 - 500;
    }
    Eigen::VectorXd scalarPotential = solve_poisson(rho1);
    Eigen::VectorXd vectorPotential = solve_poisson(rho2);

    // Compute per-face field
    std::map<Face, Vector3> field;
    Eigen::VectorXd random_v = Eigen::VectorXd::Random(3);
    Vector3 up = {0.,1.,0.};
    for (Face f : gc_mesh->faces()) {
        double A = geometry->faceArea(f);
        Vector3 N = geometry->faceNormal(f);
        Vector3 C = centroid(f);
        field[f] = {0, 0, 0};

        // Add exact and coexact components
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t i = he.next().tipVertex().getIndex();
            Vector3 e = geometry->halfedgeVector(he);
            Vector3 eT = cross(N, e);
            field[f] += eT * scalarPotential[i] / (2 * A);
            field[f] += e * vectorPotential[i] / (2 * A);
        }

        // Add harmonic component.
        Vector3 e2 = {random_v[0], random_v[1], random_v[2]};
        e2 = e2.normalize();
        Vector3 e1 = cross(up, e2).normalize();
        Vector3 e3 = cross(e1, e2).normalize();

        C -= e2 * dot(C, e2); // project
        C = e1 * (-dot(C, e3)) + e3 * (dot(C, e1)); // rotate
        C -= N * dot(C, N);
        C = C.normalize();
        field[f] += C * 0.3;
    }

    // Compute 1-form
    Eigen::VectorXd form = Eigen::VectorXd::Zero(ps_mesh->nEdges());
    for (Edge e : gc_mesh->edges()) {
        Halfedge he = e.halfedge(); // should always be an interior halfedge
        Vector3 f1 = field[he.face()];
        Vector3 f2 = he.twin().isInterior() ? field[he.twin().face()] : Vector3{0, 0, 0};
        form[e.getIndex()] = 0.5 * dot(f1 + f2, geometry->halfedgeVector(he));
    }

    omega = form;
    omega_initialized = true;

    auto vis = ps_mesh->addOneFormTangentVectorQuantity("Omega", omega, polyscopeEdgeOrientations(*gc_mesh));
    vis->setEnabled(true);
    vis->setVectorLengthScale(VECTOR_SCALE_DEFAULT);
    
}