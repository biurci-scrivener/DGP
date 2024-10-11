#include "include/ex5.h"

/*

    PART 1 below. Implement the following functions:
    - compare_dense_solvers
    - analyze_sparsity_pattern
    - compare_sparse_solvers

*/

void Exercise5::compare_dense_solvers() {
    
    // TODO: Problem 1.1 is here.

}

void Exercise5::analyze_sparsity_pattern() {

    // TODO: Problem 1.2.1 is here.

}

void Exercise5::compare_sparse_solvers() {

    // TODO: Problem 1.2.2 is here.

}

/*

    PART 2 below. Implement the following functions:
    - deflate_matrix
    - find_max_eigenvec_eigenval
    - find_all_eigenvectors_eigenvalues
    - find_all_eigenvectors_eigenvalues_sparse

*/

Eigen::Matrix3d deflate_matrix(Eigen::Matrix3d A, Eigen::Vector3d eigenvec, double eigenval) {

    // TODO: Implement matrix deflation.

}

std::tuple<Eigen::Vector3d, double, std::vector<Eigen::Vector3d>> find_max_eigenvec_eigenval(Eigen::Matrix3d A) {

    Eigen::Vector3d v_current = Eigen::MatrixXd::Random(3,1);
    double w_current = 0;
    std::vector<Eigen::Vector3d> first_iterates(5, {0.,0.,0.});

    // TODO: Implement the power method.
    // Store first 5 iterates in first_iterates, such that
    // each iterate is a unit eigenvector (estimate) scaled by its corresponding eigenvalue (estimate).

    return std::make_tuple(v_current, w_current, first_iterates);

}

void Exercise5::find_all_eigenvectors_eigenvalues() {

    Eigen::Matrix3d A;
    A <<    0.8266,     -0.3370,    0.0899,
            -0.3370,    1.0687,     -0.3629,
            0.0899,     -0.3629,    0.8547;

    Eigen::Vector3d v_1;
    Eigen::Vector3d v_2;
    Eigen::Vector3d v_3;
    double lambda_1;
    double lambda_2;
    double lambda_3;
    std::vector<Eigen::Vector3d> v1_iterates;
    std::vector<Eigen::Vector3d> v2_iterates;
    std::vector<Eigen::Vector3d> v3_iterates;

    // TODO: Problem 2.1 is here.
    
    p2_generate_iterate_plot(v1_iterates, v2_iterates, v3_iterates);
    generate_ellipse("Final", v_1, v_2, v_3, lambda_1, lambda_2, lambda_3);
    generate_ellipsoid_axes(v_1 * lambda_1, v_2 * lambda_2, v_3 * lambda_3);
    generate_world_axes();

}

/*

    PART 2.2 below. Implement the following functions:
    - inverse_iteration
    - find_closest_eigenvalues

*/

std::tuple<std::vector<double>, std::vector<Eigen::VectorXd>> inverse_iteration(Exercise5::SMatd &matrix, double target_eigval) {
    
    std::vector<double> iters_val;
    std::vector<Eigen::VectorXd> iters_vec;

    // TODO: Implement inverse power method with shift.
    // Store iterates in iters_val, iters_vec.
    
    return std::make_tuple(iters_val, iters_vec);
}

void Exercise5::find_closest_eigenvalues() {

    // TODO: Problem 2.2 is here.

}

/*

    PART 3 below. Implement the following functions:
    - func
    - grad
    - hess
    - gradient_descent
    - newtons_method

*/

double func(double x, double y, int id) {

    // TODO: define f and g.
    // if id = 0, the function being queried is f.
    // if id = 1, the function being queried is g.

    switch (id) {
        case 0:
            return 0;
        case 1:
            return 0;
        default:
            return 0;
    }
}

Eigen::Vector2d grad(double x, double y, int id) {

    // TODO: compute gradients for f and g.
    // if id = 0, the function being queried is f.
    // if id = 1, the function being queried is g.

    double dx = 0;
    double dy = 0;

    switch (id) { 
        case 0: 
            dx = 0;
            dy = 0;
            break;
        case 1: 
            dx = 0;
            dy = 0;
            break;
    }
    return {dx, dy};
}

Eigen::Matrix2d hess(double x, double y, int id) {

    // TODO: compute Hessian for f and g.
    // if id = 0, the function being queried is f.
    // if id = 1, the function being queried is g.

    Eigen::Matrix2d h;

    double a = 0;
    double b = 0;
    double c = 0;
    double d = 0;

    switch(id) {
        case 0:
            a = 0;
            b = 0;
            c = 0;
            d = 0;
            break;
        case 1:
            a = 0;
            b = 0;
            c = 0;
            d = 0;
            break;
    }

    h << a, b, c, d;
    return h;
}

void Exercise5::gradient_descent(int func_id, double step_size) {

    int MAX_ITERS = 5000;
    std::vector<Eigen::Vector3d> iters;
    
    Eigen::Vector2d curr = Eigen::MatrixXd::Random(2,1) * 5;
    double f_curr = func(curr[0], curr[1], func_id);

    // TODO: Problem 3.1 is here.

    p3_generate_iterate_plot(iters, func_id);

}

void Exercise5::newtons_method(int func_id) {

    // TODO: Problem 3.2 is here.

    int MAX_ITERS = 5000;
    std::vector<Eigen::Vector3d> iters;
    
    Eigen::Vector2d curr = Eigen::MatrixXd::Random(2,1) * 5;
    double f_curr = func(curr[0], curr[1], func_id);

    // TODO: Problem 3.1 is here.

    p3_generate_iterate_plot(iters, func_id);
    
}

// === UTILITY FUNCTIONS BELOW: do not change ===

void Exercise5::polyscope_plot_vector(std::string name, Point pos, Eigen::Vector3d dir, Point color) {
    std::vector<Point> p;
    p.push_back(pos);
    auto _pc = polyscope::registerPointCloud(name, p);
    _pc->setPointColor({0.,0.,0.});
    _pc->setPointRadius(0.005);
    std::vector<Eigen::Vector3d> q;
    q.push_back(dir);
    auto _vec = _pc->addVectorQuantity(std::to_string(1), q, VectorType::AMBIENT);
    _vec->setVectorColor(color);
    _vec->setVectorRadius(0.01);
    _vec->setEnabled(true);
}

void Exercise5::polyscope_plot_vectors(std::string name, Point pos, std::vector<Eigen::Vector3d> dirs, Point base_color) {
    std::vector<Point> p;
    p.push_back(pos);
    Point blend;
    auto _pc = polyscope::registerPointCloud(name, p);
    _pc->setPointColor({0.,0.,0.});
    _pc->setPointRadius(0.005);
    if (((base_color.x + base_color.y + base_color.z) / 3.) <= 0.5 ) {
        blend = {0.8,0.8,0.8};
    } else {
        blend = {0.,0.,0.};
    }
    for (size_t i = 0; i < dirs.size(); i++) {
        std::vector<Eigen::Vector3d> q;
        q.push_back(dirs[i]);
        auto _vec = _pc->addVectorQuantity(std::to_string(i), q, VectorType::AMBIENT);
        _vec->setVectorColor((((float)i + 1) / dirs.size()) * base_color + (1 - (((float)i + 1) / dirs.size())) * blend);
        _vec->setVectorRadius(0.01);
        _vec->setEnabled(true);
    }
}

void Exercise5::polyscope_plot_vectors(std::string name, Point pos, std::vector<Eigen::Vector3d> dirs, std::vector<Point> colors) {
    std::vector<Point> p;
    p.push_back(pos);
    auto _pc = polyscope::registerPointCloud(name, p);
    _pc->setPointColor({0.,0.,0.});
    _pc->setPointRadius(0.005);
    for (size_t i = 0; i < dirs.size(); i++) {
        std::vector<Eigen::Vector3d> q;
        q.push_back(dirs[i]);
        auto _vec = _pc->addVectorQuantity(std::to_string(i), q, VectorType::AMBIENT);
        _vec->setVectorColor(colors[i]);
        _vec->setVectorRadius(0.01);
        _vec->setEnabled(true);
    }
}

void Exercise5::polyscope_plot_vectors(std::string name, std::vector<Vector3> pos, std::vector<Eigen::Vector3d> dirs, Point color) {
    auto _pc = polyscope::registerPointCloud(name, pos);
    _pc->setPointColor({0.,0.,0.});
    _pc->setPointRadius(0.005);
    auto _vec = _pc->addVectorQuantity("GD Iters", dirs, VectorType::AMBIENT);
    _vec->setVectorColor(color);
    _vec->setVectorRadius(0.01);
    _vec->setEnabled(true);
}

void Exercise5::generate_world_axes() {
    std::vector<Point> colors;
    colors.push_back({0.3,0.,0.});
    colors.push_back({0.,0.3,0.});
    colors.push_back({0.,0.,0.3});
    std::vector<Eigen::Vector3d> dirs;
    dirs.push_back({0.25,0.,0.});
    dirs.push_back({0.,0.25,0.});
    dirs.push_back({0.,0.,0.25});
    polyscope_plot_vectors("World axes", {1.,0.,0.}, dirs, colors);
}

void Exercise5::generate_ellipsoid_axes(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c) {
    std::vector<Point> colors;
    colors.push_back({1.,0.,0.});
    colors.push_back({0.,1.,0.});
    colors.push_back({0.,0.,1.});
    std::vector<Eigen::Vector3d> dirs;
    dirs.push_back(a);
    dirs.push_back(b);
    dirs.push_back(c);
    polyscope_plot_vectors("Ellipsoid axes", {0.,0.,0.}, dirs, colors);
}

void Exercise5::p2_generate_iterate_plot(std::vector<Eigen::Vector3d> v1_it, std::vector<Eigen::Vector3d> v2_it, std::vector<Eigen::Vector3d> v3_it) {
    polyscope_plot_vectors("v1 iterates", {0.,0.,0.}, v1_it, {1.,0.1,0.1});
    polyscope_plot_vectors("v2 iterates", {0.,0.,0.}, v2_it, {0.1,1.,0.1});
    polyscope_plot_vectors("v3 iterates", {0.,0.,0.}, v3_it, {0.1,0.1,1.});
}

void Exercise5::generate_sparsity_plot(SMatd A) {
    int width = A.cols();
    int height = A.rows();
    std::vector<Point> imageColor(width * height, {0.,0.,0.});
    for (int k = 0; k < A.outerSize(); k++) {
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            imageColor[it.row() * width + it.col()] = {1.,1.,1.};
        }
    }
    
    auto im = polyscope::addColorImageQuantity("Sparsity Plot", width, height, imageColor, polyscope::ImageOrigin::UpperLeft);
    im->setEnabled(true);
    
}

void Exercise5::plot_function(int func_id) {

    polyscope::view::setUpDir(UpDir::ZUp);

    std::vector<Vector3> points_new;

    double PLOTTING_SCALE = (func_id == 0) ? PLOTTING_SCALE_f : PLOTTING_SCALE_g;

    for (Vertex v: gc_mesh->vertices()) {
        Vector3 pos = geometry->vertexPositions[v];
        double f_val = func(pos.x, pos.y, func_id);
        points_new.push_back({pos.x, pos.y, f_val * PLOTTING_SCALE});
    }
    
    for (size_t j = 0; j < gc_mesh->nVertices(); j++) {
        geometry->vertexPositions[j] = points_new[j];
    }

    ps_mesh->updateVertexPositions(geometry->vertexPositions);
    ps_mesh->setTransparency(0.33);

}

void Exercise5::p3_generate_iterate_plot(std::vector<Eigen::Vector3d> iterates, int func_id) {

    std::vector<Eigen::Vector3d> dirs;
    std::vector<Vector3> pos;

    if (iterates.size() < 2) {
        return;
    }

    double PLOTTING_SCALE = (func_id == 0) ? PLOTTING_SCALE_f : PLOTTING_SCALE_g;

    for (size_t i = 0; i < iterates.size() - 1; i++) {
        Eigen::Vector3d new_dir = iterates[i + 1] - iterates[i];
        new_dir[2] *= PLOTTING_SCALE;
        Vector3 new_pos = {iterates[i][0], iterates[i][1], iterates[i][2] * PLOTTING_SCALE};
        dirs.push_back(new_dir);
        pos.push_back(new_pos);
    }

    polyscope_plot_vectors("Optimization iterates", pos, dirs, {1.,0.,0.});

    polyscope::state::lengthScale = 28;
    polyscope::state::boundingBox = std::make_tuple(Point{-5, -5, 0}, Point{5, 5, 30});

}

void Exercise5::generate_ellipse(std::string name, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, double a_len, double b_len, double c_len) {

    // if ((abs(a.dot(b)) > 1e-3) || (abs(a.dot(c)) > 1e-3) || (abs(b.dot(c)) > 1e-3)) {
    //     std::cout << "ERROR: axes inputs to create_ellipse_with_axes not pairwise orthogonal." << std::endl;
    // }

    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<int,3>> idxs;

    Eigen::Matrix3d A;
    A.col(0) << a;
    A.col(1) << b;
    A.col(2) << c;

    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
            double phi = (i/31.) * M_PI;
            double theta = (j/32.) * 2 * M_PI;
            Eigen::Vector3d this_pt = {
                a_len * sin(theta) * cos(phi),
                b_len * sin(theta) * sin(phi),
                c_len * cos(theta)
            };
            vertices.push_back(A * this_pt);
        }
    }

    for (int i = 0; i < 31; i++) {
        for (int j = 0; j < 32; j++) {
            int next_stack = (i + 1);
            int next_slice = (j + 1) % 32;
            idxs.push_back({
                i * 32 + j,
                next_stack * 32 + j,
                i * 32 + next_slice
            });
            idxs.push_back({
                i * 32 + next_slice,
                next_stack * 32 + j,
                next_stack * 32 + next_slice
            });
        }
    }

    auto ellipse = polyscope::registerSurfaceMesh(name, vertices, idxs);
    ellipse->setTransparency(0.25);

}