#include "include/ex3.h"


void Exercise3::generate_curve_simple() {
    std::cout << "Generating simple curve" << std::endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);

    Point center = {0., 0., 0.};
    double radius = 1;

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    for (size_t i = 0; i < num_vertices; ++i) {
        Point pt;
        double frac = static_cast<double>(i) / static_cast<double>(num_vertices);
        pt[0] = center[0] + radius * (cos(2. * M_PI * frac) + distribution(generator));
        pt[1] = center[0] + radius * (sin(2. * M_PI * frac) + distribution(generator));
        pt[2] = 0.;
        points.push_back(pt);
    }

    for (size_t i = 0; i < num_vertices; ++i) {
        edges.push_back({i, (i + 1) % num_vertices});
    }
    Curve = polyscope::registerCurveNetwork("curve", points, edges);
}

void Exercise3::generate_curve_figure_eight() {
    std::cout << "Generating figure-8" << std::endl;

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    // TODO: ADD YOUR CODE HERE

    Curve = polyscope::registerCurveNetwork("curve", points, edges);
}

void Exercise3::generate_curve_limacon() {
    std::cout << "Generating limaçon" << std::endl;

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    // TODO: ADD YOUR CODE HERE

    Curve = polyscope::registerCurveNetwork("curve", points, edges);
}

void Exercise3::generate_curve_3d() {
    std::cout << "Generating 3D curve" << std::endl;

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    // TODO: ADD YOUR CODE HERE

    Curve = polyscope::registerCurveNetwork("curve", points, edges);
}

void Exercise3::laplacian_smoothing() {

    // TODO: ADD YOUR CODE HERE

    Curve->updateNodePositions(points);
}

void Exercise3::osculating_circle() {

    // TODO: ADD YOUR CODE HERE

    Curve->updateNodePositions(points);
}

glm::vec3 Exercise3::CircumscribedCircle(const Point& a, const Point& b, const Point& c) {

    Point center = {0., 0., 0.};

    // TODO: ADD YOUR CODE HERE

    return center;
}