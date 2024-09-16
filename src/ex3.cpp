#include "include/ex3.h"


void Exercise3::generate_curve_simple() {
    std::cout << "Generating simple curve" << std::endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

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
    // generate_world_axes(); // uncomment to visualize world axes
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

void Exercise3::generate_world_axes() {
    PCPointCloud * ax_x = new PCPointCloud(1);
    PointPositionGeometry * geo_x = new PointPositionGeometry(*ax_x);
    geo_x->positions[0] = {0.,0.,-1.};
    PointData<Point> ax_x_q(*ax_x);
    ax_x_q[0] = {0.25,0.,0.};
    auto pc_x = polyscope::registerPointCloud("World x-axis", geo_x->positions);
    pc_x->setPointColor({0.,0.,0.});
    auto vis_x = pc_x->addVectorQuantity("Axes", ax_x_q, VectorType::AMBIENT);
    vis_x->setVectorColor({1.,0.,0.});
    vis_x->setVectorRadius(0.01);
    vis_x->setEnabled(true);

    PCPointCloud * ax_y = new PCPointCloud(1);
    PointPositionGeometry * geo_y = new PointPositionGeometry(*ax_y);
    geo_y->positions[0] = {0.,0.,-1.};
    PointData<Point> ax_y_q(*ax_y);
    ax_y_q[0] = {0.,0.25,0.};
    auto pc_y = polyscope::registerPointCloud("World y-axis", geo_y->positions);
    pc_y->setPointColor({0.,0.,0.});
    auto vis_y = pc_y->addVectorQuantity("Axes", ax_y_q, VectorType::AMBIENT);
    vis_y->setVectorColor({0.,1.,0.});
    vis_y->setVectorRadius(0.01);
    vis_y->setEnabled(true);

    PCPointCloud * ax_z = new PCPointCloud(1);
    PointPositionGeometry * geo_z = new PointPositionGeometry(*ax_z);
    geo_z->positions[0] = {0.,0.,-1.};
    PointData<Point> ax_z_q(*ax_z);
    ax_z_q[0] = {0.,0.,0.25};
    auto pc_z = polyscope::registerPointCloud("World z-axis", geo_z->positions);
    pc_z->setPointColor({0.,0.,0.});
    auto vis_z = pc_z->addVectorQuantity("Axes", ax_z_q, VectorType::AMBIENT);
    vis_z->setVectorColor({0.,0.,1.});
    vis_z->setVectorRadius(0.01);
    vis_z->setEnabled(true);
}