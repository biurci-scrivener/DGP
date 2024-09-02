#include "include/ex2.h"
#include "include/utils.h"

bool near_zero(double x) {
    return std::abs(x) < 1e-9;
}

double Exercise2::iso_value(const Point& _pt) const {

    double iso = 0.;

    // TODO: ADD YOUR CODE HERE

    return iso;
}

void Exercise2::compute_segment_points() {
    segment_points.clear();

    // TODO: ADD YOUR CODE HERE

}

void Exercise2::show_isovalue_and_level_set() {
    show_isovalue();
    create_level_set0_segments();
}

void Exercise2::show_isovalue() {

    std::vector<double> iso_values;
    for (auto v : gc_mesh->vertices()) {
        Point point = geometry->vertexPositions[v];
        auto iv = iso_value(point);
        iso_values.push_back(iv);
    }

    auto max_iv = *std::max_element(iso_values.begin(), iso_values.end());
    auto min_iv = *std::min_element(iso_values.begin(), iso_values.end());

    const double range = max_iv - min_iv;
    std::vector<std::array<double, 3>> psColor;

    for (auto vh : gc_mesh->vertices()) {
        auto t = (iso_values[vh.getIndex()] - min_iv) / range;
        auto color = util::map_val2color(t, 0, 1.0);
        psColor.push_back({color.x, color.y, color.z});
    }

    auto cq = ps_mesh->addVertexColorQuantity("Vertex Colors", psColor);
    cq->setEnabled(true);
}

void Exercise2::create_level_set0_segments() {
    compute_segment_points();
    edges.clear();
    vertices.clear();

    if (segment_points.empty()) {
        std::cerr << "segment_points is empty" << std::endl;
        return;
    }

    if (segment_points.size() % 2 != 0) {
        std::cerr << "Size of segment_points is not even" << std::endl;
        return;
    }

    for (size_t i = 0; i < segment_points.size(); i += 2) {
        auto p0 = segment_points[i];
        auto p1 = segment_points[i + 1];
        vertices.push_back({p0[0], p0[1], p0[2]});
        vertices.push_back({p1[0], p1[1], p1[2]});
        edges.push_back({i, i + 1});
    }

    polyscope::registerCurveNetwork("Level Set 0", vertices, edges);
}