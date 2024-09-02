#ifndef CURVESMOOTHING_HH
#define CURVESMOOTHING_HH

#include "polyscope/curve_network.h"

using namespace polyscope;

class Exercise3 {
  public:
    using Point = glm::vec3;

    Exercise3() : Curve(nullptr) {}
    ~Exercise3() {}
    size_t num_vertices = 30;
    float epsilon = 0.01;
    int smooth_type = 1;
    int num_iter = 10;

    void apply_smooth() {
        // std::cout << "Smoothing type: " << smooth_type << std::endl;
        if (smooth_type == 0) {
            laplacian_smoothing();
        } else if (smooth_type == 1) {
            osculating_circle();
        }
    };
    void generate_curve_simple();
    void generate_curve_figure_eight();
    void generate_curve_limacon();
    void generate_curve_3d();

    void laplacian_smoothing();
    void osculating_circle();
    Point CircumscribedCircle(const Point& a, const Point& b, const Point& c);

  private:
    CurveNetwork* Curve;
    std::vector<Point> points;
    std::vector<std::array<size_t, 2>> edges;
};

#endif