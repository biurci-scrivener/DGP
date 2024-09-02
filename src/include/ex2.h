#ifndef ISOCONTOURING_HH
#define ISOCONTOURING_HH

#include <geometrycentral/surface/meshio.h>
#include <memory>
#include <polyscope/curve_network.h>
#include <polyscope/surface_mesh.h>
#include "util.h"
#include "FileSystem.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Exercise2 {
  public:
    using Point = geometrycentral::Vector3;

    Exercise2() : geometry(nullptr), gc_mesh(nullptr), ps_mesh(nullptr) {}
    ~Exercise2() {
        vertices.clear();
        edges.clear();
        segment_points.clear();
        geometry.reset();
        gc_mesh.reset();
        ps_mesh = nullptr;
        // std::cout << "Exercise2 destructor" << std::endl;
        // std::cout <<"geometry.usecount  "<< geometry.use_count() << std::endl;
        // std::cout <<"gc_mesh.usecount  "<< gc_mesh.use_count() << std::endl;
    }

    void preparemesh(PolygonInstance& polygon_instance) {
        gc_mesh = polygon_instance.mesh;
        geometry = polygon_instance.geometry;
        ps_mesh = polygon_instance.ps_mesh;
    }

    void show_isovalue_and_level_set();

    int function_id = 0;

  private:
    double iso_value(const Point& _pt) const;

    void show_isovalue();

    void compute_segment_points();

    void create_level_set0_segments();

    void create_segment(const std::vector<Point>& _points);


  private:
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<size_t, 2>> edges;
    std::shared_ptr<VertexPositionGeometry> geometry;
    std::shared_ptr<surface::SurfaceMesh> gc_mesh;
    polyscope::SurfaceMesh* ps_mesh;
    std::vector<Point> segment_points;
};

#endif
