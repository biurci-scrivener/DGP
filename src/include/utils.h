#ifndef UTILS_H
#define UTILS_H

#include "geometrycentral/utilities/vector3.h"
#include <glm/glm.hpp>

namespace util {
    using Point = geometrycentral::Vector3;
    glm::vec4 map_val2color(float iso, float val0, float val4);
    double cal_edge_length(const Point &a, const Point &b);
    Point circumscribedCircle(const Point &a, const Point &b, const Point &c);
}

#endif // UTILS_H
