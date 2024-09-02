#include "include/utils.h"


double util::cal_edge_length(const util::Point &a, const util::Point &b) {
  return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2));
}

glm::vec4 util::map_val2color(float iso, float val0, float val4) {

  float val2 = (val0 + val4) / 2;
  float val1 = (val0 + val2) / 2;
  float val3 = (val2 + val4) / 2;

  if (val4 <= val0) {
    return glm::vec4(0, 0, 1, 1);
  }

  if (iso < val0)
    return glm::vec4(0, 0, 1, 1);
  if (iso > val4)
    return glm::vec4(1, 0, 0, 1);

  if (iso <= val2) {
    if (iso <= val1) {
      float u = float(1.0 * (iso - val0) / (val1 - val0));
      return glm::vec4(0, u, 1, 1);
    } else {
      float u = float(1.0 * (iso - val1) / (val2 - val1));
      return glm::vec4(0, 1, 1 - u, 1);
    }

  } else {
    if (iso <= val3) {
      float u = float(1.0 * (iso - val2) / (val3 - val2));
      return glm::vec4(u, 1, 0, 1);
    }

    else {
      float u = float(1.0 * (iso - val3) / (val4 - val3));
      return glm::vec4(1, 1 - u, 0, 1);
    }
  }
}

util::Point util::circumscribedCircle(const Point &a, const Point &b,
                               const Point &c) {
  double ax = a[0];
  double ay = a[1];

  double bx = b[0];
  double by = b[1];

  double cx = c[0];
  double cy = c[1];

  double D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

  double Ux =
      ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) +
       (cx * cx + cy * cy) * (ay - by)) /
      D;
  double Uy =
      ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) +
       (cx * cx + cy * cy) * (bx - ax)) /
      D;

  Point center = {Ux, Uy, 0};
  return center;
}