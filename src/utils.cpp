#include "include/utils.h"

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