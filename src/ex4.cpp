#include "include/ex4.h"
#include "include/utils.h"

// ========================================================================
// PROBLEM 1.1
// ========================================================================
void Exercise4::compute_normals_with_constant_weights() {
    // ------------- IMPLEMENT HERE ---------
    // TODO: Compute the normals for each vertex v in the mesh using the constant 
    // weights technique (see handout) and store inside vertex_normals
    // ------------- IMPLEMENT HERE ---------

}

// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Exercise4::compute_normals_by_area_weights() {
    // ------------- IMPLEMENT HERE ---------
    // TODO: Compute the normals for each vertex v in the mesh using weights 
    // proportional to the dual area (see handout) 
    // and store inside vertex_normals
    // ------------- IMPLEMENT HERE ---------

}

// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Exercise4::compute_normals_with_angle_weights() {
    // ------------- IMPLEMENT HERE ---------
    // TODO: Compute the normals for each vertex v in the mesh using weights 
    // proportional to incident angles (see handout) 
    // and store inside vertex_normals
    // ------------- IMPLEMENT HERE ---------

}

// ========================================================================
// EXERCISE 2.1
// ========================================================================
void Exercise4::calc_uniform_laplacian() {
    // // ------------- IMPLEMENT HERE ---------
    // // For each non-boundary vertex, compute uniform Laplacian operator vector
    // // and store the dot product in vertex_curvature_weights.
    // // Store min and max values in min_curvature and max_curvature.

    // // ------------- IMPLEMENT HERE ---------


    std::cout << "Min. uniform value is: " << min_curvature << std::endl;
    std::cout << "Max. uniform value is: " << max_curvature << std::endl;

}

// ========================================================================
// EXERCISE 2.2
// ========================================================================
void Exercise4::calc_mean_curvature() {
    // ------------- IMPLEMENT HERE ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the dot product of the Laplace-Beltrami approximation.
    // ------------- IMPLEMENT HERE ---------

    std::cout << "Min. Laplace-Beltrami value is: " << min_curvature << std::endl;
    std::cout << "Max. Laplace-Beltrami value is: " << max_curvature << std::endl;
}

// ========================================================================
// EXERCISE 2.3
// ========================================================================
void Exercise4::calc_gauss_curvature() {
    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property vertex_curvature_weights.
    // Hint: When calculating angles out of dot products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.

    // ------------- IMPLEMENT HERE ---------

    std::cout << "Min. Gaussian value is: " << min_curvature << std::endl;
    std::cout << "Max. Gaussian value is: " << max_curvature << std::endl;
}

//====================================================================================================================//

void Exercise4::show_normal() {
    switch (normal_type) {
    case 0:
        compute_normals_with_constant_weights();
        break;
    case 1:
        compute_normals_by_area_weights();
        break;
    case 2:
        compute_normals_with_angle_weights();
        break;
    }
    auto vvq = ps_mesh->addVertexVectorQuantity("normal", vertex_normals);
    vvq->setEnabled(true);
}

void Exercise4::show_curvature() {
    switch (curvature_type) {
    case 0:
        calc_uniform_laplacian();
        break;
    case 1:
        calc_mean_curvature();
        break;
    case 2:
        calc_gauss_curvature();
        break;
    }
    color_coding_d(vertex_curvature_weights);
}

void Exercise4::color_coding_d(const std::vector<double> data, const double _min_value,
                                         const double _max_value, const int _bound) {
    std::vector<double> datavalue(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        datavalue[i] = data[i];
    }
    auto min_value = _min_value;
    auto max_value = _max_value;

    if (min_value == 0 && max_value == 0) {
        // discard upper and lower bound
        auto n = data.size() - 1;
        auto i = n / _bound;

        std::sort(datavalue.begin(), datavalue.end());
        min_value = datavalue[i];
        max_value = datavalue[n - 1 - i];
    }

    const auto range = max_value - min_value;
    std::vector<std::array<double, 3>> psColor;

    for (auto vh : gc_mesh->vertices()) {
        auto t = ((data[vh.getIndex()]) - min_value) / range;
        auto color = util::map_val2color(t, 0, 1.0);
        psColor.push_back({color.x, color.y, color.z});
    }
    auto vcq = ps_mesh->addVertexColorQuantity("Vertex Colors", psColor);
    vcq->setEnabled(true);
}
