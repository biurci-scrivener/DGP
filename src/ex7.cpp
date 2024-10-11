#include "include/ex7.h"
#include <iostream>
#include "geometrycentral/surface/halfedge_element_types.h"

// ======================================================================
// Part 1.1 Mapping to circle
// ========================================================================
void Exercise7::map_suface_boundary_to_circle() {

    // initialize all points to origin

    for (Vertex v: tex_mesh->vertices()) {
        tex_geometry->vertexPositions[v] = _origin;
    }

    // TODO: ADD YOUR CODE HERE

    // register mesh with polyscope

    tex_geometry->refreshQuantities();
    ps_tex_mesh = polyscope::registerSurfaceMesh("tex_coords_circle", tex_geometry->vertexPositions, tex_mesh->getFaceVertexList());
    ps_tex_mesh->setEdgeWidth(1.);
    ps_tex_mesh->setSurfaceColor({1.,1.,1.});
    is_mapped = true;
}

// ======================================================================
// Part 1.2 Interactively smoothing the texture
// ========================================================================
void Exercise7::explicitly_smooth_texture() {

    // TODO: ADD YOUR CODE HERE

    geometry->refreshQuantities();
    auto qParam = ps_mesh->addVertexParameterizationQuantity("tex_coords", _tex_coords);
    qParam->setCheckerSize(1.);
    qParam->setEnabled(true);
    ps_tex_mesh->updateVertexPositions(tex_geometry->vertexPositions);

}

// ======================================================================
// Part 1.3 Implicitly smoothing the texture
// ========================================================================
void Exercise7::implicitly_smooth_texture() {
    
    // TODO: ADD YOUR CODE HERE

    geometry->refreshQuantities();
    auto qParam = ps_mesh->addVertexParameterizationQuantity("tex_coords", _tex_coords);
    qParam->setCheckerSize(1.);
    qParam->setEnabled(true);
    ps_tex_mesh->updateVertexPositions(tex_geometry->vertexPositions);

}

// ======================================================================
// Part 2 Minimal Surfaces
// ======================================================================
void Exercise7::minimal_surface() {

    // TODO: ADD YOUR CODE HERE

    geometry->refreshQuantities();
    ps_mesh->updateVertexPositions(geometry->vertexPositions);

}