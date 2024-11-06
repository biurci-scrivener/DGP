#ifndef FILESYSTEM_HH
#define FILESYSTEM_HH


// A new class to manage the file system
// It will load the obj file and store the mesh and geometry




#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
// #include "polyscope/messages.h"
#include "polyscope/messages.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace polyscope;

class PolygonInstance {
   public:
    PolygonInstance(std::string name, int index) : name(name), index(index) {
        try {
            std::tie(mesh, geometry) = readSurfaceMesh(name);
            ps_mesh = polyscope::registerSurfaceMesh(
                name, geometry->vertexPositions, mesh->getFaceVertexList());
            // fixes bug in https://github.com/nmwsharp/polyscope/issues/235
            ps_mesh->setEdgePermutation(polyscopePermutations(*mesh)[2].first, mesh->nEdges());
        } catch (const std::runtime_error& e) {
            polyscope::warning("Error loading mesh: {}", e.what());
            polyscope::error("Error loading mesh: {}");
            std::cerr << "Error loading mesh: " << e.what() << std::endl;
        }
    }
    ~PolygonInstance() {
        //  std::cout << "PolygonInstance destructor" << std::endl;
        // std::cout <<"geometry.usecount  "<< geometry.use_count() <<
        // std::endl; std::cout <<"trimesh.usecount  "<< mesh.use_count() <<
        // std::endl;
    }
    void clear() { polyscope::removeStructure(ps_mesh); }

    std::string name;
    int index;
    std::shared_ptr<geometrycentral::surface::SurfaceMesh> mesh;
    std::shared_ptr<geometrycentral::surface::VertexPositionGeometry> geometry;
    polyscope::SurfaceMesh* ps_mesh;
};

class FileSystem {
   public:
    FileSystem() = default;
    ~FileSystem() { clearAll(); }

    void loadObj(const std::string& filePath) {
        std::string fileName =
            filePath.substr(filePath.find_last_of("/\\") + 1);
        fileName = fileName.substr(0, fileName.find_last_of('.'));
        try {
            PolygonInstance* PI = new PolygonInstance(filePath, objnum);
            objlst.emplace_back(PI);
            objnamelist.push_back(fileName);
            setcurpoly(objnum);
            objnum++;
        } catch (const std::runtime_error& e) {
            // std::cerr << "Error loading object: " << e.what() << std::endl;
            return;  
        }
    }

    void clearAll() {
        objlst.clear();
        objnamelist.clear();
        objnum = 0;
    }

    const std::vector<std::unique_ptr<PolygonInstance>>* getObjList() const {
        return &objlst;
    }
    int objnum = 0;
    std::vector<std::string> objnamelist;
    void setcurpoly(int id) { cur_polygon = objlst[id].get(); }

    PolygonInstance* getCurrentPolygon() const { return cur_polygon; }

   private:
    PolygonInstance* cur_polygon = nullptr;
    std::vector<std::unique_ptr<PolygonInstance>> objlst;
};

#endif