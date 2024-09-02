#include "include/Panel.h"

#include "imgui.h"

// For deformation picking interface
std::vector<int> picked_verts_deform;
std::vector<int> picked_verts_fixed;
std::vector<Point> picked_pts_deform;
std::vector<Point> picked_pts_fixed;
Point displacement_vector = Point{0., 0., 0.};

static bool Panel::VectorOfStringGetter(void* vec, int idx,
                                        const char** out_text) {
    const auto& vector = *static_cast<std::vector<std::string>*>(vec);
    if (idx < 0 || idx >= static_cast<int>(vector.size())) return false;
    *out_text = vector[idx].c_str();
    return true;
}

void Panel::DGPPanel() {
    ImGuiStyle& style = GetStyle();

    style.WindowRounding = 5.0f;  
    style.FrameRounding = 3.0f;  
    style.PopupRounding = 3.0f; 
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.1f, 0.1f, 0.1f, 1.0f);

    float margin = 10;
    float right_w = 500;

    float w = GetIO().DisplaySize.x;
    float h = GetIO().DisplaySize.y;
    static bool firstTime = true;

    if (firstTime) {
        SetNextWindowPos(ImVec2(w - right_w - margin, margin));
        SetNextWindowSize(ImVec2(right_w, h * 0.5));
        firstTime = false;
    }

    static bool open = ImGuiTreeNodeFlags_DefaultOpen;

    ImVec2 full_block_size = ImVec2(480, 25);
    ImVec2 half_button_size = ImVec2(240, 25);
    double spacing = 10;
    double second_offset = half_button_size.x * 0.6 + spacing;

    Begin("Digital Geometry Processing", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    // PushItemWidth(half_button_size.x);

    if (ImGui::Button("Import mesh", half_button_size)) {
        IGFD::FileDialogConfig config;
        config.path = ".";
        ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey",
                                                "Choose File", ".*", config);
    }

    if (ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey")) {
        if (ImGuiFileDialog::Instance()->IsOk()) {  // action if OK
            std::string fullPath =
                ImGuiFileDialog::Instance()->GetFilePathName();
            // std::string fileName = extractBeforeLastDot(
            //     ImGuiFileDialog::Instance()->GetCurrentFileName());
            fs.loadObj(fullPath);
        }
        ImGuiFileDialog::Instance()->Close();
    }

    SameLine(250);

    if (ImGui::Button("Clear all meshes", half_button_size)) {
        fs.clearAll();
        ex7.is_mapped = false;
        polyscope::removeAllStructures();
    }

    Text("Select target mesh");
    SameLine(second_offset);
    static int cur_mesh_id = 0;
    if (Combo("##mesh", &cur_mesh_id, VectorOfStringGetter, &fs.objnamelist,
              fs.objnum)) {
        fs.setcurpoly(cur_mesh_id);
    }

    if (CollapsingHeader("Exercise 1", open)) {
        if (Button("Solve sparse linear system", full_block_size)) {
            ex1.solve_sparse_linear_system();
        };

        static int current_item = 0;

        std::vector<const char*> fuctions = {
            "Function 1: Height, Vertices",
            "Function 2: Length, Edges",
            "Function 3: Normal -> Color, Faces",
        };
        AlignTextToFramePadding();

        Text("Current Plot");
        SameLine(second_offset);
        Combo("##", &current_item, fuctions.data(), fuctions.size());

        if (Button("Generate plot", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                ex1.set_mesh(*cur_polygon);
                ex1.generate_plot(current_item);
            }
        }
    };

    if (CollapsingHeader("Exercise 2", open)) {
        static int current_item = 0;

        std::vector<const char*> fuctions = {
            "sqrt(x^2+y^2)-1=0",
            "y^2-sin(x^2)-1=0",
            "sin(2x+2y)-cos(4xy)=0",
            "(3x^2 - y^2)^2 * y^2 - (x^2 + y^2)^4"
        };
        AlignTextToFramePadding();

        Text("Function");
        SameLine(second_offset);
        Combo("##", &current_item, fuctions.data(), fuctions.size());

        if (Button("Show isovalue and level set", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                ex2.preparemesh(*cur_polygon);
                ex2.function_id = current_item;
                ex2.show_isovalue_and_level_set();
            }
        }
    };

    if (CollapsingHeader("Exercise 3", open)) {
        static int ex3_num_points = 30;
        static int curve_gen_id = 0;

        std::vector<const char*> curve_types = {
            "Simple",
            "Figure-8",
            "Limaçon",
            "3D"
        };
        AlignTextToFramePadding();

        Text("Number of points:"); 
        SameLine(second_offset);
        if (SliderInt("##numPoints", &ex3_num_points, 0, 100)) {
            ex3.num_vertices = ex3_num_points;
        };
        Combo("##Curve type", &curve_gen_id, curve_types.data(), curve_types.size());
        if (Button("Init. random curve", full_block_size)) {
            switch (curve_gen_id) {
                case 0:
                    ex3.generate_curve_simple();
                    break;
                case 1:
                    ex3.generate_curve_figure_eight();
                    break;
                case 2:
                    ex3.generate_curve_limacon();
                    break;
                case 3:
                    ex3.generate_curve_3d();
                    break;
                default:
                    break;
            }
        }

        static int ex3_iternum = 1;
        static float ex3_eps = ex3.epsilon;
        static int cur_smooth_id = 0;
        static const char* smooth_types = "Laplacian\0Osculating circle\0";

        Text("Iteration");
        SameLine(second_offset);
        if (SliderInt("##ex3_iternum ", &ex3_iternum, 1, 100)) {
            ex3.num_iter = ex3_iternum;
        };

        Text("Epsilon");
        SameLine(second_offset);
        if (SliderFloat("##ex3_eps", &ex3_eps, 1e-5, 1e-2, "%.5f")) {
            ex3.epsilon = ex3_eps;
        };

        Text("Smooth type");
        SameLine(second_offset);
        if (Combo("##smooth_types", &cur_smooth_id, smooth_types, 2)) {
            ex3.smooth_type = cur_smooth_id;
        };

        if (Button("Apply smoothing", full_block_size)) {
            ex3.apply_smooth();
        };
    };

    if (CollapsingHeader("Exercise 4", open)) {
        AlignTextToFramePadding();
        static int ex4_normal_type = 0;
        const char* normal_types = "Uniform\0Area-weighted\0Angle-weighted\0";
        static int ex4_curvature_type = 0;
        const char* curvature_types =
            "Uniform\0Mean curvature\0Gaussian curvature\0";

        Text("Normal type");
        SameLine(second_offset);
        if (Combo("##ex4_normal", &ex4_normal_type, normal_types, 3)) {
            ex4.normal_type = ex4_normal_type;
        };
        if (Button("Compute Normals", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");

            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                ex4.set_mesh(*cur_polygon);
                ex4.show_normal();
            }
        };

        Text("Curvature type");
        SameLine(second_offset);
        if (Combo("##ex4_curvature", &ex4_curvature_type, curvature_types, 3)) {
            ex4.curvature_type = ex4_curvature_type;
        };
        if (Button("Compute curvature", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");

            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                ex4.set_mesh(*cur_polygon);
                ex4.show_curvature();
            }
        };
    };

    if (CollapsingHeader("Exercise 6", open)) {
        AlignTextToFramePadding();
        static int ex6_cur_laplacian = 0;
        const char* laplacian_types = "Uniform\0Cotangent\0Implicit\0";
        static int ex6_iteration = 10;
        static float ex6_timestep = 1e-5;
        static float ex6_feature_enhancement_coefficient = 2.0;


        SetWindowFontScale(1.5f);
        Text("Smoothing");
        SetWindowFontScale(1.0f);

        Text("Operator type");
        SameLine(second_offset);
        if (Combo("##smoothtype", &ex6_cur_laplacian, laplacian_types, 3)) {
            ex6.smooth_type = ex6_cur_laplacian;
        }

        Text("Iterations (explicit)");
        SameLine(second_offset);
        if (SliderInt("##iteration", &ex6_iteration, 1., 1000.)) {
            ex6.num_iters = ex6_iteration;
        };

        Text("Timestep (implicit)");
        SameLine(second_offset);
        if (InputFloat("##Timestep", &ex6_timestep, 1e-5, 1e-1, "%.5f")) {
            ex6.timestep = ex6_timestep;
        };

        Text("F.E. coefficient");
        SameLine(second_offset);
        if (SliderFloat("##Feature Enhancement Coefficient",
                        &ex6_feature_enhancement_coefficient, 0., 10.)) {
            ex6.coefficient = ex6_feature_enhancement_coefficient;
        };

        if (Button("Apply smoothing", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                ex6.prepare(*cur_polygon);
                ex6.smooth();
            }
        }

        if (Button("Enhance feature", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");

            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                ex6.prepare(*cur_polygon);
                ex6.enhance_feature();
            }
        };
        
        Text("");
        SetWindowFontScale(1.5f);
        Text("Deformation");
        SetWindowFontScale(1.0f);

        if (Button("Pick vertex for deformation", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                int selected = cur_polygon->ps_mesh->selectVertex();
                picked_verts_deform.push_back(selected);
                picked_pts_deform.push_back(cur_polygon->geometry->vertexPositions[selected]);
                auto pc = polyscope::registerPointCloud("Selected points (deform)", picked_pts_deform);
                pc->setPointColor({0.1,0.1,1.});
            }
        }

        if (Button("Pick fixed vertex", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                int selected = cur_polygon->ps_mesh->selectVertex();
                picked_verts_fixed.push_back(selected);
                picked_pts_fixed.push_back(cur_polygon->geometry->vertexPositions[selected]);
                auto pc = polyscope::registerPointCloud("Selected points (fixed)", picked_pts_fixed);
                pc->setPointColor({1.,1.,0.1});
            }
        }

        if (Button("Clear vertices", full_block_size)) {
            picked_verts_deform.clear();
            picked_pts_deform.clear();
            picked_verts_fixed.clear();
            picked_pts_fixed.clear();
            
            polyscope::removeStructure("Selected points (deform)");
            polyscope::removeStructure("Selected points (fixed)");
            polyscope::removeStructure("Deformation vectors");

        }

        Text("Displacement vector");
        static float x, y, z = 0.;
        if (InputFloat("##X coord.", &x) || InputFloat("##Y coord.", &y) || InputFloat("##Z coord.", &z)) {
            if (!(cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) && !(picked_verts_deform.size() == 0)) {
                displacement_vector = Point{x, y, z};
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                std::vector<Point> deformation_viz(cur_polygon->mesh->nVertices(), Point{0.,0.,0.}); 
                for (auto v : picked_verts_deform){
                    deformation_viz[v] = displacement_vector;
                }
                auto deform_viz_ps = cur_polygon->ps_mesh->addVertexVectorQuantity("Deformation vectors", deformation_viz, VectorType::AMBIENT);
                deform_viz_ps->setVectorColor({1., 0., 0.});
                deform_viz_ps->setEnabled(true);
            }
        }

        if (Button("Run deformation", full_block_size)){

            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();
                ex6.prepare(*cur_polygon);
                ex6.deform_surface(picked_verts_fixed, picked_verts_deform, displacement_vector);
                
                //update the point cloud for the deformed vertices
                picked_pts_deform.clear(); 
                for (int v : picked_verts_deform){
                    picked_pts_deform.push_back(cur_polygon->geometry->vertexPositions[v]);
                    auto pc = polyscope::registerPointCloud("Selected points (deform)", picked_pts_deform);
                    pc->setPointColor({0.1,0.1,1.});
                }
            }
        }   
            
    }

    if (CollapsingHeader("Exercise 7", open)) {
        AlignTextToFramePadding();

        static int ex7_iteration = 10;
        static int ex7_cur_laplacian = 0;
        const char* laplacian_types = "Uniform\0Cotangent\0";

        SetWindowFontScale(1.5f);
        Text("Parameterization");
        SetWindowFontScale(1.0f);


        if (Button("Map boundary to circle", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();

                if (!cur_polygon->mesh->isManifold()) {
                    std::cout << "ERROR: Original mesh is not manifold" << std::endl;
                } else {
                    ex7.prepare_param(*cur_polygon);
                    ex7.map_suface_boundary_to_circle();
                }
                
            }
        }

        Text("Iteration");
        SameLine(second_offset);
        if (SliderInt("## ex7_iteration", &ex7_iteration, 1, 100)) {
            ex7.num_iter = ex7_iteration;
        };

        if (Button("Explicit smooth", full_block_size)) {
            if (ex7.is_mapped) {
                ex7.explicitly_smooth_texture();
            } else {
                polyscope::warning("Click \"Map Boundary to Circle\" first.");
            }
        };
        if (Button("Implicit smooth", full_block_size)) {
            if (ex7.is_mapped) {
                ex7.implicitly_smooth_texture();
            } else {
                polyscope::warning("Click \"Map Boundary to Circle\" first.");
            }
        };

        Text("");
        SetWindowFontScale(1.5f);
        Text("Minimal Surfaces");
        SetWindowFontScale(1.0f);

        Text("Operator type");
        SameLine(second_offset);
        if (Combo("##smoothtype", &ex7_cur_laplacian, laplacian_types, 3)) {
            ex7.smooth_type = ex7_cur_laplacian;
        }

        if (Button("Compute minimal surface", full_block_size)) {
            if (cur_mesh_id < 0 || cur_mesh_id >= static_cast<int>(fs.objnum)) {
                polyscope::warning("Please select a mesh first");
            } else {
                PolygonInstance* cur_polygon = fs.getCurrentPolygon();

                ex7.prepare_minimal(*cur_polygon);
                ex7.minimal_surface();
            }
            
        };

    }

    End();
}
