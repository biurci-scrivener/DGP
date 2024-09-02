#ifndef IMGUI_DEFINE_MATH_OPERATORS
#define IMGUI_DEFINE_MATH_OPERATORS

#include "ImGuiFileDialog.h"
#include "imgui.h"
#include "FileSystem.h"
#include "ex1.h"
#include "ex2.h"
#include "ex3.h"
#include "ex4.h"
#include "ex6.h"
#include "ex7.h"

// 全局 FileSystem 对象
static FileSystem fs;
static Exercise1 ex1;
static Exercise2 ex2;
static Exercise3 ex3;
static Exercise4 ex4;
static Exercise6 ex6;
static Exercise7 ex7;
using namespace ImGui;
using Point = Vector3;

namespace Panel {
    static bool VectorOfStringGetter(void* vec, int idx, const char** out_text);
    void DGPPanel();
}

#endif