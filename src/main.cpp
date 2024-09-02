#include "include/Panel.h"
#include "polyscope/polyscope.h"

int main() {
    polyscope::init();
    //polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    //polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
    polyscope::options::groundPlaneHeightFactor = 0.1;
    // polyscope::state::userCallback = [&]() { Panel::DGPPanel(); };
    polyscope::state::userCallback = Panel::DGPPanel;
    polyscope::buildPickGui();
    polyscope::show();
    return EXIT_SUCCESS;
}