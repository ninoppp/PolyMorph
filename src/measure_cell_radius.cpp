#include "ensembleController.h"

int main() {
    double diameter = 0;
    for (int width : {6, 8, 10, 13, 17, 20, 25, 30, 35, 40}) {
        Domain domain(-20, -20, 20, 20);
        std::string filename = "ensemble/tissues_varwidth/" + std::to_string(width) + "_2.off";
        Ensemble ensemble(filename.c_str(), domain);
        // measure average area
        double total_area = 0;
        for (int p = 0; p < ensemble.polygons.size(); ++p) {
            total_area += ensemble.polygons[p].area();
        }
        double average_area = total_area / ensemble.polygons.size();
        double average_radius = std::sqrt(average_area / M_PI);
        double average_diameter = 2 * average_radius;
        diameter += average_diameter;
        std::cout << "Width=" << width << " Average Diameter=" << average_diameter << std::endl;
    }
    std::cout << "Full Average Diameter=" << diameter / 10 << std::endl;
    return 0;
}

// average diameter = 1.31
// radius = 0.655