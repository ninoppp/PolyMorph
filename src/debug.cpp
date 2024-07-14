#include "../include/ensembleController.h"

int main() {
    std::cout << "DEBUGGING NOW" << std::endl;
    double L = 60;
    double W = 50;
    int seed = 32;
    Domain domain(-L/2, -W/2, L/2, W/2);
    Ensemble ensemble("ensemble/tissues_varwidth/30_0.off", domain, seed); 
    Solver solver(domain, dx, Reactions::linearDegradation);
    Interpolator interpolator(ensemble, solver);
    //ensemble.is_producing = [](const Polygon& p) { return std::vector<bool> {p.vertices[0].p == Nr}; };
    ensemble.step();
    // calculate average number of vertices per box
    int num_vert_tot = 0;
    int non_empty_boxes = 0;
    int max_vert = 0;
    for (int bi = 0; bi < ensemble.Nx; bi++) {
        for (int bj = 0; bj < ensemble.Ny; bj++) {
            int b = bi + bj * ensemble.Ny;
            int num_vert = 0; 
            Vertex* v = ensemble.first[b];
            while (v != nullptr) {
                num_vert++;
                v = v->next;
            }
            if (num_vert != 0) { // only consider box if not empty
                num_vert_tot += num_vert;
                ++non_empty_boxes;
            }
            if (num_vert > max_vert) {
                max_vert = num_vert;
            }
        }
    }

    double avg_vert_per_box = 1.0 * num_vert_tot / non_empty_boxes;
    double avg_vert_all_boxes = 1.0 * num_vert_tot / (ensemble.Nx * ensemble.Ny);
    std::cout << "number of boxes: " << ensemble.Nx * ensemble.Ny << std::endl;
    std::cout << "average vertices per non-empty box: " << avg_vert_per_box << std::endl;
    std::cout << "average vertices over all boxes: " << avg_vert_all_boxes << std::endl;
    std::cout << "max vertices in a box: " << max_vert << std::endl;
    return 0;
}