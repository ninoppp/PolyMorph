#include "ensembleController.h"

/** [Usage example]
*
* @brief Generate tissues with number of polygons = 1, 10, 100, 1000, 10000, 100000
*/

int main() {
    welcome();
    write_config("generate_by_num");
    assert(kh != 0); // adhesion stiffness should be non-zero
    for (int num_polygons = 1e4; num_polygons < 1e6; num_polygons *= 10) {
        Ensemble ensemble = EnsembleController::grow_tissue_by_num(num_polygons);
        ensemble.write_OFF("ensemble/tissues_bynum/" + std::to_string(num_polygons) + ".off");
        ensemble.output(num_polygons);
    }
    return 0;
}