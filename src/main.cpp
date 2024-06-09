#include "experiments.h"
#include "ensembleController.h"
#include "utils.h"

int main(int argc, char* argv[]) {
  welcome();
  validate_parameters();
  write_config();
  //default_testrun();
  const char* nodeID_str = getenv("SLURM_NODEID");
  if (!nodeID_str) {
      std::cerr << "SLURM_NODEID is not set." << std::endl;
      return 1;
  }
  int nodeID = std::atoi(nodeID_str);
  std::cout << "Node ID: " << nodeID << std::endl;
}

