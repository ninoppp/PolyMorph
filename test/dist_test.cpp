#include <random>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>

const double mu = 3;
const double CV = 0.1;
const double lnCV = std::log(1 + CV * CV);
std::mt19937 rng;
std::lognormal_distribution<double> dist(std::log(mu) - lnCV/2, std::sqrt(lnCV));
std::vector<double> samples;

int main() {
  rng.seed(7200);
  std::cout << "lnCV: " << lnCV << " mu-lnCV/2: " << mu - lnCV/2 << std::endl;

  for (int i = 0; i < 100; i++) {
    double sample = dist(rng);
    samples.push_back(sample);
    printf("Sample %d: %f\n", i, sample);
  }

  double sum = 0.0;
  for (double sample : samples) {
    sum += sample;
  }
  double average = sum / samples.size();
  printf("Average: %f\n", average);
}