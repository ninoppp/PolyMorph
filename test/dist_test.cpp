#include <random>
#include <cstdio>
#include <cmath>
#include <vector>

std::mt19937 rng;
std::lognormal_distribution<double> dist(std::log(50), 0.2);
std::vector<double> samples;

int main() {
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