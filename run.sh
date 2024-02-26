
#!/bin/bash
g++ -fopenmp -O3 -o polyhoop polyhoop.cpp
perl coplanar_circles.pl 10 1 5 1 # 10 vert/circle, 1 radius, 5 circles, 1 distance 
OMP_NUM_THREADS=8 ./polyhoop

# Create a new directory for the output
mkdir -p out/$(date +%Y-%m-%d_%H-%M-%S)
mv *.vtp out/$(date +%Y-%m-%d_%H-%M-%S)/