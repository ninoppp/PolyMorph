
#!/bin/bash
g++ -fopenmp -O3 -o polyhoop polyhoop.cpp
perl coplanar_circles.pl 10 1 1 1 # 10 vert, 1 radius, 1 circles, 1 distance 
OMP_NUM_THREADS=8 ./polyhoop

# Create a new directory for the output and move all files
mkdir -p out/$(date +%Y-%m-%d_%H-%M)
mv *.vtp out/$(date +%Y-%m-%d_%H-%M)/
mv *.vts out/$(date +%Y-%m-%d_%H-%M)/