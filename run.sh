
#!/bin/bash
echo "cleaning up leftover output files ..."
rm *.vtp
rm *.vts

echo "compiling polyhoop ... "
g++ -fopenmp -O3 -o polyhoop.out polyhoop.cpp
echo "generating ensemble ... "
perl coplanar_circles.pl 10 1 1 1 > ensemble.off # 10 vert, 1 radius, 1 circles, 1 distance 

echo " running polyhoop ... "
OMP_NUM_THREADS=8 ./polyhoop.out

# Create a new directory for the output and move all files
mkdir -p out/$(date +%Y-%m-%d_%H-%M)
mv *.vtp out/$(date +%Y-%m-%d_%H-%M)/
mv *.vts out/$(date +%Y-%m-%d_%H-%M)/
echo "all done."