
#!/bin/bash

move_files() {
    echo "moving output files to desktop ..."
    local timestamp=$(date +%Y-%m-%d_%H-%M)
    local path="/mnt/c/Users/muell/Desktop/PolymorphOutput" # change this to your desktop path
    mkdir $path/$timestamp
    mv *.vtp $path/$timestamp/
    mv *.vts $path/$timestamp/ 
}

cleanup_files() {
    echo "cleaning up leftover output files ..."
    rm *.vtp
    rm *.vts
}

compile_polyhoop() {
    echo "compiling polyhoop ... "
    g++ -fopenmp -O3 -o polyhoop.out polyhoop.cpp
}

generate_ensemble() {
    echo "generating ensemble ... "
    perl coplanar_circles.pl 10 1 1 1 > ensemble.off # 10 vert, 1 radius, 1 circles, 1 distance 
}

run_polyhoop() {
    echo "running polyhoop in parallel... "
    OMP_NUM_THREADS=8 ./polyhoop.out
}

run_test() {
    g++ -fopenmp -O3 -o test/a.out test/unit_tests.cpp 
    ./test/a.out
}

plot_steady_state() {
    echo "plotting steady state ..."
    python3 test/plot.py
}

if [ "$1" = "c" ]; then
    echo "running cleanup only ..."
    cleanup_files

elif [ "$1" = "t" ]; then
    echo "running test only ..."
    run_test
    plot_steady_state
    
else
    echo "running full pipeline ..."
    compile_polyhoop
    #generate_ensemble
    run_polyhoop
    move_files
fi

echo "all done."