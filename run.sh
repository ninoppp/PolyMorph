
#!/bin/bash

move_files() {
    echo "moving output files to desktop ..."
    local timestamp=$(date +%Y-%m-%d_%H-%M)
    local path="/mnt/c/Users/Nico/Desktop/PolymorphOutput" # change this to your desktop path (muell)
    mkdir $path/$timestamp
    mv *.vtp $path/$timestamp/
    mv *.vts $path/$timestamp/ 
    mv *.cfg $path/$timestamp/
}

cleanup_files() {
    echo "cleaning up leftover output files ..."
    rm *.vtp
    rm *.vts
    rm *.cfg
}

compile_polyhoop() {
    echo "compiling polymorph ... "
    rm -f polymorph.out
    g++ -fopenmp -O3 -o polymorph.out src/main.cpp
}

generate_ensemble() {
    echo "generating ensemble ... "
    perl coplanar_circles.pl 10 1 1 1 > ensemble.off # 10 vert, 1 radius, 1 circles, 1 distance 
}

run_polyhoop() {
    OMP_NUM_THREADS=8
    export OMP_NUM_THREADS
    echo "running polymorph with $OMP_NUM_THREADS threads... "
    ./polymorph.out
}

compile_and_run_test() {
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
    compile_and_run_test
    plot_steady_state
    
elif [ "$1" = "m" ]; then
    move_files

else
    echo "running full pipeline ..."
    compile_polyhoop
    #generate_ensemble
    run_polyhoop
    move_files
fi

echo "all done."