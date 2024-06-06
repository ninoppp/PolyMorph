#!/bin/bash

#output="./out"
output="/mnt/c/Users/muell/Desktop/PolymorphOutput" # nico or muell
OMP_NUM_THREADS=8

move_files() {
    echo "moving output files to folder ..."
    local timestamp=$(date +%Y-%m-%d_%H-%M)
    mkdir $output/$timestamp
    mv *.vtp *.vts *.cfg "$output/$timestamp/"
}

cleanup_files() {
    echo "cleaning up leftover output files ..."
    rm *.vtp *.vts *.cfg
}

run() {
    echo "compiling src/$1.cpp ..."
    make clean
    make $1
    export OMP_NUM_THREADS
    echo "running $1 with $OMP_NUM_THREADS threads... "
    ./polymorph
}

if [ "$1" = "c" ]; then
    cleanup_files
    
elif [ "$1" = "m" ]; then
    move_files

elif [ -n "$1" ]; then
    run "$1"
    move_files

else
    run "main"
    move_files
fi

echo "all done."