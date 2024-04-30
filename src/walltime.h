#include <sys/time.h>
#include <iostream>

double walltime(){
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        std::cout << "Unexpected error in get_wall_time" << std::endl;
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}