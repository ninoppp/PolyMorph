#ifndef GRID_H
#define GRID_H

#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>

struct Index {
  int i; 
  int j;
  Index(int i, int j): i(i), j(j) {}
};

template<typename T>
struct Grid {
    std::vector<std::vector<T>> data;
    Grid (size_t Nx, size_t Ny) : data(Nx, std::vector<T>(Ny)) {}
    Grid (size_t Nx, size_t Ny, T value) : data(Nx, std::vector<T>(Ny, value)) {}
    Grid () { Grid(100, 100); } // ToDo: get rid of magic numbers
    T& operator()(int i, int j) { return data[i][j]; }
    T& operator()(Index idx) { return data[idx.i][idx.j]; }
    size_t sizeX() const { return data.size(); }
    size_t sizeY() const { return data.empty() ? 0 : data[0].size(); }
    size_t sizeZ() const; // only defined for T=vector
    std::string to_vtk(std::string name);
};

template<typename T>
size_t Grid<T>::sizeZ() const {
    return 0;
}

template<>
size_t Grid<std::vector<double>>::sizeZ() const {
    return data[0][0].size();
}

template<typename T>
std::string Grid<T>::to_vtk(std::string name) { // arr_size ignored  
    std::stringstream xml;
    xml << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">" << std::endl;
    for (int i = 0; i < sizeX(); i++) {
        for (int j = 0; j < sizeY(); j++) {
            xml << data[i][j] << " ";
        }
        xml << std::endl;
    }
    xml << "</DataArray>" << std::endl;
    return xml.str();
}

template<>  // TODO use vector components in vtk
std::string Grid<std::vector<double>>::to_vtk(std::string name) {
    std::stringstream xml;
    for (int k = 0; k < sizeZ(); k++) {
        xml << "<DataArray type=\"Float64\" Name=\"" << name << std::to_string(k) << "\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < sizeX(); i++) {
            for (int j = 0; j < sizeY(); j++) {
                xml << data[i][j][k] << " ";
            }
            xml << std::endl;
        }
        xml << "</DataArray>" << std::endl;
    }
    return xml.str();
}

#endif