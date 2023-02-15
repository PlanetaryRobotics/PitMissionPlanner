#pragma once
#ifndef TERRAINMAP_H
#define TERRAINMAP_H

#include "tinyexr.h"
#include <fmt/format.h>
#include <cmath>
#include <fstream>
#include <vector>

template <typename T>
class TerrainMap {
public:
    size_t rows;
    size_t cols;
    double pitch;

    TerrainMap(double mapWidth, double mapHeight, double mapPitch) :
     rows(std::ceil(mapHeight/mapPitch)), cols(std::ceil(mapWidth/mapPitch)), pitch(mapPitch), _data(rows*cols, 0) {}

    double width() const {
        return cols*pitch;
    }
    double height() const {
        return rows*pitch;
    }

    T& atIJ(int i, int j) { return operator()(i,j); }
    T& atXY(double x, double y) {
        int j = xCoordToGridIndex(x);
        int i = yCoordToGridIndex(y);
        return operator()(i,j);
    }
    const T& atIJ(int i, int j) const { return operator()(i,j); }
    const T& atXY(int x, int y) const {
        int j = xCoordToGridIndex(x);
        int i = yCoordToGridIndex(y);
        return operator()(i,j);
    }

    double gridIndexToXCoord(size_t j) const { return pitch/2.0 + j*pitch; }
    double gridIndexToYCoord(size_t i) const { return pitch/2.0 + i*pitch; }

    size_t xCoordToGridIndex(double x) const { return x/pitch; }
    size_t yCoordToGridIndex(double y) const { return y/pitch; }

    T& operator()(size_t i, size_t j) { return _data[i*cols+j]; }
    const T& operator()(size_t i, size_t j) const { return _data[i*cols+j]; }

    const std::vector<T>& data() const { return _data; }

private:
    std::vector<T> _data;
};

using TerrainMapFloat = TerrainMap<float>;
using TerrainMapU8 = TerrainMap<uint8_t>;

template<typename T>
void drawCircle(TerrainMap<T>& map, double x, double y, double val, double rad) {
    int cj = map.xCoordToGridIndex(x);
    int ci = map.yCoordToGridIndex(y);
    int R = std::ceil(rad / map.pitch);

    for(int i=ci-R; i<=ci+R; ++i) {
        for(int j=cj-R; j<=cj+R; ++j) {
            if( i < 0 || i > map.rows ) { continue; }
            if( j < 0 || j > map.cols ) { continue; }
            if( (i-ci)*(i-ci)+(j-cj)*(j-cj) < R*R ) {
                map(i,j) = static_cast<T>(val);
            }
        }
    }
}

template <typename T>
void saveEXR(const TerrainMap<T>& map, const std::string& filename) {
    std::vector<float> flipped(map.data().size());
    for(int i=0; i<map.rows; ++i) {
        for(int j=0; j<map.cols; ++j) {
            flipped[i*map.cols+j] = static_cast<float>(map.data()[(map.rows-i-1)*map.cols+j]);
        }
    }
    const char *err;
    int success = SaveEXR(flipped.data(), map.cols, map.rows, 1, 0, filename.c_str(), &err); 

    if( success != TINYEXR_SUCCESS ) {
        fmt::print("TINYEXR ERROR: {}\n", err);
    }
}

template <typename T>
void savePFM(const TerrainMap<T>& map, const std::string& filename) {
    std::ofstream os(filename, std::ios::out);
    os << fmt::format("Pf\n");
    os << fmt::format("{} {}\n", map.cols, map.rows);
    os << fmt::format("-1\n");
    for(int i=0; i<map.rows; ++i) {
        for(int j=0; j<map.cols; ++j) {
            float f = static_cast<float>(map.data()[i*map.cols+j]);
            os.write(reinterpret_cast<const char*>(&f), sizeof(float));
        }
    }
    os.close();
}


#endif // TERRAINMAP_H
