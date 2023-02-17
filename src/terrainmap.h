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
        int j = x2j(x);
        int i = y2i(y);
        return operator()(i,j);
    }
    const T& atIJ(int i, int j) const { return operator()(i,j); }
    const T& atXY(int x, int y) const {
        int j = x2j(x);
        int i = y2i(y);
        return operator()(i,j);
    }

    double j2x(size_t j) const { return (j+0.5)*pitch; }
    double i2y(size_t i) const { return (rows-i-0.5)*pitch; }

    size_t x2j(double x) const { return x/pitch - 0.5; }
    size_t y2i(double y) const { return rows - 0.5 - y/pitch; }

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
    int cj = map.x2j(x);
    int ci = map.y2i(y);
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
    const char *err;
    int success = SaveEXR(map.data().data(), map.cols, map.rows, 1, 0, filename.c_str(), &err); 

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
