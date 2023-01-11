#pragma once
#ifndef TERRAINMAP_H
#define TERRAINMAP_H

#include <fmt/format.h>
#include <cmath>
#include <fstream>

class TerrainMap {
public:
    size_t rows;
    size_t cols;
    double pitch;

    TerrainMap(double mapWidth, double mapHeight, double mapPitch) :
     rows(std::ceil(mapHeight/mapPitch)), cols(std::ceil(mapWidth/mapPitch)), pitch(mapPitch), data(rows*cols, 0) {}

    double width() const {
        return cols*pitch;
    }
    double height() const {
        return rows*pitch;
    }

    double gridIndexToXCoord(size_t j) const { return pitch/2.0 + j*pitch; }
    double gridIndexToYCoord(size_t i) const { return pitch/2.0 + i*pitch; }

    size_t xCoordToGridIndex(double x) const { return x/pitch; }
    size_t yCoordToGridIndex(double y) const { return y/pitch; }

    float& operator()(size_t i, size_t j) { return data[i*cols+j]; }
    const float& operator()(size_t i, size_t j) const { return data[i*cols+j]; }

    void savePFM(const std::string& filename) const {
        std::ofstream os(filename, std::ios::out);
        os << fmt::format("Pf\n");
        os << fmt::format("{} {}\n", cols, rows);
        os << fmt::format("-1\n");
        for(int i=0; i<rows; ++i) {
            for(int j=0; j<cols; ++j) {
                float f = static_cast<float>(data[i*cols+j]);
                os.write(reinterpret_cast<const char*>(&f), sizeof(float));
            }
        }
        os.close();
    }

private:
    std::vector<float> data;
};

#endif // TERRAINMAP_H
