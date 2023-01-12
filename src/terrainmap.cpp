#include "terrainmap.h"
#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"
#include <fmt/format.h>

void TerrainMap::saveEXR(const std::string& filename) const {
    std::vector<float> flipped(data.size());
    for(int i=0; i<rows; ++i) {
        for(int j=0; j<cols; ++j) {
            flipped[i*cols+j] = data[(rows-i-1)*cols+j];
        }
    }
    const char *err;
    int success = SaveEXR(flipped.data(), cols, rows, 1, 0, filename.c_str(), &err); 

    if( success != TINYEXR_SUCCESS ) {
        fmt::print("TINYEXR ERROR: {}\n", err);
    }
}

void TerrainMap::savePFM(const std::string& filename) const {
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

void TerrainMap::drawCircle(double x, double y, double val, double rad) {
    int cj = xCoordToGridIndex(x);
    int ci = yCoordToGridIndex(y);
    int R = std::ceil(rad / pitch);

    for(int i=ci-R; i<=ci+R; ++i) {
        for(int j=cj-R; j<=cj+R; ++j) {
            if( i < 0 || i > rows ) { continue; }
            if( j < 0 || j > cols ) { continue; }
            if( (i-ci)*(i-ci)+(j-cj)*(j-cj) < R*R ) {
                this->operator()(i,j) = val;
            }
        }
    }
}

