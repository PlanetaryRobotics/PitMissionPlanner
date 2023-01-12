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
