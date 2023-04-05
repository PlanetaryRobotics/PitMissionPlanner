#pragma once
#ifndef PATH_H
#define PATH_H

#include "terrainmap.h"
#include "config.h"

extern PlannerConfiguration config;
enum class Direction : int8_t { N = 0, NE, E, SE, S, SW, W, NW };

struct Path {
    struct State {
        int i = 0;
        int j = 0;
        Direction d = Direction::N;

        bool operator==(const State& rhs) const {
            return rhs.i==i && rhs.j==j && d==rhs.d;
        }
        bool operator!=(const State& rhs) const {
            return !(*this == rhs);
        }

        bool sameLocation(const State& rhs) { return rhs.i==i && rhs.j==j; }
        bool sameDirection(const State& rhs) { return d==rhs.d; }
    };

    std::vector<State> states;

    double cost = 0.0;
    double dist = 0.0;
};

struct SlopeAtlas {
    TerrainMapFloat absolute;
    TerrainMapFloat north;
    TerrainMapFloat northeast;
    TerrainMapFloat east;
    TerrainMapFloat southeast;

    SlopeAtlas(double mapWidth, double mapHeight, double mapPitch) :
        absolute(mapWidth, mapHeight, mapPitch), 
        north(mapWidth, mapHeight, mapPitch), 
        northeast(mapWidth, mapHeight, mapPitch), 
        east(mapWidth, mapHeight, mapPitch), 
        southeast(mapWidth, mapHeight, mapPitch) {}

        double j2x(size_t j) const { return absolute.j2x(j); }
        double i2y(size_t i) const { return absolute.i2y(i); }

        size_t x2j(double x) const { return absolute.x2j(x); }
        size_t y2i(double y) const { return absolute.y2i(y); }

        double width() const { return absolute.width(); }
        double height() const { return absolute.height(); }

        size_t rows() const { return absolute.rows; }
        size_t cols() const { return absolute.cols; }
        double pitch() const { return absolute.pitch; }
};

namespace std {
    template <>
    struct hash<Path::State> {
        std::size_t operator()(const Path::State& p) const {
            // pack two int32 into one uint64
            std::size_t res;
            res  = (std::size_t)p.i << 32;
            res |= (std::size_t)p.j & 0x00000000FFFFFFFF;
            // use boost hash_combine to mix in the direction as well.
            res ^= (std::size_t)p.d + 0x9e3779b9 + (res << 6) + (res >> 2);
            return res;
        }
    };
}

double directionToDegrees(const Direction& d);

std::string directionToString(const Direction& d);

Path assembleRoute(const std::vector<int>& route, const std::vector<std::vector<Path>> paths);

Path append(const Path& a, const Path& b);

Path reverse(const Path& path);

uint8_t setBit(uint8_t byte, uint8_t bit);
uint8_t clearBit(uint8_t byte, uint8_t bit);
bool checkBit(uint8_t byte, uint8_t bit);
bool checkMapBit(const TerrainMapU8& map, const Path::State& state);
void setMapBit(TerrainMapU8& map, const Path::State& state);
std::vector<Path::State> getSuccessors(const Path::State& p, const SlopeAtlas& slopeAtlas, const TerrainMapFloat& commsMap);
double computeLongitudinalSlope(const Path::State& s, const SlopeAtlas& slopeAtlas);
double computeLateralSlope(const Path::State& s, const SlopeAtlas& slopeAtlas);
double octileDistance(const Path::State& a, const Path::State& b);
Direction directionFromTo(const Path::State& from, const Path::State& to);

#endif 
