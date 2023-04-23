#pragma once
#include "path.h"
#include <vector>

#include "terrainmap.h"
#include "terrainmesh.h"
struct Probe {
    double x, y, z;
    double priority;
};
struct Vantage {
    double x, y, z;
    Direction dir = Direction::N;
    std::vector<bool> coverage;
    double totalCoverage = 0;
};
std::pair<std::vector<Probe>, TerrainMapFloat> generateVisibilityProbes(const TerrainMapFloat &priorityMap,
                                                                        const TerrainMapFloat &elevationMap);
std::pair<std::vector<Vantage>, TerrainMapFloat> generateVantageCandidates(const TerrainMesh &tmesh, const TerrainMapU8 &reachMap,
                                                                           const TerrainMapFloat &elevationMap,
                                                                           const std::vector<Probe> probes);
std::vector<Vantage> selectVantages(const std::vector<Vantage> &candidates, const std::vector<Probe> &probes);
std::vector<Vantage> sortCCW(const std::vector<Vantage> vantages, double siteX, double siteY);