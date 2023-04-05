#pragma once
#include "priority_queue.h"
#include "terrainmap.h"
#include "terrainmesh.h"
#include <vector>

#include "path.h"
#include <Eigen/Dense>

std::pair<std::vector<Path>, int> multipathplan(const SlopeAtlas &slopeAtlas, const TerrainMapFloat &commsMap,
                                                const Path::State &start, const std::vector<Path::State> &goals,
                                                const int threadID);

std::pair<std::vector<int>, int64_t> solveTSP(const Eigen::MatrixXd &costs, int depot);
std::vector<int> routeplan(const Eigen::MatrixXd &costs);
std::vector<std::vector<Path>> planAllPairsSLOW(const std::vector<Path::State> &sites, const SlopeAtlas &slopeAtlas,
                                                const TerrainMapFloat &commsMap);
std::vector<std::vector<Path>> planAllPairs(const std::vector<Path::State> &sites, const SlopeAtlas &slopeAtlas,
                                            const TerrainMapFloat &commsMap);
