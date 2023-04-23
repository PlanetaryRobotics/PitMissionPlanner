#pragma once
#ifndef BUILDMAPS_H
#define BUILDMAPS_H
#include "config.h"
#include "path.h"
#include "terrainmap.h"
#include "terrainmesh.h"
#include "vantage.h"
#include <tuple>

extern PlannerConfiguration config;

std::tuple<TerrainMapFloat, TerrainMapFloat, SlopeAtlas> buildTerrainMaps(const TerrainMesh &tmesh, const double mapPitch);
TerrainMapU8 buildReachabilityMap(const TerrainMapFloat &commsMap, const SlopeAtlas &slopeAtlas, const Path::State &landingSite);
TerrainMapFloat buildCommsMap(const TerrainMesh &tmesh, const TerrainMapFloat &elevationMap, const Path::State &landingSite);
TerrainMapFloat buildCoverageMap(const TerrainMesh &tmesh, const TerrainMapFloat &elevationMap, const std::vector<Vantage> &vantages);

#endif // BUILDMAPS_H
