#pragma once
#ifndef MAP_H
#define MAP_H
#include <tuple>
#include "terrainmap.h"
#include "terrainmesh.h"
#include "path.h"
#include "config.h"

    extern PlannerConfiguration config;
    std::tuple<TerrainMapFloat,TerrainMapFloat,SlopeAtlas>
        buildTerrainMaps(const TerrainMesh& tmesh, const double mapPitch);
    TerrainMapU8 buildReachabilityMap(const TerrainMapFloat& commsMap,
                                  const SlopeAtlas& slopeAtlas,
                                  const Path::State& landingSite);
    TerrainMapFloat buildCommsMap(const TerrainMesh& tmesh,
                              const TerrainMapFloat& elevationMap,
                              const Path::State& landingSite);
#endif