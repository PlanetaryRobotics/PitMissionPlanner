#include "terrainmesh.h"
#include "terrainmap.h"
#include "argh.h"
#include <fmt/format.h>
#include <optional>


struct PlannerConfiguration {
    std::string meshfile = "../meshes/lmp_science.ply";
    double mapPitch = 2.0;
};

std::optional<PlannerConfiguration> parseCommandLine(int argc, char* argv[]) {
    PlannerConfiguration cfg;

    argh::parser cmdl(argc, argv);

    if( cmdl[{ "-h", "--help" }] ) {
        fmt::print("usage: planranger [-pitch p] meshfile\n");
        fmt::print("\tmeshfile: a .ply format map of terrain surrounding a lunar pit.\n");
        fmt::print("\tpitch: the grid spacing (in meters) to use for generated maps.\n");
        return {};
    }
    cmdl(1, cfg.meshfile) >> cfg.meshfile;
    cmdl("pitch", cfg.mapPitch) >> cfg.mapPitch;
    return std::make_optional(cfg);
}

std::tuple<TerrainMap<float>,TerrainMap<float>,TerrainMap<float>>
buildTerrainMaps(const TerrainMesh& tmesh, const double mapPitch) {
    const double mapX = tmesh.maxX()-tmesh.minX();
    const double mapY = tmesh.maxY()-tmesh.minY();
    
    TerrainMap<float> elevationMap(mapX, mapY, mapPitch);
    TerrainMap<float> slopeMap(mapX, mapY, mapPitch);
    TerrainMap<float> priorityMap(mapX, mapY, mapPitch);

    const double maxZ = tmesh.maxZ();
    const double minZ = tmesh.minZ();

    for(int i=0; i<slopeMap.rows; ++i) {
        fmt::print("Generating Maps {}/{}\n", i, slopeMap.rows);
    #pragma omp parallel for
        for(int j=0; j<slopeMap.cols; ++j) {
            TerrainMesh::Ray ray;
            ray.oX = slopeMap.gridIndexToXCoord(j);
            ray.oY = slopeMap.gridIndexToYCoord(i);
            ray.oZ = maxZ + 10.0;
            ray.dX = 0.0; ray.dY = 0.0; ray.dZ = -1.0;
            const auto hit = tmesh.raytrace(ray);
            if( hit ) {
                elevationMap(i,j) = hit->z;
                slopeMap(i,j) = 180.0/M_PI * std::acos(0.0*hit->nx+0.0*hit->ny+1.0*hit->nz);
                priorityMap(i,j) = hit->priority;
            } else {
                elevationMap(i,j) = minZ;
                slopeMap(i,j) = 0.0;
                priorityMap(i,j) = 0.0;
            }
        }
    }
    return std::make_tuple(elevationMap, slopeMap, priorityMap);
}

int main(int argc, char* argv[]) {

    const auto config = parseCommandLine(argc, argv);
    if( !config ) { return 0; }

    TerrainMesh tmesh(config->meshfile);

    auto [elevationMap, slopeMap, priorityMap] = buildTerrainMaps(tmesh, config->mapPitch);
    elevationMap.savePFM("elevation.pfm");
    slopeMap.savePFM("slope.pfm");
    priorityMap.savePFM("priority.pfm");

    //TerrainMap commsMap = buildCommsMap(tmesh, landingSite, landerHeight, roverHeight);
    //TerrainMap reachMap = buildReachabilityMap(slopeMap, landingSite, roverMaxSlope);

    return 0;
}
