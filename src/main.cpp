#include "terrainmesh.h"
#include "terrainmap.h"
#include "argh.h"
#include <fmt/format.h>
#include <optional>
#include <unordered_set>
#include <iostream>


struct PlannerConfiguration {
    std::string meshfile = "../meshes/lmp_science.ply";
    double mapPitch = 1.0;
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
        fmt::print("Building Terrain Maps [{}/{}]\n", i, slopeMap.rows);
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

TerrainMap<char> buildReachabilityMap(const TerrainMap<float>& slopeMap,
                                       double siteX, double siteY,
                                       double roverMaxSlope) {
    TerrainMap<char> reachMap(slopeMap.width(), slopeMap.height(), slopeMap.pitch);

    int rows = slopeMap.rows;
    int cols = slopeMap.cols;

    std::vector<bool> visited(rows*cols, false);
    int ri = slopeMap.yCoordToGridIndex(siteY);
    int rj = slopeMap.xCoordToGridIndex(siteX);

    std::vector<int> open;
    open.push_back(ri*cols+rj);

    int iterations = 0;
    while( !open.empty() ) {
        int curr = open.back();
        open.pop_back();

        if( iterations++ % (1<<16) ) {
            fmt::print("Building Reach Map {}\n", iterations);
        }

        if( visited[curr] ) { continue; }
        visited[curr] = true;

        int ci = curr / cols;
        int cj = curr % cols;
        reachMap(ci, cj) = 1;

        // Walk east. Add north/south neighbors to the open stack if you can.
        {
            int dj=1;
            while( cj+dj < cols ) {
                int n = ci*cols+cj+dj;
                if( visited[n] || slopeMap(ci, cj+dj) > roverMaxSlope ) {
                    break;
                }
                // Can you add the north neighbor to the open set?
                if( ci+1 < rows && !visited[(ci+1)*cols+cj+dj] && slopeMap(ci+1, cj+dj) <= roverMaxSlope) {
                    open.push_back((ci+1)*cols+cj+dj);
                }
                // Can you add the south neighbor to the open set?
                if( ci-1 >= 0 && !visited[(ci-1)*cols+cj+dj] && slopeMap(ci-1, cj+dj) <= roverMaxSlope) {
                    open.push_back((ci-1)*cols+cj+dj);
                }
                reachMap(ci, cj+dj) = 1;
                visited[n] = true;
                dj++;
            }
        }

        // Walk west. Add north/south neighbors to the open stack if you can.
        {
            int dj=1;
            while( cj-dj >= 0 ) {
                int n = ci*cols+cj-dj;
                if( visited[n] || slopeMap(ci, cj-dj) > roverMaxSlope ) {
                    break;
                }
                // Can you add the north neighbor to the open set?
                if( ci+1 < rows && !visited[(ci+1)*cols+cj-dj] && slopeMap(ci+1, cj-dj) <= roverMaxSlope) {
                    open.push_back((ci+1)*cols+cj-dj);
                }
                // Can you add the south neighbor to the open set?
                if( ci-1 >= 0 && !visited[(ci-1)*cols+cj-dj] && slopeMap(ci-1, cj-dj) <= roverMaxSlope) {
                    open.push_back((ci-1)*cols+cj-dj);
                }
                reachMap(ci, cj-dj) = 1;
                visited[n] = true;
                dj++;
            }
        }
    }
    return reachMap;
}

TerrainMap<char> buildCommsMap(const TerrainMesh& tmesh,
                               const TerrainMap<float>& elevationMap,
                               double siteX, double siteY,
                               double landerHeight, double roverHeight) {

    TerrainMap<char> commsMap(elevationMap.width(),
                              elevationMap.height(),
                              elevationMap.pitch);

    const int li = elevationMap.yCoordToGridIndex(siteY);
    const int lj = elevationMap.xCoordToGridIndex(siteX);
    const double landerX = siteX;
    const double landerY = siteY;
    const double landerZ = elevationMap(li, lj) + landerHeight;

    for(int ri=0; ri<commsMap.cols; ++ri) {
        for(int rj=0; rj<commsMap.rows; ++rj) {
            TerrainMesh::Ray ray;
            const double roverX = elevationMap.gridIndexToXCoord(rj);
            const double roverY = elevationMap.gridIndexToYCoord(ri);
            const double roverZ = elevationMap(ri, rj) + roverHeight;
            ray.oX = landerX;
            ray.oY = landerY;
            ray.oZ = landerZ;
            ray.dX = roverX-landerX;
            ray.dY = roverY-landerY;
            ray.dZ = roverZ-landerZ;
            const auto hit = tmesh.raytrace(ray);
            if( !hit ) {
                commsMap(ri, rj) = 1;
            } 
            double roverDist2 = (landerX-roverX)*(landerX-roverX)+
                                (landerY-roverY)*(landerY-roverY)+
                                (landerZ-roverZ)*(landerZ-roverZ);
            double hitDist2 = (ray.oX-hit->x)*(ray.oX-hit->x) + (ray.oY-hit->y)*(ray.oY-hit->y) + (ray.oZ-hit->z)*(ray.oZ-hit->z);
            if ( hitDist2 >= 0.98*roverDist2 ) {
                commsMap(ri, rj) = 1;
            }
        }
    }

    return commsMap;
}

int main(int argc, char* argv[]) {

    const auto config = parseCommandLine(argc, argv);
    if( !config ) { return 0; }

    TerrainMesh tmesh(config->meshfile);

    auto [elevationMap, slopeMap, priorityMap] = buildTerrainMaps(tmesh, config->mapPitch);
    elevationMap.savePFM("elevation.pfm");
    slopeMap.savePFM("slope.pfm");
    priorityMap.savePFM("priority.pfm");

    double siteX = 825;
    double siteY = 1000;
    double landerHeight = 3.0;
    double roverHeight = 1.0;

    {
        TerrainMap<float> siteMap = slopeMap;
        auto siteI = siteMap.yCoordToGridIndex(siteY);
        auto siteJ = siteMap.xCoordToGridIndex(siteX);
        int radius = 5.0 / siteMap.pitch;
        for(int i=siteI-radius; i<siteI+radius; ++i) {
            for(int j=siteJ-radius; j<siteJ+radius; ++j) {
                if( (i-siteI)*(i-siteI)+(j-siteJ)*(j-siteJ) < radius*radius ) {
                    siteMap(i,j) = 180;
                }
            }
        }
        siteMap.savePFM("site.pfm");
    }

    TerrainMap reachMap = buildReachabilityMap(slopeMap, siteX, siteY, 13);
    reachMap.savePFM("reach.pfm");

    TerrainMap commsMap = buildCommsMap(tmesh, elevationMap,
                                            siteX, siteY,
                                            landerHeight, roverHeight);
    commsMap.savePFM("comms.pfm");

    return 0;
}
