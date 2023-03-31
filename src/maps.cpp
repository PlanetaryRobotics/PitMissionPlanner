#include "maps.h"
#include <tuple>
//#include "terrainmap.h"
//#include "terrainmesh.h"
std::tuple<TerrainMapFloat,TerrainMapFloat,SlopeAtlas>
buildTerrainMaps(const TerrainMesh& tmesh, const double mapPitch) {
    const double mapX = tmesh.maxX()-tmesh.minX();
    const double mapY = tmesh.maxY()-tmesh.minY();
    
    TerrainMapFloat elevationMap(mapX, mapY, mapPitch);
    TerrainMapFloat priorityMap(mapX, mapY, mapPitch);
    SlopeAtlas slopeAtlas(mapX, mapY, mapPitch);

    const double maxZ = tmesh.maxZ();
    const double minZ = tmesh.minZ();

    for(int i=0; i<slopeAtlas.rows(); ++i) {
        fmt::print("[{}/{}] Building terrain maps.\n", i, slopeAtlas.rows());
        #pragma omp parallel for
        for(int j=0; j<slopeAtlas.cols(); ++j) {
            TerrainMesh::Ray ray;
            ray.oX = slopeAtlas.j2x(j);
            ray.oY = slopeAtlas.i2y(i);
            ray.oZ = maxZ + 10.0;
            ray.dX = 0.0; ray.dY = 0.0; ray.dZ = -1.0;
            const auto hit = tmesh.raytrace(ray);
            if( hit ) {
                elevationMap(i,j) = hit->z;
                priorityMap(i,j) = hit->priority;

                double slope = 180.0/M_PI * std::acos(0.0*hit->nx+0.0*hit->ny+1.0*hit->nz);
                if (std::isnan(slope)) { slope = 0.0f; }
                slopeAtlas.absolute(i,j) = slope;

                // Project the normal onto the YZ plane and normalize it.
                // Then compute its angle with the +Y axis.
                {
                    // Set x to zero to project onto the YZ plane.
                    double nx_on_yz = 0.0;
                    double ny_on_yz = hit->ny;
                    double nz_on_yz = hit->nz;
                    double nz_on_yz_norm = std::sqrt(nx_on_yz*nx_on_yz + ny_on_yz*ny_on_yz + nz_on_yz*nz_on_yz);
                    // Dot product with +Z axis.
                    double dot = 0*nx_on_yz + 0*ny_on_yz + 1*nz_on_yz;
                    // Compute the angle with the XY plane.
                    slopeAtlas.north(i,j) = 180.0/M_PI * std::acos(dot/nz_on_yz_norm) * ((hit->ny >= 0) ? -1.0 : 1.0);
                }

                // Project the normal onto the XZ plane and normalize it.
                // Then compute its angle with the +X axis.
                {
                    // Set y to zero to project onto the XZ plane.
                    double nx_on_xz = hit->nx;
                    double ny_on_xz = 0.0;
                    double nz_on_xz = hit->nz;
                    double nz_on_xz_norm = std::sqrt(nx_on_xz*nx_on_xz + ny_on_xz*ny_on_xz + nz_on_xz*nz_on_xz);
                    // Dot product with +Z axis.
                    double dot = 0*nx_on_xz + 0*ny_on_xz + 1*nz_on_xz;
                    // Compute the angle with the XY plane.
                    slopeAtlas.east(i,j) = 180.0/M_PI * std::acos(dot/nz_on_xz_norm) * ((hit->nx >= 0) ? -1.0 : 1.0);
                }

                // Compute the directional slope in the northeast and southeast directions.
                {
                    double eSlope = slopeAtlas.east(i,j); double nSlope = slopeAtlas.north(i,j);
                    slopeAtlas.northeast(i,j) = 90 - 180/M_PI * std::acos(0.5 * (std::sin(M_PI/180 * eSlope) + std::sin(M_PI/180 * nSlope)));
                    slopeAtlas.southeast(i,j) = 90 - 180/M_PI * std::acos(0.5 * (std::sin(M_PI/180 * eSlope) - std::sin(M_PI/180 * nSlope)));
                }

            } else {
                elevationMap(i,j) = minZ;
                priorityMap(i,j) = 0.0;

                slopeAtlas.north(i,j)     = 0.0;
                slopeAtlas.northeast(i,j) = 0.0;
                slopeAtlas.east(i,j)      = 0.0;
                slopeAtlas.southeast(i,j) = 0.0;
                slopeAtlas.absolute(i,j)  = 0.0;
            }
        }
    }
    return std::make_tuple(elevationMap, priorityMap, slopeAtlas);
}
TerrainMapU8 buildReachabilityMap(const TerrainMapFloat& commsMap,
                                  const SlopeAtlas& slopeAtlas,
                                  const Path::State& landingSite) {

    TerrainMapU8 reachMap(commsMap.width(), commsMap.height(), commsMap.pitch);

    std::vector<Path::State> open;
    open.push_back(landingSite);

    int iterations = 0;
    while( !open.empty() ) {
        if( iterations++ % (1<<14) == 0 ) {
            fmt::print("[{}] Building reach map.\n", iterations);
        }

        // Pop currState from open.
        const auto currState = open.back();
        open.pop_back();

        // Otherwise, mark it as reached.
        setMapBit(reachMap, currState);

        // Add successors of the current state to the open stack.
        const auto& successors = getSuccessors(currState, slopeAtlas, commsMap);
        for(const auto& succ : successors) {
            if( !checkMapBit(reachMap, succ) ) {
                open.push_back(succ);
            }
        }
    }
    return reachMap;
}

TerrainMapFloat buildCommsMap(const TerrainMesh& tmesh,
                              const TerrainMapFloat& elevationMap,
                              const Path::State& landingSite) {

    TerrainMapFloat commsMap(elevationMap.width(),
                         elevationMap.height(),
                         elevationMap.pitch);

    const int li = landingSite.i;
    const int lj = landingSite.j;
    const double landerX = elevationMap.j2x(lj);
    const double landerY = elevationMap.i2y(li);
    const double landerZ = elevationMap(li, lj) + config.landerHeight;

    for(int ri=0; ri<commsMap.rows; ++ri) {
        for(int rj=0; rj<commsMap.cols; ++rj) {
            TerrainMesh::Ray ray;
            const double roverX = elevationMap.j2x(rj);
            const double roverY = elevationMap.i2y(ri);
            const double roverZ = elevationMap(ri, rj) + config.roverHeight;
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
