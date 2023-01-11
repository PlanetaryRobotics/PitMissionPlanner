#include "terrainmesh.h"
#include "terrainmap.h"
#include "argh.h"
#include <fmt/format.h>
#include <optional>
#include <unordered_set>
#include <iostream>


struct PlannerConfiguration {
    std::string meshfile = "../meshes/lmp_science.ply";
    double mapPitch      = 1.0;
    double landingSiteX  = 825;
    double landingSiteY  = 1000;
    double landerHeight  = 3.0;
    double roverHeight   = 1.0;
    double roverMaxSlope = 20.0;
    int    numProbes     = 50000;
    int    numCandidates = 100000;
    int    numVantages   = 30;
    double visAngle      = 55;
};

std::optional<PlannerConfiguration> parseCommandLine(int argc, char* argv[]) {
    PlannerConfiguration cfg;

    argh::parser cmdl(argc, argv);

    if( cmdl[{ "-h", "--help" }] ) {
        fmt::print("usage: planranger [-p][-x][-y][-rs][-rh][-lh] meshfile\n");
        fmt::print("\tmeshfile: a .ply format map of terrain surrounding a lunar pit.\n");
        fmt::print("\tp: the grid spacing (in meters) to use for generated maps.\n");
        fmt::print("\tx: the x-coordinate of the landing site.\n");
        fmt::print("\ty: the y-coordinate of the landing site.\n");
        fmt::print("\trs: the rover's max slope capability in degrees.\n");
        fmt::print("\trh: the height of the rover antenna.\n");
        fmt::print("\tlh: the height of the lander antenna.\n");
        fmt::print("\tnp: the number of visibility probes.\n");
        fmt::print("\tnc: the number of candidate vantages to evaluate.\n");
        fmt::print("\tnv: the number of vantages to select from the candidates.\n");
        fmt::print("\tva: the visibility angle [deg] beyond which a view is not counted.\n");
        return {};
    }
    cmdl(1, cfg.meshfile) >> cfg.meshfile;
    cmdl("p", cfg.mapPitch) >> cfg.mapPitch;
    cmdl("x", cfg.landingSiteX) >> cfg.landingSiteX;
    cmdl("y", cfg.landingSiteY) >> cfg.landingSiteY;
    cmdl("rs", cfg.roverMaxSlope) >> cfg.roverMaxSlope;
    cmdl("rh", cfg.roverHeight) >> cfg.roverHeight;
    cmdl("lh", cfg.landerHeight) >> cfg.landerHeight;
    cmdl("np", cfg.numProbes) >> cfg.numProbes;
    cmdl("nc", cfg.numCandidates) >> cfg.numCandidates;
    cmdl("nv", cfg.numVantages) >> cfg.numVantages;
    return std::make_optional(cfg);
}

std::tuple<TerrainMap,TerrainMap,TerrainMap>
buildTerrainMaps(const TerrainMesh& tmesh, const double mapPitch) {
    const double mapX = tmesh.maxX()-tmesh.minX();
    const double mapY = tmesh.maxY()-tmesh.minY();
    
    TerrainMap elevationMap(mapX, mapY, mapPitch);
    TerrainMap slopeMap(mapX, mapY, mapPitch);
    TerrainMap priorityMap(mapX, mapY, mapPitch);

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

TerrainMap buildReachabilityMap(const TerrainMap& safeMap,
                                double siteX, double siteY,
                                double roverMaxSlope) {
    TerrainMap reachMap(safeMap.width(), safeMap.height(), safeMap.pitch);

    int rows = safeMap.rows;
    int cols = safeMap.cols;

    std::vector<bool> visited(rows*cols, false);
    int ri = safeMap.yCoordToGridIndex(siteY);
    int rj = safeMap.xCoordToGridIndex(siteX);

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
                if( visited[n] || safeMap(ci, cj+dj) == 0 ) {
                    break;
                }
                // Can you add the north neighbor to the open set?
                if( ci+1 < rows && !visited[(ci+1)*cols+cj+dj] && safeMap(ci+1, cj+dj)>0 ) {
                    open.push_back((ci+1)*cols+cj+dj);
                }
                // Can you add the south neighbor to the open set?
                if( ci-1 >= 0 && !visited[(ci-1)*cols+cj+dj] && safeMap(ci-1, cj+dj)>0 ) {
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
                if( visited[n] || safeMap(ci, cj-dj) == 0 ) {
                    break;
                }
                // Can you add the north neighbor to the open set?
                if( ci+1 < rows && !visited[(ci+1)*cols+cj-dj] && safeMap(ci+1, cj-dj)>0 ) {
                    open.push_back((ci+1)*cols+cj-dj);
                }
                // Can you add the south neighbor to the open set?
                if( ci-1 >= 0 && !visited[(ci-1)*cols+cj-dj] && safeMap(ci-1, cj-dj)>0 ) {
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

TerrainMap buildCommsMap(const TerrainMesh& tmesh,
                               const TerrainMap& elevationMap,
                               double siteX, double siteY,
                               double landerHeight, double roverHeight) {

    TerrainMap commsMap(elevationMap.width(),
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

TerrainMap buildSafeMap(const TerrainMap& commsMap, const TerrainMap& slopeMap, double roverMaxSlope) {
    TerrainMap safeMap = commsMap;
    for(int i=0; i<safeMap.rows; ++i) {
        for(int j=0; j<safeMap.cols; ++j) {
            safeMap(i,j) = (commsMap(i,j) > 0 && slopeMap(i,j) <= roverMaxSlope) ? 1 : 0;
        }
    }
    return safeMap;
}

struct Probe {
    double x,y,z;
    double priority;
};

std::pair<std::vector<Probe>, TerrainMap>
generateVisibilityProbes(const TerrainMap& priorityMap, const TerrainMap& elevationMap, int numProbes, double roverHeight) {
    std::vector<Probe> probes;
    TerrainMap probeMap(priorityMap.width(), priorityMap.height(), priorityMap.pitch);
    while(probes.size() < numProbes) {
        int i = priorityMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        int j = priorityMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        if( priorityMap(i,j) > 0 ) {
            Probe p;
            p.x = priorityMap.gridIndexToXCoord(j);
            p.y = priorityMap.gridIndexToYCoord(i);
            p.z = elevationMap(i,j);
            p.priority = priorityMap(i,j);
            probeMap(i,j) = p.priority;
            probes.push_back(p);
        }
    }
    return std::make_pair(probes, probeMap);
}

struct Vantage {
    double x, y, z;
    std::vector<bool> coverage;
    double totalCoverage = 0;
};

std::pair<std::vector<Vantage>, TerrainMap>
generateVantageCandidates(const TerrainMesh& tmesh,
                          const TerrainMap& reachMap,
                          const TerrainMap& elevationMap,
                          const std::vector<Probe> probes,
                          int numCandidates, double roverHeight, double visAngle) {
    std::vector<Vantage> candidates;
    while(candidates.size() < numCandidates) {
        int i = reachMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        int j = reachMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        if( reachMap(i,j) > 0 ) {
            Vantage v;
            v.x = reachMap.gridIndexToXCoord(j);
            v.y = reachMap.gridIndexToYCoord(i);
            v.z = elevationMap(i,j) + roverHeight;
            candidates.push_back(v);
        }
    }

    auto candidateMap = reachMap;

    // For every candidate, trace a ray to all view coverage probes.
    #pragma omp parallel for
    for(int ci=0; ci<candidates.size(); ++ci) {
        fmt::print("Computing Visibility {:6d}/{:6d}\n", ci, candidates.size());
        auto& candidate = candidates[ci];
        candidate.coverage.resize(probes.size());

        for(int pi = 0; pi<probes.size(); ++pi) {
            auto& p = probes[pi];
            TerrainMesh::Ray ray;
            ray.oX = candidate.x;
            ray.oY = candidate.y;
            ray.oZ = candidate.z + roverHeight;
            ray.dX = p.x-candidate.x;
            ray.dY = p.y-candidate.y;
            ray.dZ = p.z-candidate.z;
            double rayNorm = std::sqrt(ray.dX*ray.dX+ray.dY*ray.dY+ray.dZ*ray.dZ);
            ray.dX /= rayNorm; ray.dY /= rayNorm; ray.dZ /= rayNorm;
            const auto hit = tmesh.raytrace(ray);
            if( hit ) {
                double hitAngle = 180/M_PI * std::acos(-ray.dX*hit->nx-ray.dY*hit->ny-ray.dZ*hit->nz);
                double hitDist = std::sqrt((ray.oX-hit->x)*(ray.oX-hit->x) +
                                           (ray.oY-hit->y)*(ray.oY-hit->y) +
                                           (ray.oZ-hit->z)*(ray.oZ-hit->z));
                if( hitAngle < visAngle && std::abs(rayNorm-hitDist) < 0.05*rayNorm ) {
                    candidate.coverage[pi] = true;
                    candidate.totalCoverage += p.priority;
                }
            }
        }
        int i = candidateMap.yCoordToGridIndex(candidate.y);
        int j = candidateMap.xCoordToGridIndex(candidate.x);
        candidateMap(i,j) = candidate.totalCoverage;
    }
    return std::make_pair(candidates, candidateMap);
}

std::vector<Vantage> selectVantages(const std::vector<Vantage>& candidates,
                                    const std::vector<Probe>& probes,
                                    int numVantages) {

    std::vector<Vantage> vantages;
    std::set<int> taken;
    std::vector<unsigned char> visCounters(probes.size(), 0);

    // Make k selections. Choose the candidate that produces the greatest *new* coverage.
    for(int k = 0; k<numVantages; ++k) {
        fmt::print("Selecting vantages {}/{}\n", k, numVantages);

        // Assign a score to every candidate.
        std::vector<float> scores(candidates.size(), 0.0f);
        #pragma omp parallel for
        for(int ci=0; ci<candidates.size(); ++ci) {
            if( taken.contains(ci) ) { continue; }
            for(int pi=0; pi<probes.size(); ++pi) {
                if( !candidates[ci].coverage[pi] ) { continue; }

                if( visCounters[pi] == 0 ) {
                    scores[ci] += 1.2 * probes[pi].priority;
                } else if( visCounters[pi] == 1 ) {
                    scores[ci] += 1.0 * probes[pi].priority;
                } else if( visCounters[pi] == 2 ) {
                    scores[ci] += 0.01 * probes[pi].priority;
                }
            }
        }

        // Select the candidate with the highest score.
        int best_ci = 0;
        for(int ci=0; ci<candidates.size(); ++ci) {
            if( scores[ci] > scores[best_ci] ) {
                best_ci = ci;
            }
        }
        taken.insert(best_ci);
        fmt::print("Score: {}\n", scores[best_ci]);

        // Add new visibility to the visCounters.
        for(int pi=0; pi<probes.size(); ++pi) {
            if( candidates[best_ci].coverage[pi] ) {
                visCounters[pi]++;
            }
        }
        vantages.push_back(candidates[best_ci]);
    }
    return vantages;
}

TerrainMap buildCoverageMap(const TerrainMesh& tmesh, const TerrainMap& elevationMap,
                            const std::vector<Vantage>& vantages, double roverHeight, double visAngle) {

    TerrainMap coverageMap(elevationMap.width(), elevationMap.height(), elevationMap.pitch);

    for(int vi=0; vi<vantages.size(); ++vi) {
        auto& v = vantages[vi];
        fmt::print("Generating Coverage Map: {}/{}\n", vi, vantages.size());
        for(int ri=0; ri<coverageMap.cols; ++ri) {
            for(int rj=0; rj<coverageMap.rows; ++rj) {
                TerrainMesh::Ray ray;
                const double terrainX = elevationMap.gridIndexToXCoord(rj);
                const double terrainY = elevationMap.gridIndexToYCoord(ri);
                const double terrainZ = elevationMap(ri, rj);
                ray.oX = v.x;
                ray.oY = v.y;
                ray.oZ = v.z+roverHeight;
                ray.dX = terrainX-ray.oX;
                ray.dY = terrainY-ray.oY;
                ray.dZ = terrainZ-ray.oZ;
                double rayNorm = std::sqrt(ray.dX*ray.dX + ray.dY*ray.dY + ray.dZ*ray.dZ);
                ray.dX /= rayNorm;
                ray.dY /= rayNorm;
                ray.dZ /= rayNorm;
                const auto hit = tmesh.raytrace(ray);
                if( hit ) {
                    double hitAngle = 180/M_PI * std::acos(-ray.dX*hit->nx-ray.dY*hit->ny-ray.dZ*hit->nz);
                    double hitDist = std::sqrt((ray.oX-hit->x)*(ray.oX-hit->x) +
                                               (ray.oY-hit->y)*(ray.oY-hit->y) +
                                               (ray.oZ-hit->z)*(ray.oZ-hit->z));
                    if( hitAngle < visAngle && hitDist < 500 && std::abs(rayNorm-hitDist) < rayNorm*0.05 ) { 
                        coverageMap(ri, rj)++;
                    }
                } 
            }
        }
    }
    for(const auto& v : vantages) {
        int i = coverageMap.yCoordToGridIndex(v.y);
        int j = coverageMap.xCoordToGridIndex(v.x);
        coverageMap(i,j) = vantages.size();
    }
    return coverageMap;
}

int main(int argc, char* argv[]) {

    // Configure the planner.
    const auto config = parseCommandLine(argc, argv);
    if( !config ) { return 0; }

    // Read the terrain mesh.
    TerrainMesh tmesh(config->meshfile);

    // Construct elevation, slope, and priority maps.
    auto [elevationMap, slopeMap, priorityMap] = buildTerrainMaps(tmesh, config->mapPitch);
    elevationMap.savePFM("elevation.pfm");
    slopeMap.savePFM("slope.pfm");
    priorityMap.savePFM("priority.pfm");

    // Construct lander communications map.
    TerrainMap commsMap = buildCommsMap(tmesh, elevationMap,
                                            config->landingSiteX, config->landingSiteY,
                                            config->landerHeight, config->roverHeight);
    commsMap.savePFM("comms.pfm");

    // Map low-slope, communicable terrain.
    const auto safeMap = buildSafeMap(commsMap, slopeMap, config->roverMaxSlope);
    safeMap.savePFM("safe.pfm");

    // Flood-fill reachable safe terrain.
    TerrainMap reachMap = buildReachabilityMap(safeMap, config->landingSiteX, config->landingSiteY, 13);
    reachMap.savePFM("reach.pfm");

    // Generate visibility probes.
    auto [probes, probeMap] = generateVisibilityProbes(priorityMap, elevationMap, config->numProbes, config->roverHeight);
    probeMap.savePFM("probes.pfm");

    // Generate candidate vantages.
    auto [candidates, candidateMap] = generateVantageCandidates(tmesh, reachMap, elevationMap, probes,
                                                                config->numCandidates, config->roverHeight, config->visAngle);
    candidateMap.savePFM("candidates.pfm");

    // Select the best vantages from all of the candidates.
    auto vantages = selectVantages(candidates, probes, config->numVantages);
    auto vantageMap = reachMap;
    for(const auto& v : vantages) {
        int j = vantageMap.xCoordToGridIndex(v.x);
        int i = vantageMap.yCoordToGridIndex(v.y);
        vantageMap(i,j) = v.totalCoverage;
    }
    vantageMap.savePFM("vantages.pfm");

    auto coverageMap = buildCoverageMap(tmesh, elevationMap, vantages, config->roverHeight, config->visAngle);
    coverageMap.savePFM("coverage.pfm");

    for(int vi=0; vi<vantages.size(); ++vi) {
        std::vector<Vantage> tmp;
        tmp.push_back(vantages[vi]);
        auto coverageMap = buildCoverageMap(tmesh, elevationMap, tmp, config->roverHeight, config->visAngle);
        int j = vantageMap.xCoordToGridIndex(vantages[vi].x);
        int i = vantageMap.yCoordToGridIndex(vantages[vi].y);
        coverageMap(i,j) = 2;
        coverageMap.savePFM(fmt::format("coverage_{:02}.pfm",vi));
    }

    // Compute paths between all pairs of k vantages plus the landing site.

    // Compute exploration route.

    // Draw pretty pictures.

    return 0;
}
