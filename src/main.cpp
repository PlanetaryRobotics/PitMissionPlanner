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
    double landingSiteX = 825;
    double landingSiteY = 1000;
    double landerHeight = 3.0;
    double roverHeight = 1.0;
    double roverMaxSlope = 15.0;
    double numCandidates = 100000;
    double visAngle = 55;
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
        fmt::print("\tnc: the number of candidate vantages to evaluate.\n");
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
    cmdl("nc", cfg.numCandidates) >> cfg.numCandidates;
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

TerrainMap<int> buildReachabilityMap(const TerrainMap<int>& safeMap,
                                      double siteX, double siteY,
                                      double roverMaxSlope) {
    TerrainMap<int> reachMap(safeMap.width(), safeMap.height(), safeMap.pitch);

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

TerrainMap<int> buildCommsMap(const TerrainMesh& tmesh,
                               const TerrainMap<float>& elevationMap,
                               double siteX, double siteY,
                               double landerHeight, double roverHeight) {

    TerrainMap<int> commsMap(elevationMap.width(),
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

TerrainMap<int> buildSafeMap(const TerrainMap<int>& commsMap, const TerrainMap<float>& slopeMap, double roverMaxSlope) {
    TerrainMap<int> safeMap = commsMap;
    for(int i=0; i<safeMap.rows; ++i) {
        for(int j=0; j<safeMap.cols; ++j) {
            safeMap(i,j) = (commsMap(i,j) > 0 && slopeMap(i,j) <= roverMaxSlope) ? 1 : 0;
        }
    }
    return safeMap;
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

    // Generate candidate vantages.
    struct Vantage {
        double x, y, z;
        std::vector<bool> coverage;
        double totalCoverage = 0;
    };
    std::vector<Vantage> candidates;

    auto candidateMap = reachMap;
    while(candidates.size() < config->numCandidates) {
        int i = reachMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        int j = reachMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        if( reachMap(i,j) > 0 ) {
            candidateMap(i,j) = 2;
            Vantage v;
            v.x = reachMap.gridIndexToXCoord(j);
            v.y = reachMap.gridIndexToYCoord(i);
            v.z = elevationMap(i,j) + config->roverHeight;
            candidates.push_back(v);
        }
    }
    candidateMap.savePFM("candidates.pfm");

    // Compute visibility coverage for each candidate.
    // Generate visibility probes.
    struct Probe {
        double x,y,z;
        double priority;
    };
    int numProbes = 10000;
    std::vector<Probe> probes;
    auto probeMap = priorityMap;
    while(probes.size() < numProbes) {
        int i = reachMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        int j = reachMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        if( priorityMap(i,j) > 0 ) {
            probeMap(i,j) = 2;
            Probe p;
            p.x = reachMap.gridIndexToXCoord(j);
            p.y = reachMap.gridIndexToYCoord(i);
            p.z = elevationMap(i,j) + config->roverHeight;
            p.priority = priorityMap(i,j);
            probes.push_back(p);
        }
    }
    probeMap.savePFM("probes.pfm");

    // For every candidate, trace a ray to all view coverage probes.
    #pragma omp parallel for
    for(int ci=0; ci<candidates.size(); ++ci) {
        fmt::print("Computing Visibility {}/{}\n", ci, candidates.size());
        auto& candidate = candidates[ci];
        candidate.coverage.resize(probes.size());

        TerrainMesh::Ray ray;
        ray.oX = candidate.x;
        ray.oY = candidate.y;
        ray.oZ = candidate.z + config->roverHeight;
        for(int pi = 0; pi<probes.size(); ++pi) {
            auto& p = probes[pi];
            ray.dX = p.x-candidate.x;
            ray.dY = p.y-candidate.y;
            ray.dZ = p.z-candidate.z;
            double rayNorm = std::sqrt(ray.dX*ray.dX+ray.dY*ray.dY+ray.dZ*ray.dZ);
            ray.dX /= rayNorm; ray.dY /= rayNorm; ray.dZ /= rayNorm;
            const auto hit = tmesh.raytrace(ray);
            if( hit ) {
                double hitAngle = 180/M_PI * std::acos(-ray.dX*hit->nx-ray.dY*hit->ny-ray.dZ*hit->nz);
                if( hitAngle < config->visAngle ) {
                    candidate.coverage[pi] = true;
                    candidate.totalCoverage += p.priority;
                }
            }
        }
        int i = candidateMap.yCoordToGridIndex(candidate.y);
        int j = candidateMap.xCoordToGridIndex(candidate.x);
        candidateMap(i,j) = candidate.totalCoverage;
    }
    candidateMap.savePFM("candidates_coverage.pfm");

    // Select k vantages with maximum weighted coverage.
    int numVantages = 30;
    std::vector<unsigned char> visCounters(probes.size());
    std::vector<Vantage> finalVantages;

    candidateMap = reachMap;
    std::set<int> taken;
    // Make k selections. Choose the candidate that produces the greatest *new* coverage.
    for(int k = 0; k<numVantages; ++k) {
        int best_ci = 0;
        for(int ci=0; ci<candidates.size(); ++ci) {
            if( candidates[ci].totalCoverage > candidates[best_ci].totalCoverage && !taken.contains(ci) ) {
                best_ci = ci;
            }
        }
        taken.insert(best_ci);
        Vantage newVantage;
        newVantage.x = candidates[best_ci].x;
        newVantage.y = candidates[best_ci].y;
        newVantage.z = candidates[best_ci].z;
        finalVantages.push_back(newVantage);
        int i = candidateMap.yCoordToGridIndex(newVantage.y);
        int j = candidateMap.xCoordToGridIndex(newVantage.x);
        candidateMap(i,j) = candidates[best_ci].totalCoverage;
    }
    candidateMap.savePFM("vantages.pfm");

    // Compute paths between all pairs of k vantages plus the landing site.

    // Compute exploration route.

    // Draw pretty pictures.

    return 0;
}
