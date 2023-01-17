#include "terrainmesh.h"
#include "terrainmap.h"
#include "argh.h"
#include "flat_hash_map.hpp"
#include "priority_queue.h"
#include <Eigen/Dense>
#include <fmt/format.h>
#include <optional>
#include <unordered_set>
#include <iostream>


struct PlannerConfiguration {
    std::string meshfile   = "../meshes/lmp.ply";
    std::string outputDir  = "./";
    double mapPitch        = 1.0;       // meters
    double landingSiteX    = 750;       // meters
    double landingSiteY    = 400;       // meters
    double landerHeight    = 3.0;       // meters
    double roverHeight     = 1.0;       // meters
    double roverMaxSlope   = 20.0;      // degrees
    double roverSpeed      = 0.01;      // m/s
    int    numProbes       = 10000;     // 
    int    numCandidates   = 10000;     // 
    int    numVantages     = 15;        // 
    double visAngle        = 55;        // degrees
};

std::optional<PlannerConfiguration> parseCommandLine(int argc, char* argv[]) {
    PlannerConfiguration cfg;

    argh::parser cmdl(argc, argv);

    if( cmdl[{ "-h", "--help" }] ) {
        fmt::print("usage: planranger [-p][-x][-y][-rs][-rh][-lh] meshfile outdir\n");
        fmt::print("\tmeshfile: a .ply format map of terrain surrounding a lunar pit.\n");
        fmt::print("\toutdir: a directory in which to place the output.\n");
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
    cmdl(2, cfg.outputDir) >> cfg.outputDir;
    cfg.outputDir = cfg.outputDir + "/";

    cmdl("p", cfg.mapPitch) >> cfg.mapPitch;
    cmdl("x", cfg.landingSiteX) >> cfg.landingSiteX;
    cmdl("y", cfg.landingSiteY) >> cfg.landingSiteY;
    cmdl("rs", cfg.roverMaxSlope) >> cfg.roverMaxSlope;
    cmdl("rh", cfg.roverHeight) >> cfg.roverHeight;
    cmdl("lh", cfg.landerHeight) >> cfg.landerHeight;
    cmdl("np", cfg.numProbes) >> cfg.numProbes;
    cmdl("nc", cfg.numCandidates) >> cfg.numCandidates;
    cmdl("nv", cfg.numVantages) >> cfg.numVantages;
    cmdl("va", cfg.visAngle) >> cfg.visAngle;
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
                if (std::isnan(slopeMap(i,j))) { slopeMap(i,j) = 0.0f; }
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

    for(int ri=0; ri<commsMap.rows; ++ri) {
        for(int rj=0; rj<commsMap.cols; ++rj) {
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
                if( hitAngle < visAngle && hitDist < 500 && std::abs(rayNorm-hitDist) < 0.05*rayNorm ) {
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
                                    int numVantages, double minSeparation=10.0) {

    std::vector<Vantage> vantages;
    std::unordered_set<int> taken;
    std::vector<unsigned char> visCounters(probes.size(), 0);

    // Make k selections. Choose the candidate that produces the greatest *new* coverage.
    for(int k = 0; k<numVantages; ++k) {
        fmt::print("Selecting vantages {}/{}\n", k, numVantages);

        // Assign a score to every candidate.
        std::vector<float> scores(candidates.size(), 0.0f);
        #pragma omp parallel for
        for(int ci=0; ci<candidates.size(); ++ci) {
            if( taken.contains(ci) ) { continue; }
            const auto& c = candidates[ci];

            // If this candidate is too close to one
            // that has already been selected, skip it.
            bool tooClose = false;
            for(const auto& ti : taken) {
                const auto& t = candidates[ti];
                double d2 = (t.x-c.x)*(t.x-c.x)+(t.y-c.y)*(t.y-c.y)+(t.z-c.z);
                if( d2 < minSeparation*minSeparation ) {
                    tooClose = true;
                    break;
                }
            }
            if( tooClose ) { continue; }

            for(int pi=0; pi<probes.size(); ++pi) {
                if( !c.coverage[pi] ) { continue; }

                if( visCounters[pi] == 0 ) {
                    scores[ci] += 1.0 * probes[pi].priority;
                } else if( visCounters[pi] == 1 ) {
                    scores[ci] += 0.8 * probes[pi].priority;
                } else if( visCounters[pi] == 2 ) {
                    scores[ci] += 0.05 * probes[pi].priority;
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
        #pragma omp parallel for
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
        int R = 2.0 / coverageMap.pitch;
        for(int ii=-R; ii<=R; ++ii) {
            for(int jj=-R; jj<=R; ++jj) {
                if( ii*ii + jj*jj <= R*R ) {
                    coverageMap(i+ii,j+jj) = 2*vantages.size();
                }
            }
        }
    }
    return coverageMap;
}

struct Path {
    struct Point {
        int i = 0;
        int j = 0;
        bool operator==(const Point& rhs) const {
            return rhs.i==i && rhs.j==j;
        }
    };
    std::vector<Point> waypoints;

    double cost = 0.0;
    double dist = 0.0;
};

Path reverse(const Path& path) {
    Path rPath = path;
    std::reverse(rPath.waypoints.begin(), rPath.waypoints.end());
    return rPath;
}

Path pathplan(const TerrainMap& slopeMap, const TerrainMap& reachMap,
              const Path::Point& start, const Path::Point& goal) {

    int ROWS = slopeMap.rows;
    int COLS = slopeMap.cols;

    // Construct planning datastructures.
    auto pointHash = [COLS](const Path::Point& p) { return p.i*COLS + p.j; };
    ska::flat_hash_set<Path::Point, decltype(pointHash)> closed(10, pointHash);
    PriorityQueue<Path::Point, double> open;

    struct NodeData {
        Path::Point pred;
        double gscore;
    };
    ska::flat_hash_map<Path::Point, NodeData, decltype(pointHash)> nodeMap(10, pointHash);

    // Define the planning heuristic.
    auto getHeuristic = [](const Path::Point& a, const Path::Point& b) -> double {
        // Branchless octile distance
        int di = std::abs(a.i-b.i);
        int dj = std::abs(a.j-b.j);
        constexpr double diagonal = std::sqrt(2);
        constexpr double twoCardinalMinusDiagonal = 2-diagonal;
        return (twoCardinalMinusDiagonal*abs(di-dj) + diagonal*(di+dj)) / 2;
    };

    // Generate the successors of a node.
    auto getSuccessors = [ROWS, COLS, &reachMap](const Path::Point& p) -> std::vector<Path::Point> {
        std::vector<Path::Point> succs; succs.reserve(8);
        if( p.i+1 < ROWS && reachMap(p.i+1, p.j) ) { succs.push_back( Path::Point { .i = p.i+1, .j = p.j } ); }
        if( p.i-1 >= 0   && reachMap(p.i+1, p.j) ) { succs.push_back( Path::Point { .i = p.i-1, .j = p.j } ); }
        if( p.j+1 < COLS && reachMap(p.i+1, p.j) ) { succs.push_back( Path::Point { .i = p.i, .j = p.j+1 } ); }
        if( p.j-1 >= 0   && reachMap(p.i+1, p.j) ) { succs.push_back( Path::Point { .i = p.i, .j = p.j-1 } ); }

        if( p.i+1 < ROWS && p.j+1 < COLS && reachMap(p.i+1, p.j+1) ) {
            succs.push_back( Path::Point { .i = p.i+1, .j = p.j+1 } );
        }
        if( p.i-1 >= 0 && p.j-1 >= 0 && reachMap(p.i-1, p.j-1) ) {
            succs.push_back( Path::Point { .i = p.i-1, .j = p.j-1 } );
        }
        if( p.i+1 < ROWS && p.j-1 >= 0 && reachMap(p.i+1, p.j-1) ) {
            succs.push_back( Path::Point { .i = p.i+1, .j = p.j-1 } );
        }
        if( p.i-1 >= 0 && p.j+1 < COLS && reachMap(p.i-1, p.j+1) ) {
            succs.push_back( Path::Point { .i = p.i-1, .j = p.j+1 } );
        }
        return succs; 
    };

    // Compute the cost to travel from one point to another.
    auto getCost = [&slopeMap, &getHeuristic](const Path::Point& a, const Path::Point& b) -> double {
        const double slope = slopeMap(b.i, b.j);
        const double octileDistance = slopeMap.pitch * getHeuristic(a, b);
        return 8 * slope / 90.0 + octileDistance;
    };

    // Initialize the search.
    double h = getHeuristic(start, goal);
    open.insert(start, h);

    nodeMap[start].pred = start;
    nodeMap[start].gscore = 0.0;

    // Perform the search.
    bool success = false;
    while( true )
    {
        // The goal is unreachable.
        if( open.empty() ) { success = false; break; }

        const Path::Point currPoint = open.top();
        open.pop();

        // If you have expanded this state before, skip it.
        if( closed.find(currPoint) != closed.end() ) { continue; }

        // Otherwise, add it to closed and expand it.
        closed.insert(currPoint);

        // If this state is the goal, quit searching!
        if( currPoint == goal ) {
            success = true;
            break;
        }

        // Since we can't end the search yet, let's look
        // up some info on the current point and keep working.
        const double currG = nodeMap[currPoint].gscore;

        // Examine each successor of currPoint.
        const auto succs = getSuccessors(currPoint);
        for(const auto& succ : succs) {
            const double succCost = getCost(currPoint, succ);
            const double tentativeG = currG + succCost;

            if( nodeMap.find(succ) == nodeMap.end() ||
                tentativeG < nodeMap.at(succ).gscore ) {
                double h = getHeuristic(succ, goal);
                double f = tentativeG + h;
                open.insert(succ, f);

                nodeMap[succ].pred = currPoint;
                nodeMap[succ].gscore = tentativeG;
            }
        }
    }

    Path path;
    if( !success ) { 
        path.cost = -1;
        path.dist = -1;
        return path;
    }

    Path::Point state = goal;
    while( !(state == start) ) {
        path.waypoints.push_back(state);
        Path::Point pred = nodeMap[state].pred;
        path.dist += std::sqrt((pred.i-state.i)*(pred.i-state.i)+(pred.j-state.j)*(pred.j-state.j));
        state = pred;
    }
    path.waypoints.push_back(start);
    path.cost = nodeMap[goal].gscore;

    path = reverse(path);

    return path;
}

std::vector<Vantage> routeplan(const std::vector<Vantage> allSites, const Eigen::MatrixXd& dists) {
    std::vector<int> routeIndices;
    std::vector<bool> visited(allSites.size(), false);
    int visitedCount = 1;

    int currentIdx = 0;
    visited[currentIdx] = true;
    routeIndices.push_back(currentIdx);
    double routeLength = 0;

    while(visitedCount < allSites.size()) {
        int closestIdx = 0;
        double closestDist = std::numeric_limits<double>::infinity();
        for(int i=0; i<dists.rows(); ++i) {
            const auto dist = dists(currentIdx, i);
            if( !visited[i] && dist >= 0 && dist < closestDist ) {
                closestIdx = i;
                closestDist = dist;
            }
        }
        currentIdx = closestIdx;
        routeIndices.push_back(closestIdx);
        visited[closestIdx] = true;
        routeLength += closestDist;
        visitedCount++;
    }

    // Two-Opt Route Improvement
    int iterations = 0;
    bool improved = false;
    do {
        fmt::print("Improving route [{}] ...\n", iterations);
        fmt::print("Route Length: {}\n", routeLength);
        improved = false;
        for(int i=1; i<routeIndices.size()-1; ++i) {

            // Don't rewire the segment 0<->1 because that segment
            // connects the landing site to the closest vantage.
            if( i== 1 ) { continue; }

            // NOTES (Jordan):
            // There is a tension here.
            // On the one hand, we want a short total path.
            // On the other hand, we want to visit vantages *early*.
            // I think we should apply a discount factor that values
            // vantages less if they are visited later in time.
            // Then we want a short route that captures the most
            // discounted coverage. I guess the discount factor will be
            // a configurable parameter.

            for(int j=i+1; j<routeIndices.size(); ++j) {
                std::vector<int> newIndices;

                // Take route[0] to route[i-1] and add them to the new route.
                for(int k=0; k<=i-1; ++k) {
                    newIndices.push_back(routeIndices[k]);
                }
                // Take route[i] to route[j] and add them in reverse order.
                for(int k=j; k>=i; --k) {
                    newIndices.push_back(routeIndices[k]);
                }
                // Take route[j+1] to the end and add them to the new route.
                for(int k=j+1; k<routeIndices.size(); ++k) {
                    newIndices.push_back(routeIndices[k]);
                }
                // Calculate the length of the new route.
                double dist = 0;
                for(int i=0; i<newIndices.size()-1; ++i) {
                    dist += dists(newIndices[i], newIndices[i+1]);
                }
                // This route is shorter! Keep it.
                if( dist < routeLength ) {
                    fmt::print("[{}] Swap {:2}<->{:2} Length: {}\n", iterations, i,j, dist);
                    fmt::print("Found a better answer!\n");
                    routeLength = dist;
                    routeIndices = newIndices;
                    improved = true;
                }
            }
        }
    } while( improved && iterations++ < 100 );

    std::vector<Vantage> route;
    for(const auto& ri : routeIndices) { route.push_back(allSites[ri]); }

    return route;
}

// Sort vantages by their angle relative to (siteX, siteY).
std::vector<Vantage> sortCCW(const std::vector<Vantage> vantages, double siteX, double siteY) {
    std::vector<Vantage> sorted = vantages;
    auto angle = [siteX, siteY](const Vantage& v0, const Vantage& v1) {
        double aX = v0.x-siteX;
        double aY = v0.y-siteY;
        double bX = v1.x-siteX;
        double bY = v1.y-siteY;
        return std::atan2(aY, aX) < std::atan2(bY, bX);
    };
    std::sort(sorted.begin(), sorted.end(), angle);
    return sorted;
}

int main(int argc, char* argv[]) {

    // Configure the planner.
    const auto config = parseCommandLine(argc, argv);
    if( !config ) { return 0; }

    // Read the terrain mesh.
    TerrainMesh tmesh(config->meshfile);

    // Construct elevation, slope, and priority maps.
    auto [elevationMap, slopeMap, priorityMap] = buildTerrainMaps(tmesh, config->mapPitch);

    // Compute landing site coordinates.
    double landingSiteX = config->landingSiteX;
    double landingSiteY = config->landingSiteY;
    double landingSiteZ = elevationMap.atXY(landingSiteX, landingSiteY);

    int landingSiteI = elevationMap.yCoordToGridIndex(landingSiteY);
    int landingSiteJ = elevationMap.xCoordToGridIndex(landingSiteX);

    if( landingSiteI < 0 || landingSiteI >= elevationMap.cols ||
        landingSiteJ < 0 || landingSiteJ >= elevationMap.rows ) {
        throw std::runtime_error(
            fmt::format("Landing site at ({}, {}) is outside of map boundaries.",
                        landingSiteX, landingSiteY));
    }

    // Construct lander communications map.
    TerrainMap commsMap = buildCommsMap(tmesh, elevationMap,
                                            landingSiteX, landingSiteY,
                                            config->landerHeight, config->roverHeight);

    // Map low-slope, communicable terrain.
    const auto safeMap = buildSafeMap(commsMap, slopeMap, config->roverMaxSlope);

    // Flood-fill reachable safe terrain.
    TerrainMap reachMap = buildReachabilityMap(safeMap, landingSiteX, landingSiteY, 13);

    // Generate visibility probes.
    auto [probes, probeMap] = generateVisibilityProbes(priorityMap, elevationMap, config->numProbes, config->roverHeight);

    // Generate candidate vantages.
    auto [candidates, candidateMap] = generateVantageCandidates(tmesh, reachMap, elevationMap, probes,
                                                                config->numCandidates, config->roverHeight, config->visAngle);

    // Select the best vantages from all of the candidates.
    auto vantages = selectVantages(candidates, probes, config->numVantages);

    // Construct the vantageMap.
    auto vantageMap = slopeMap;
    {
        double maxCoverage = 0;
        for(const auto& v : vantages) { maxCoverage = std::max(maxCoverage, v.totalCoverage); }
        for(const auto& v : vantages) {
            int j = vantageMap.xCoordToGridIndex(v.x);
            int i = vantageMap.yCoordToGridIndex(v.y);
            double markerValue = 120 + 10 * v.totalCoverage / maxCoverage;
            int R = 2.0 / config->mapPitch;
            for(int ii=-R; ii<=R; ++ii) {
                for(int jj=-R; jj<=R; ++jj) {
                    if( ii*ii + jj*jj <= R*R ) {
                        vantageMap(i+ii,j+jj) = markerValue;
                    }
                }
            }
        }
    }

    // Map combined coverage from all vantages.
    auto coverageMap = buildCoverageMap(tmesh, elevationMap, vantages, config->roverHeight, config->visAngle);

    // Save landing site and vantages to xyz file.
    {
        std::ofstream file;
        file.open(config->outputDir+"sites.xyz");
        file << fmt::format("{} {} {}\n", landingSiteX, landingSiteY, landingSiteZ);
        for(const auto& v : vantages) {
            file << fmt::format("{} {} {}\n", v.x, v.y, v.z + config->roverHeight);
        }
        file.close();
    }

    // Save maps.
    {
        {
            auto map = elevationMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"elevation.exr");
        }
        {
            auto map = slopeMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"slope.exr");
        }
        {
            auto map = priorityMap;
            map.drawCircle(landingSiteX, landingSiteY, 10, 4.0);
            map.saveEXR(config->outputDir+"priority.exr");
        }
        {
            auto map = reachMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"reach.exr");
        }
        {
            auto map = commsMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"comms.exr");
        }
        {
            auto map = safeMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"safe.exr");
        }
        {
            auto map = probeMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"probes.exr");
        }
        {
            auto map = candidateMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"candidates.exr");
        }
        {
            auto map = vantageMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"vantages.exr");
        }
        {
            auto map = coverageMap;
            map.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
            map.saveEXR(config->outputDir+"coverage.exr");
        }
        for(int vi=0; vi<vantages.size(); ++vi) {
            std::vector<Vantage> tmp;
            tmp.push_back(vantages[vi]);
            auto coverageMap = buildCoverageMap(tmesh, elevationMap, tmp, config->roverHeight, config->visAngle);
            int j = vantageMap.xCoordToGridIndex(vantages[vi].x);
            int i = vantageMap.yCoordToGridIndex(vantages[vi].y);
            coverageMap.drawCircle(vantages[vi].x, vantages[vi].y, vantages.size()+1, 3.0);
            coverageMap.drawCircle(landingSiteX, landingSiteY, vantages.size()+10, 4.0);
            coverageMap.saveEXR(config->outputDir+fmt::format("coverage_{:02}.exr",vi));
        }
    }

    // Compute paths between all pairs of k vantages plus the landing site.
    Eigen::MatrixXd costs(vantages.size()+1, vantages.size()+1); costs.fill(-1);
    Eigen::MatrixXd dists(vantages.size()+1, vantages.size()+1); dists.fill(-1);

    using PointPair = std::pair<Path::Point, Path::Point>;
    auto hashPointPair = [&elevationMap](const PointPair& pp) {
        int ROWS = elevationMap.rows;
        int COLS = elevationMap.cols;
        int h0 = pp.first.i*COLS + pp.first.j;
        int h1 = pp.second.i*COLS + pp.second.j;
        return h0*ROWS*COLS + h1;
    };
    ska::flat_hash_map<PointPair, Path, decltype(hashPointPair)> pathLookup(10, hashPointPair);
    Vantage landingSite {.x = landingSiteX, .y = landingSiteY, .z = landingSiteZ, .totalCoverage = 0.0};

    // Sort vantages counterclockwise around their centroid.
    double centroidX = 0; double centroidY = 0;
    for(const auto& v : vantages) { centroidX += v.x; centroidY += v.y; }
    centroidX /= vantages.size(); centroidY /= vantages.size();
    vantages = sortCCW(vantages, centroidX, centroidY);

    std::vector<Vantage> allSites;
    allSites.push_back(landingSite);
    allSites.insert(std::end(allSites), std::begin(vantages), std::end(vantages));

    // Generate all combinations of two indices into the allSites vector.
    std::vector<std::pair<int,int>> allPairIndices;
    for(int a=0; a<allSites.size(); ++a) {
        for(int b=a+1; b<allSites.size(); ++b) {
            allPairIndices.push_back(std::make_pair(a,b));
        }
    }

    #pragma omp parallel for
    for(int i=0; i<allPairIndices.size(); ++i) {
        auto [a, b] = allPairIndices[i];
        fmt::print("Planning path from {:2} to {:2} [{}/{}]\n", a, b, i, allPairIndices.size());

        Path::Point start, goal;
        start.i = slopeMap.yCoordToGridIndex(allSites[a].y);
        start.j = slopeMap.xCoordToGridIndex(allSites[a].x);
        goal.i = slopeMap.yCoordToGridIndex(allSites[b].y);
        goal.j = slopeMap.xCoordToGridIndex(allSites[b].x);

        const auto path = pathplan(slopeMap, reachMap, start, goal);

        pathLookup[std::make_pair(start, goal)] = path;
        pathLookup[std::make_pair(goal, start)] = reverse(path);

        costs(a,b) = path.cost;
        costs(b,a) = path.cost;
        dists(a,b) = path.dist;
        dists(b,a) = path.dist;
    }

    // Compute exploration route.
    std::vector<Vantage> route = routeplan(allSites, dists);

    // Draw the final route!
    TerrainMap routeMap = vantageMap;
    routeMap.drawCircle(landingSiteX, landingSiteY, 100, 4.0);
    for(int i=0; i<route.size()-1; ++i) {
        Vantage v0 = route[i]; 
        Vantage v1 = route[i+1]; 
        Path::Point p0, p1;
        p0.i = routeMap.yCoordToGridIndex(v0.y);
        p0.j = routeMap.xCoordToGridIndex(v0.x);
        p1.i = routeMap.yCoordToGridIndex(v1.y);
        p1.j = routeMap.xCoordToGridIndex(v1.x);
        Path path = pathLookup[std::make_pair(p0, p1)];
        for(const auto& p : path.waypoints) {
            routeMap(p.i, p.j) = 100;
        }
    }
    routeMap.saveEXR(config->outputDir+"route.exr"); 

    // Save route to xyz file.
    {
        std::ofstream file;
        file.open(config->outputDir+"route.xyz");
        for(int i=0; i<route.size()-1; ++i) {
            Vantage v0 = route[i]; 
            Vantage v1 = route[i+1]; 
            Path::Point p0, p1;
            p0.i = routeMap.yCoordToGridIndex(v0.y);
            p0.j = routeMap.xCoordToGridIndex(v0.x);
            p1.i = routeMap.yCoordToGridIndex(v1.y);
            p1.j = routeMap.xCoordToGridIndex(v1.x);
            Path path = pathLookup[std::make_pair(p0, p1)];
            for(const auto& p : path.waypoints) {
                Vantage v;
                v.x = elevationMap.gridIndexToXCoord(p.j);
                v.y = elevationMap.gridIndexToXCoord(p.i);
                v.z = elevationMap(p.i, p.j) + config->roverHeight;
                file << fmt::format("{} {} {}\n", v.x, v.y, v.z);
            }
        }
        file.close();
    }

    return 0;
}
