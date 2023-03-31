#include "image.h"
#include "terrainmesh.h"
#include "terrainmap.h"
#include "argh.h"
#include "flat_hash_map.hpp"
#include "toml.hpp"
#include "priority_queue.h"
#include "config.h"
#include "path.h"
#include "vantage.h"
#include "maps.h"
#include <Eigen/Dense>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <array>
#include <algorithm>
#include <optional>
#include <unordered_set>
#include <iostream>
#include <filesystem>
#include <numeric>
#include "ortools/constraint_solver/routing.h"
#include "ortools/constraint_solver/routing_enums.pb.h"
#include "ortools/constraint_solver/routing_index_manager.h"
#include "ortools/constraint_solver/routing_parameters.h"

#include <chrono>

// Global configuration structure
PlannerConfiguration config;




std::pair<std::vector<Path>, int>
multipathplan(const SlopeAtlas& slopeAtlas,
              const TerrainMapFloat& commsMap,
              const Path::State& start,
              const std::vector<Path::State>& goals,
              const int threadID = 0) {
    int ROWS = slopeAtlas.rows();
    int COLS = slopeAtlas.cols();

    // Shove all the goals into a set.
    ska::flat_hash_set<Path::State> goalSet;
    for(const auto& g : goals) { goalSet.insert(g); }

    // Construct planning datastructures.
    ska::flat_hash_set<Path::State> closed;
    PriorityQueue<Path::State, double> open;

    struct NodeData {
        Path::State pred;
        double gscore;
    };
    ska::flat_hash_map<Path::State, NodeData> nodeMap;

    // Define the planning heuristic.
    auto getHeuristic = [&slopeAtlas](const Path::State& a, const Path::State& b) -> double {
        const auto goalDir = directionFromTo(a, b);
        const double turn = 0.25 * std::abs(std::abs(static_cast<int>(goalDir)-static_cast<int>(a.d))-4);
        const double dist = slopeAtlas.pitch() * octileDistance(a, b);
        return config.heuristicMultiplier * (config.distCostMultiplier*dist + config.turnCostMultiplier*turn);
    };
    auto getMinHeuristic = [&getHeuristic](const Path::State& a, const std::vector<Path::State>& goals) -> double {
        double minH = std::numeric_limits<double>::infinity();
        for(const auto& g : goals) {
            const auto newH = getHeuristic(a, g);
            minH = std::min<double>(newH, minH);
        }
        return minH;
    };

    // Compute the cost to travel from one point to another.
    auto getCost = [&slopeAtlas](const Path::State& a, const Path::State& b) -> double {
        const double slope = slopeAtlas.absolute(b.i, b.j);
        const double dist = slopeAtlas.pitch() * octileDistance(a, b);
        const double turn = 0.25 * std::abs(std::abs(static_cast<int>(b.d)-static_cast<int>(a.d))-4);

        // Slightly penalize driving backward.
        double back = ( dist > 0.0 && directionFromTo(a,b) != a.d ) ? 1.0 : 0.0;

        return config.slopeCostMultiplier * (slope/90.0) +
               config.distCostMultiplier * dist +
               config.turnCostMultiplier * (turn + back);
    };

    // Initialize the search.
    double h = getMinHeuristic(start, goals);
    open.insert(start, h);

    nodeMap[start].pred = start;
    nodeMap[start].gscore = 0.0;

    // Perform the search.
    int expansions = 0;
    while( true )
    {
        // The goal is unreachable.
        if( open.empty() ) {
            fmt::print("[Thread {:2}] Failed! Found only {}/{} goals.\n", threadID, goals.size()-goalSet.size(), goals.size());
            break;
        }

        const Path::State currState = open.top();
        open.pop();

        // Otherwise, add it to closed and expand it.
        closed.insert(currState);
        expansions++;

        // If this state is a goal, remove it from the goal set!
        auto it = goalSet.find(currState);
        if( it != goalSet.end() ) {
            goalSet.erase(it);
            if( goalSet.empty() ) {
                fmt::print("[Thread {:2}] Finished! Found {}/{} goals.\n", threadID, goals.size(), goals.size());
                break;
            } else {
                fmt::print("[Thread {:2}] Found {}/{} goals.\n", threadID, goals.size()-goalSet.size(), goals.size());
            }
        }

        // Since we can't end the search yet, let's look
        // up some info on the current point and keep working.
        const double currG = nodeMap.at(currState).gscore;

        // Examine each successor of currState.
        const auto succs = getSuccessors(currState, slopeAtlas, commsMap);
        for(const auto& succ : succs) {
            // If you have expanded this state before, skip it.
            if( closed.find(succ) != closed.end() ) { continue; }

            const double succCost = getCost(currState, succ);
            const double tentativeG = currG + succCost;

            if( nodeMap.find(succ) == nodeMap.end() ||
                tentativeG < nodeMap.at(succ).gscore ) {
                double h = getMinHeuristic(succ, goals);
                double f = tentativeG + h;
                open.insert(succ, f);

                nodeMap[succ].pred = currState;
                nodeMap[succ].gscore = tentativeG;
            }
        }
    }

    std::vector<Path> paths;

    // Backtrack paths from all goals.
    for(const auto& goal : goals) {
        Path path;
        if( nodeMap.find(goal) == nodeMap.end() ) {
            paths.push_back(path);
            continue;
        }
        Path::State state = goal;
        while( !(state == start) ) {
            path.states.push_back(state);
            Path::State pred = nodeMap[state].pred;
            path.dist += slopeAtlas.pitch() * octileDistance(pred, state);
            state = pred;
        }
        path.states.push_back(start);
        path.cost = nodeMap[goal].gscore;

        path = reverse(path);
        paths.push_back(path);
    }

    return std::make_pair(paths, expansions);
}

std::pair<std::vector<int>,int64_t> solveTSP(const Eigen::MatrixXd& costs, int depot=0) {
    using namespace operations_research;

    struct DataModel {
        std::vector<std::vector<int64_t>> distance_matrix;
        const int num_vehicles = 1;
        RoutingIndexManager::NodeIndex depotNode;
    };

    DataModel data;
    data.depotNode = RoutingIndexManager::NodeIndex{depot};

    data.distance_matrix.resize(costs.rows());
    for(int i=0; i<costs.rows(); ++i) {
        data.distance_matrix[i].resize(costs.cols());
        for(int j=0; j<costs.cols(); ++j) {
            data.distance_matrix[i][j] = static_cast<int64_t>(costs(i,j));
            if( i == j) {
                data.distance_matrix[i][j] = 0;
            }
        }
    }

    RoutingIndexManager manager(data.distance_matrix.size(),
                                data.num_vehicles,
                                data.depotNode);
    RoutingModel routing(manager);

    const int transit_callback_index = routing.RegisterTransitCallback(
    [&data, &manager](int64_t from_index, int64_t to_index) -> int64_t {
      // Convert from routing variable Index to distance matrix NodeIndex.
      auto from_node = manager.IndexToNode(from_index).value();
      auto to_node = manager.IndexToNode(to_index).value();
      return data.distance_matrix[from_node][to_node];
    });

    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);

    RoutingSearchParameters searchParameters = DefaultRoutingSearchParameters();
    searchParameters.set_first_solution_strategy(FirstSolutionStrategy::PATH_CHEAPEST_ARC);

    const Assignment* solution = routing.SolveWithParameters(searchParameters);

    // Inspect solution.
    std::vector<int> routeIndices;
    int64_t distance{0};
    int64_t index = routing.Start(0);
    while (routing.IsEnd(index) == false) {
      routeIndices.push_back(index);
      int64_t previous_index = index;
      index = solution->Value(routing.NextVar(index));
      distance += routing.GetArcCostForVehicle(previous_index, index, int64_t{0});
    }
    return std::make_pair(routeIndices, distance);
}

std::vector<int> routeplan(const Eigen::MatrixXd& costs) {
    int SITES = costs.rows();

    // NOTE(Jordan):
    // For some reason, or-tools doesn't like it when you use int64_t::max()
    // as a cost in the traveling salesman cost matrix.
    // Dividing by 1024 makes it work without overflowing.
    int64_t INF = std::numeric_limits<int64_t>::max() / 1024;

    std::vector<int> minRoute;
    int64_t minCost = INF;

    for(int endpoint=1; endpoint<SITES; ++endpoint) {
        Eigen::MatrixXd augCosts(SITES+1, SITES+1);
        augCosts.fill(INF);
        augCosts.block(0,0, costs.rows(), costs.cols()) = costs;
        augCosts(SITES, SITES) = 0.0;

        // START
        augCosts(SITES, 0) = 0.0;
        augCosts(0, SITES) = 0.0;

        // END
        augCosts(SITES, endpoint) = 0.0;
        augCosts(endpoint, SITES) = 0.0;

        auto [route, cost] = solveTSP(augCosts, SITES); 
        if( cost < minCost ) {
            minRoute = route;
            minCost = cost;
        }
        // fmt::print("COST {} ROUTE: ", cost);
        // for(const auto& x : route) { fmt::print("{} ", x); }
        // fmt::print("\n");
    }

    // Remove the imaginary depot city from the beginning of the route.
    minRoute.erase(minRoute.begin());

    // Because of symmetry, the route may traverse
    // from the endpoint to the landing site.
    // If the route is backward, flip it!
    if( minRoute[0] != 0 ) {
        std::reverse(minRoute.begin(), minRoute.end());
    }
    assert(minRoute[0] == 0);

    return minRoute;
}
std::vector<std::vector<Path>> planAllPairsSLOW(const std::vector<Path::State>& sites,
                                                const SlopeAtlas& slopeAtlas,
                                                const TerrainMapFloat& commsMap) {
    std::vector<std::vector<Path>> allPaths;
    allPaths.resize(sites.size());
    for(auto& pathList : allPaths) { pathList.resize(sites.size()); }

    int totalExpansions = 0;

    // Generate all combinations of two indices into the allSites vector.
    #pragma omp parallel for
    for(int a=0; a<sites.size()-1; ++a) {
        const Path::State start = sites[a];

        for(int b=a+1; b<sites.size(); ++b) {
            std::vector<Path::State> goals;
            goals.push_back(sites[b]);
            auto [paths, expansions] = multipathplan(slopeAtlas, commsMap, start, goals, a);
            #pragma omp critical
            {
                totalExpansions += expansions;
            }

            for( const auto& path : paths ) {
                if( path.states.size() == 0 ) { continue; }
                const auto& start = path.states[0];
                const auto& goal = path.states[path.states.size()-1];

                // We forgot which (a,b) pair this goal came from.
                // We have to find this goal in the goals list so we know where
                // to record things in the costs and dists matrices.
                {
                    const auto it = std::find(goals.begin(), goals.end(), goal);
                    assert( it != goals.cend() );
                    b += std::distance(goals.begin(), it);
                }
                allPaths[a][b] = path;
                allPaths[b][a] = reverse(path);
            }
        }
    }
    fmt::print("Total Expansions (Separate): {}\n", totalExpansions);

    return allPaths;
}

std::vector<std::vector<Path>> planAllPairs(const std::vector<Path::State>& sites,
                                            const SlopeAtlas& slopeAtlas,
                                            const TerrainMapFloat& commsMap) {
    std::vector<std::vector<Path>> allPaths;
    allPaths.resize(sites.size());
    for(auto& pathList : allPaths) { pathList.resize(sites.size()); }

    int totalExpansions = 0;

    // Generate all combinations of two indices into the allSites vector.
    #pragma omp parallel for
    for(int a=0; a<sites.size()-1; ++a) {
        const Path::State start = sites[a];

        std::vector<Path::State> goals;
        for(int b=a+1; b<sites.size(); ++b) {
            goals.push_back(sites[b]);
        }

        auto [paths, expansions] = multipathplan(slopeAtlas, commsMap, start, goals, a);
        #pragma omp critical
        {
            totalExpansions += expansions;
        }

        for( const auto& path : paths ) {
            if( path.states.size() == 0 ) { continue; }
            const auto& start = path.states[0];
            const auto& goal = path.states[path.states.size()-1];

            // We forgot which (a,b) pair this goal came from.
            // We have to find this goal in the goals list so we know where
            // to record things in the costs and dists matrices.
            int b = a+1;
            {
                const auto it = std::find(goals.begin(), goals.end(), goal);
                assert( it != goals.cend() );
                b += std::distance(goals.begin(), it);
            }
            allPaths[a][b] = path;
            allPaths[b][a] = reverse(path);
        }
    }
    fmt::print("Total Expansions (Fused): {}\n", totalExpansions);

    // Draw all paths on their own map.
    for(int a=0; a<allPaths.size(); ++a) {
        for(int b=0; b<allPaths[0].size(); ++b) {
            if( a >= b ) { continue; }
            const auto& path = allPaths[a][b];

            TerrainMapFloat pathMap = slopeAtlas.absolute;

            // Draw the start.
            {
                const auto& start = sites[a];
                double sx = pathMap.j2x(start.j);
                double sy = pathMap.i2y(start.i);
                drawCircle(pathMap, sx,sy, 100, 2.0);
                double dX = std::cos(M_PI/180.0*directionToDegrees(start.d));
                double dY = std::sin(M_PI/180.0*directionToDegrees(start.d));
                drawCircle(pathMap, sx+3.0*dX, sy+3.0*dY, 100, 1.0);
            }

            // Draw the goal.
            {
                const auto&  goal = sites[b];
                double gx = pathMap.j2x(goal.j);
                double gy = pathMap.i2y(goal.i);
                drawCircle(pathMap, gx, gy, 100, 2.0);
                double dX = std::cos(M_PI/180.0*directionToDegrees(goal.d));
                double dY = std::sin(M_PI/180.0*directionToDegrees(goal.d));
                drawCircle(pathMap, gx+3.0*dX, gy+3.0*dY, 100, 1.0);
            }

            // Draw the path.
            for(const auto& p : path.states) {
                pathMap(p.i, p.j) = 100;
            }

            // Save the map.
            saveEXR(pathMap, config.outputDir+fmt::format("paths_{:02}_{:02}.exr", a, b)); 
        }
    }

    // Draw all of the paths on a single map.
    TerrainMapFloat allPathsMap = slopeAtlas.absolute;
    for(int a=0; a<allPaths.size(); ++a) {
        for(int b=0; b<allPaths[0].size(); ++b) {
            if( a >= b ) { continue; }
            const auto& path = allPaths[a][b];

            // Draw the start.
            const auto& start = sites[a];
            double sx = allPathsMap.j2x(start.j);
            double sy = allPathsMap.i2y(start.i);
            drawCircle(allPathsMap, sx,sy, 100, 2.0);

            // Draw the goal.
            const auto& goal = sites[b];
            double gx = allPathsMap.j2x(goal.j);
            double gy = allPathsMap.i2y(goal.i);
            drawCircle(allPathsMap, gx, gy, 100, 2.0);

            // Draw the path.
            for(const auto& p : path.states) {
                allPathsMap(p.i, p.j) = 100;
            }
        }
    }
    // Save the map.
    saveEXR(allPathsMap, config.outputDir+"paths.exr"); 
    
    return allPaths;
}



// Sort vantages by their angle relative to (siteX, siteY).

int main(int argc, char* argv[]) {
    namespace cm = tinycolormap;
    
    // Parse the command line and populate the global config struct.
    parseCommandLine(config, argc, argv);

    // Create the output directory if it doesn't exist already.
    if( !std::filesystem::exists(config.outputDir) ) {
        std::filesystem::create_directory(config.outputDir);
    }

    // Read the terrain mesh.
    TerrainMesh tmesh(config.meshFile);

    // Construct elevation, priority, and slope maps.
    auto [elevationMap, priorityMap, slopeAtlas] = buildTerrainMaps(tmesh, config.mapPitch);

    // Compute landing site coordinates.
    double landingSiteX = config.landingSiteX;
    double landingSiteY = config.landingSiteY;
    double landingSiteZ = elevationMap.atXY(landingSiteX, landingSiteY);

    int landingSiteI = elevationMap.y2i(landingSiteY);
    int landingSiteJ = elevationMap.x2j(landingSiteX);

    if( landingSiteI < 0 || landingSiteI >= elevationMap.rows ||
        landingSiteJ < 0 || landingSiteJ >= elevationMap.cols ) { throw std::runtime_error(
            fmt::format("Landing site at ({}, {}) is outside of map boundaries ({}, {}).",
                        landingSiteX, landingSiteY, elevationMap.height(), elevationMap.width()));
    }

    Path::State landingState;
    landingState.i = landingSiteI;
    landingState.j = landingSiteJ;
    landingState.d = (landingSiteJ > elevationMap.rows/2) ? Direction::W : Direction::E;

    // Construct lander communications map.
    TerrainMapFloat commsMap = buildCommsMap(tmesh, elevationMap, landingState);

    // Flood-fill reachable safe terrain.
    TerrainMapU8 reachMap = buildReachabilityMap(commsMap, slopeAtlas, landingState);
    {
        for(int d=0; d<8; ++d) {
            TerrainMapFloat map(reachMap.width(), reachMap.height(), reachMap.pitch);
            for(int i=0; i<map.rows; ++i) {   
                for(int j=0; j<map.cols; ++j) {
                    Path::State s {.i=i, .j=j, .d=(Direction)d};
                    map(i,j) = checkMapBit(reachMap, s) ? 1.0f : 0.0f;
                }
            }
            drawCircle(map, landingSiteX, landingSiteY, 2.0f, 4.0);
            saveEXR(map, config.outputDir+fmt::format("reach_{}.exr", d));
        }
    }

    // Generate visibility probes.
    const auto [probes, probeMap] = generateVisibilityProbes(priorityMap, elevationMap);

    // Generate candidate vantages.
    const auto [candidates, candidateMap] = generateVantageCandidates(tmesh, reachMap, elevationMap, probes);

    // Select the best vantages from all of the candidates.
    auto vantages = selectVantages(candidates, probes);

    // Sort vantages in counterclockwise order.
    {
        // Sort vantages counterclockwise around their centroid.
        double centroidX = 0; double centroidY = 0;
        for(const auto& v : vantages) { centroidX += v.x; centroidY += v.y; }
        centroidX /= vantages.size(); centroidY /= vantages.size();
        vantages = sortCCW(vantages, centroidX, centroidY);
    }

    // Map combined coverage from all vantages.
    auto coverageMap = buildCoverageMap(tmesh, elevationMap, vantages);

    // Save landing site and vantages to xyz file.
    {
        std::ofstream file;
        file.open(config.outputDir+"sites.xyz");
        file << fmt::format("{} {} {}\n", landingSiteX, landingSiteY, landingSiteZ);
        for(const auto& v : vantages) {
            file << fmt::format("{} {} {}\n", v.x, v.y, v.z + config.roverHeight);
        }
        file.close();
    }

    // Save some maps.
    {
        {
            auto map = elevationMap;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"elevation.exr");
        }
        {
            auto map = slopeAtlas.absolute;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            drawTriangle(map, landingSiteX, landingSiteY, 0, 100, 5, 10);
            saveEXR(map, config.outputDir+"slope.exr");
        }
        {
            auto map = slopeAtlas.north;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeN.exr");
        }
        {
            auto map = slopeAtlas.east;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeE.exr");
        }
        {
            auto map = slopeAtlas.northeast;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeNE.exr");
        }
        {
            auto map = slopeAtlas.southeast;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeSE.exr");
        }


        auto negate = [](const TerrainMapFloat& m) {
            TerrainMapFloat newm = m;
            for(int i=0; i<newm.rows; ++i) {
            for(int j=0; j<newm.cols; ++j) {
                newm(i,j) = -1*m(i,j);
            }}
            return newm;
        };
        {
            auto map = negate(slopeAtlas.north);
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeS.exr");
        }
        {
            auto map = negate(slopeAtlas.east);
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeW.exr");
        }
        {
            auto map = negate(slopeAtlas.northeast);
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeSW.exr");
        }
        {
            auto map = negate(slopeAtlas.southeast);
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeNW.exr");
        }
        {
            auto map = priorityMap;
            drawCircle(map, landingSiteX, landingSiteY, 10, 3.0);
            saveEXR(map, config.outputDir+"priority.exr");
        }
        {
            auto map = commsMap;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"comms.exr");

            int upsample = 4;
            double r = upsample * 2.5/slopeAtlas.pitch();
            ImageRGB img(commsMap, cm::ColormapType::Plasma, upsample);
            drawCircle(img, landingSiteI*upsample, landingSiteJ*upsample, r, cm::Color(1.0, 0.0, 0.0));
            savePNG(img, config.outputDir+"comms.png");
        }
        {
            auto map = probeMap;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"probes.exr");
        }
        {
            auto map = candidateMap;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"candidates.exr");
        }
        {
            auto map = coverageMap;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"coverage.exr");
        }
    }

    // Draw all vantages on a single map.
    auto vantageMap = slopeAtlas.absolute;
    {
        double maxCoverage = 0;
        for(const auto& v : vantages) { maxCoverage = std::max(maxCoverage, v.totalCoverage); }
        for(const auto& v : vantages) {
            double markerValue = 120 + 10 * v.totalCoverage / maxCoverage;
            drawCircle(vantageMap, v.x, v.y, markerValue, 2.0); 
            double dX = std::cos(M_PI/180.0*directionToDegrees(v.dir));
            double dY = std::sin(M_PI/180.0*directionToDegrees(v.dir));
            drawCircle(vantageMap, v.x+3.0*dX, v.y+3.0*dY, markerValue, 1.0);
        }
    }
    drawCircle(vantageMap, landingSiteX, landingSiteY, 100, 3.0);
    saveEXR(vantageMap, config.outputDir+"vantages.exr");

    // Draw separate coverage maps for each vantage.
    std::vector<TerrainMap<float>> coverageMaps;
    for(int vi=0; vi<vantages.size(); ++vi) {
        const auto& v = vantages[vi];
        std::vector<Vantage> tmp; tmp.push_back(v);
        auto coverageMap = buildCoverageMap(tmesh, elevationMap, tmp);
        coverageMaps.push_back(coverageMap);
        int j = vantageMap.x2j(v.x);
        int i = vantageMap.y2i(v.y);
        drawCircle(coverageMap, v.x, v.y, vantages.size()+1, 2.0);
        drawCircle(coverageMap, landingSiteX, landingSiteY, vantages.size()+10, 3.0);

        double dX = std::cos(M_PI/180.0*directionToDegrees(v.dir));
        double dY = std::sin(M_PI/180.0*directionToDegrees(v.dir));
        drawCircle(coverageMap, v.x+3*dX, v.y+3*dY, vantages.size()+10, 1.0);

        saveEXR(coverageMap, config.outputDir+fmt::format("coverage_{:02}.exr",vi));
    }

    // Create a vector of sites relevant for planning.
    std::vector<Path::State> allSites;

    // The first site in the vector is the landing site.
    allSites.push_back(landingState);

    // The remaining sites are the sorted vantages.
    for(const auto& v : vantages) {
        Path::State s;
        s.i = slopeAtlas.y2i(v.y);
        s.j = slopeAtlas.x2j(v.x);
        s.d = v.dir;
        allSites.push_back(s);
    }

    // Save the sites to a csv file.
    {
        std::ofstream file;
        file.open(config.outputDir+"sites.csv");
        file << fmt::format("x,y,z,i,j,d,slope,n,ne,e,se\n");
        for(const auto& s : allSites) {
            Vantage v;
            v.x = elevationMap.j2x(s.j);
            v.y = elevationMap.j2x(s.i);
            v.z = elevationMap(s.i, s.j) + config.roverHeight;
            file << fmt::format("{:0.3f},{:0.3f},{:0.3f}, {},{},{}\n",
                                v.x, v.y, v.z,
                                s.i, s.j, directionToString(s.d));
        }
        file.close();
    }

    // Compute paths between all pairs of k vantages plus the landing site.
    auto tstart = std::chrono::high_resolution_clock::now();
    const auto paths2 = planAllPairsSLOW(allSites, slopeAtlas, commsMap);
    auto tmid = std::chrono::high_resolution_clock::now();
    const auto paths = planAllPairs(allSites, slopeAtlas, commsMap);
    auto tstop = std::chrono::high_resolution_clock::now();
    fmt::print("SLOW: {} sec\n", std::chrono::duration_cast<std::chrono::milliseconds>(tmid-tstart).count()/1000.0);
    fmt::print("FAST: {} sec\n", std::chrono::duration_cast<std::chrono::milliseconds>(tstop-tmid).count()/1000.0);

    // Organize path costs into convenient matrices for route planning.
    Eigen::MatrixXd costs(allSites.size(), allSites.size()); costs.fill(-1);
    Eigen::MatrixXd dists(allSites.size(), allSites.size()); dists.fill(-1);
    for(int a=0; a<paths.size(); ++a) {
        for(int b=0; b<paths[0].size(); ++b) {
            costs(a,b) = paths[a][b].cost;
            dists(a,b) = paths[a][b].dist;
        }
    }
    // fmt::print("Costs:\n{}\n", costs);
    // fmt::print("Dists:\n{}\n", dists);

    // Compute exploration route.
    auto route = routeplan(costs);

    if( route.size() < allSites.size() ) {
        fmt::print("Oh no! I failed to plan a route to all vantages.\n");
    }

    // Chain paths together to create the final path.
    Path path = assembleRoute(route, paths);
    fmt::print("Final Cost: {:.3f}\n", path.cost);
    fmt::print("Final Dist: {:.3f}\n", path.dist);

    // Draw the final route!
    TerrainMapFloat routeMap = vantageMap;
    drawCircle(routeMap, landingSiteX, landingSiteY, 100, 3.0);
    for(const auto& p : path.states) {
        routeMap(p.i, p.j) = 100;
    }
    saveEXR(routeMap, config.outputDir+"route.exr"); 

    // Save the route to an xyz file.
    {
        std::ofstream file;
        file.open(config.outputDir+"route.xyz");
        for(const auto& p : path.states) {
            Vantage v;
            v.x = elevationMap.j2x(p.j);
            v.y = elevationMap.i2y(p.i);
            v.z = elevationMap(p.i, p.j) + config.roverHeight;
            file << fmt::format("{} {} {}\n", v.x, v.y, v.z);
        }
        file.close();
    }

    // Save the route to a csv file.
    {
        std::ofstream file;
        file.open(config.outputDir+"route.csv");
        file << fmt::format("x,y,z,dir_num,dir_name,is_vantage,"
                            "rover_roll_deg,rover_pitch_deg,rover_yaw_deg\n");
        for(const auto& p : path.states) {
            Vantage v;
            v.x = elevationMap.j2x(p.j);
            v.y = elevationMap.j2x(p.i);
            v.z = elevationMap(p.i, p.j) + config.roverHeight;
            file << fmt::format("{:.2f},{:.2f},{:.2f},{},{},", v.x, v.y, v.z, (int)p.d, directionToString(p.d));

            // Is this state a vantage?
            int is_vantage = 0;
            if( p != path.states[0] && std::find(allSites.begin(), allSites.end(), p) != allSites.end() ) {
                is_vantage = 1;
            }
            file << fmt::format("{},", is_vantage);

            double rover_roll  = computeLateralSlope(p, slopeAtlas);
            double rover_pitch = computeLongitudinalSlope(p, slopeAtlas);
            double rover_yaw   = 45 * (int)p.d;
            file << fmt::format("{:.2f},{:.2f},{:.2f}\n", rover_roll, rover_pitch, rover_yaw);
        }
        file.close();
    }

    // Animate the final path.
    if( config.drawAnimation ) {
        // Create a directory for animation frames.
        std::string animDir = config.outputDir + "/anim";
        if( !std::filesystem::exists(animDir) ) {
            std::filesystem::create_directory(animDir);
        }

        int upsample = 4;
        ImageRGB baseImg(slopeAtlas.absolute, cm::ColormapType::Plasma, upsample);

        #pragma omp parallel for
        for(int f=0; f<path.states.size(); ++f) {
            //fmt::print("Animating frame {}...\n", f);
            const auto& s = path.states[f];
            ImageRGB img = baseImg;

            // Draw the landing site.
            double r = upsample * 5.0/slopeAtlas.pitch();
            drawCircle(img, landingSiteI*upsample, landingSiteJ*upsample, r, r+4, cm::Color(1.0, 1.0, 1.0));

            // Draw all coverage the rover has already visited.
            // Don't check the first state because it will match the landing site.
            for(int ff=1; ff<=f; ++ff) {
                const auto s = path.states[ff];
                const auto vIt = std::find(allSites.begin(), allSites.end(), s);
                if( vIt != allSites.end() ) {
                    int vIdx = std::distance(allSites.begin(), vIt);

                    auto covMap = coverageMaps[vIdx-1];
                    for(int i=0; i<covMap.rows; ++i) {
                        for(int j=0; j<covMap.cols; ++j) {
                            if( covMap(i,j) ) {
                                for(int ii=0; ii<upsample; ++ii) {
                                    for(int jj=0; jj<upsample; ++jj) {
                                        const double alpha = 0.8;
                                        img(i*upsample+ii, j*upsample+jj, 0) = img(i*upsample+ii, j*upsample+jj, 0) * alpha +   0*(1-alpha);
                                        img(i*upsample+ii, j*upsample+jj, 1) = img(i*upsample+ii, j*upsample+jj, 1) * alpha + 255*(1-alpha);
                                        img(i*upsample+ii, j*upsample+jj, 2) = img(i*upsample+ii, j*upsample+jj, 2) * alpha +   0*(1-alpha);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Draw a trail behind the rover.
            for(int ff=0; ff<f; ++ff) {
                const auto& s0 = path.states[ff];
                const auto& s1 = path.states[ff+1];
                drawLine(img, s0.i*upsample,s0.j*upsample, s1.i*upsample,s1.j*upsample, cm::Color(1.0, 1.0, 1.0), 1.0*upsample);
            }

            // Draw all vantages the rover has already visited.
            for(int ff=1; ff<f; ++ff) {
                const auto s = path.states[ff];
                if( std::find(allSites.begin(), allSites.end(), s) != allSites.end() ) {
                    drawTriangle(img, s.i*upsample, s.j*upsample, (int)s.d*45, cm::Color(1.0, 1.0, 1.0), 5*upsample, 8*upsample);
                }
            }

            // Draw the rover.
            drawTriangle(img, s.i*upsample, s.j*upsample, (int)s.d*45, cm::Color(1.0, 0.0, 0.0), 5*upsample, 8*upsample);
            savePNG(img, animDir+"/"+fmt::format("anim_{:05}.png", f));
        }
    }

    return 0;
}
