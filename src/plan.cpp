#include "plan.h"
#include <vector>
#include "terrainmesh.h"
#include "terrainmap.h"
#include "path.h"
#include <Eigen/Dense>
#include "ortools/constraint_solver/routing.h"
#include "ortools/constraint_solver/routing_enums.pb.h"
#include "ortools/constraint_solver/routing_index_manager.h"
#include "ortools/constraint_solver/routing_parameters.h"
#include "flat_hash_map.hpp"

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
