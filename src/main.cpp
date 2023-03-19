#include "image.h"
#include "terrainmesh.h"
#include "terrainmap.h"
#include "argh.h"
#include "flat_hash_map.hpp"
#include "toml.hpp"
#include "priority_queue.h"
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

struct PlannerConfiguration {
    // Basic
    int numVantages       =   18; // 
    double mapPitch       =  1.0; // meters
    std::string outputDir = "./"; //

    // Lander
    double landingSiteX = 100; // meters
    double landingSiteY = 100; // meters
    double landerHeight = 1.0; // meters

    // Rover
    double roverHeight                 =  1.0; // meters
    double roverSpeed                  = 0.01; // m/s
    double roverFOV                    =  180; // degrees
    double roverLateralSlopeLimit      = 12.0; // degrees
    double roverLongitudinalSlopeLimit = 20.0; // degrees
    double roverPointTurnSlopeLimit    =  5.0; // degrees

    // Advanced
    std::string meshFile = "../meshes/lmp-interior.ply"; //

    int    numProbes     = 10000; // 
    int    numCandidates = 10000; // 

    double maxVisAngle          =   55; // degrees
    double maxVisRange          =  300; // meters
    double minVantageSeparation = 20.0; // meters

    std::array<double, 10> visMultipliers = {1.0, 1.0, 0.8, 0.5, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1}; //
 
    double distCostMultiplier  =  1.0; //
    double slopeCostMultiplier = 16.0; //
    double turnCostMultiplier  =  1.0; //
    double heuristicMultiplier =  1.0; //
};

// Global configuration structure
PlannerConfiguration config;

void parseConfigFile(const std::string& fileName) {
    const std::string fileDir = std::filesystem::path(fileName).parent_path();

    const toml::table tbl = toml::parse_file(fileName);
    config.numVantages  = tbl["Basic"]["NumVantages"].value_or(config.numVantages);
    config.mapPitch     = tbl["Basic"]["MapPitch"].value_or(config.mapPitch);
    config.outputDir    = tbl["Basic"]["OutputDir"].value_or(config.outputDir);

    config.landingSiteX = tbl["Lander"]["LocationXY"][0].value_or(config.landingSiteX);
    config.landingSiteY = tbl["Lander"]["LocationXY"][1].value_or(config.landingSiteY);
    config.landerHeight = tbl["Lander"]["Height"].value_or(config.landerHeight);

    config.roverHeight                 = tbl["Rover"]["Height"].value_or(config.roverHeight);
    config.roverSpeed                  = tbl["Rover"]["Speed"].value_or(config.roverSpeed);
    config.roverFOV                    = tbl["Rover"]["FOV"].value_or(config.roverFOV);
    config.roverLateralSlopeLimit      = tbl["Rover"]["LateralSlopeLimit"].value_or(config.roverLateralSlopeLimit);
    config.roverLongitudinalSlopeLimit = tbl["Rover"]["LongitudinalSlopeLimit"].value_or(config.roverLongitudinalSlopeLimit);
    config.roverPointTurnSlopeLimit    = tbl["Rover"]["PointTurnSlopeLimit"].value_or(config.roverPointTurnSlopeLimit);

    config.meshFile             = fileDir + "/" + tbl["Advanced"]["MeshFile"].value_or(config.meshFile);
    config.numProbes            = tbl["Advanced"]["NumProbes"].value_or(config.numProbes);
    config.numCandidates        = tbl["Advanced"]["NumCandidates"].value_or(config.numCandidates);
    config.maxVisAngle          = tbl["Advanced"]["MaxVisAngle"].value_or(config.maxVisAngle);
    config.maxVisRange          = tbl["Advanced"]["MaxVisRange"].value_or(config.maxVisRange);
    config.minVantageSeparation = tbl["Advanced"]["MinVantageSeparation"].value_or(config.minVantageSeparation);

    config.distCostMultiplier  = tbl["Advanced"]["DistCostMultiplier"].value_or(config.distCostMultiplier);
    config.slopeCostMultiplier = tbl["Advanced"]["SlopeCostMultiplier"].value_or(config.slopeCostMultiplier);
    config.turnCostMultiplier  = tbl["Advanced"]["TurnCostMultiplier"].value_or(config.turnCostMultiplier);
    config.heuristicMultiplier = tbl["Advanced"]["HeuristicMultiplier"].value_or(config.heuristicMultiplier);
}

void parseCommandLine(int argc, char* argv[]) {
    argh::parser cmdl(argc, argv);

    if( cmdl[{ "-h", "--help", "--Help" }] ) {
        fmt::print("usage: planranger [options] [config.toml]\n");
        fmt::print("  options:\n");
        fmt::print("    [Basic]\n");
        fmt::print("      -h,   --Help                    Show this help message.\n");
        fmt::print("      -nv,  --NumVantages             The number of vantages to include in the route.\n");
        fmt::print("      -p,   --MapPitch                The resolution (meters) to use for generating maps.\n");
        fmt::print("      -o,   --OutputDir               A directory to store the planner output.\n");
        fmt::print("\n");

        fmt::print("    [Lander]\n");
        fmt::print("      -x,   --LandingSiteX            The x-coordinate of the landing site.\n");
        fmt::print("      -y,   --LandingSiteY            The y-coordinate of the landing site.\n");
        fmt::print("      -lh,  --LanderHeight            The height of the comms antenna on the lander.\n");
        fmt::print("\n");

        fmt::print("    [Rover]\n");
        fmt::print("      -rh,   --RoverHeight            The height of the comms antenna on the rover.\n");
        fmt::print("      -rs,   --RoverSpeed             The rover's driving speed in m/s .\n");
        fmt::print("      -fov,  --RoverFOV               The rover camera's field-of-view in degrees.\n");
        fmt::print("\n");

        fmt::print("      -lat,  --LateralSlopeLimit      The rover's max lateral slope in degrees.\n");
        fmt::print("      -lon,  --LongitudinalSlopeLimit The rover's max longitudinal slope in degrees.\n");
        fmt::print("      -turn, --PointTurnSlopeLimit    The steepest slope (in degrees) on which the rover can turn.\n");
        fmt::print("\n");

        fmt::print("    [Advanced]\n");
        fmt::print("      -mesh,   --MeshFile             A *.ply file representing the terrain.\n");
        fmt::print("      -np,     --NumProbes            The number of visibility probes to use for vantage selection.\n");
        fmt::print("      -nc,     --NumCandidates        The number of vantage candidates to consider for vantage selection.\n");
        fmt::print("      -va,     --MaxVisAngle          The maximum incidence angle (in degrees) that is considered a valid view hit.\n");
        fmt::print("      -vr,     --MaxVisRange          The maximum distance (in meters) beyond which a view hit is not counted.\n");
        fmt::print("      -minsep, --MinVantageSeparation The minimum distance between any two vantages.\n");
        fmt::print("\n");

        fmt::print("      -dmult,  --DistCostMultiplier   The unitless multiplier used to weight the distance cost term in the planning cost function.\n");
        fmt::print("      -smult,  --SlopeCostMultiplier  The unitless multiplier used to weight the slope cost term in the planning cost function.\n");
        fmt::print("      -tmult,  --TurnCostMultiplier   The unitless multiplier used to weight the turning cost term in the planning cost function.\n");
        fmt::print("      -hmult,  --HeuristicMultiplier  The unitless multiplier used to inflate the A* heuristic for faster, less-optimal planning.\n");
        fmt::print("\n");

        exit(0);
    }

    // If the first positional argument is a toml file, use it to configure the planner.
    if( !cmdl[1].empty() && cmdl[1].find(".toml") != std::string::npos ) {
        parseConfigFile(cmdl[1]);
    }

    // Basic
    cmdl({"NumVantages",      "nv"}, config.numVantages)                 >> config.numVantages;
    cmdl({"MapPitch",          "p"}, config.mapPitch)                    >> config.mapPitch;
    cmdl({"OutputDir",         "o"}, config.outputDir)                   >> config.outputDir;

    // Lander
    cmdl({"LandingSiteX",      "x"}, config.landingSiteX)                >> config.landingSiteX;
    cmdl({"LandingSiteY",      "y"}, config.landingSiteY)                >> config.landingSiteY;
    cmdl({"LanderHeight",     "lh"}, config.landerHeight)                >> config.landerHeight;

    // Rover
    cmdl({"RoverHeight",      "rh"}, config.roverHeight)                 >> config.roverHeight;
    cmdl({"RoverSpeed",       "rs"}, config.roverSpeed)                  >> config.roverSpeed;
    cmdl({"RoverFOV",        "fov"}, config.roverFOV)                    >> config.roverFOV;

    cmdl({"LateralSlopeLimit",       "lat"}, config.roverLateralSlopeLimit)      >> config.roverLateralSlopeLimit;
    cmdl({"LongitudinalSlopeLimit",  "lon"}, config.roverLongitudinalSlopeLimit) >> config.roverLongitudinalSlopeLimit;
    cmdl({"PointTurnSlopeLimit",    "turn"}, config.roverPointTurnSlopeLimit)    >> config.roverPointTurnSlopeLimit;

    // Advanced
    cmdl({"MeshFile",               "mesh"}, config.meshFile)             >> config.meshFile;
    cmdl({"NumProbes",                "np"}, config.numProbes)            >> config.numProbes;
    cmdl({"NumCandidates",            "nc"}, config.numCandidates)        >> config.numCandidates;
    cmdl({"MaxVisAngle",              "va"}, config.maxVisAngle)          >> config.maxVisAngle;
    cmdl({"MaxVisRange",              "vr"}, config.maxVisRange)          >> config.maxVisRange;
    cmdl({"MinVantageSeparation", "minsep"}, config.minVantageSeparation) >> config.minVantageSeparation;

    cmdl({"DistCostMultiplier",       "dmult"}, config.distCostMultiplier)  >> config.distCostMultiplier;
    cmdl({"SlopeCostMultiplier",      "smult"}, config.slopeCostMultiplier) >> config.slopeCostMultiplier;
    cmdl({"TurnCostMultiplier",       "tmult"}, config.turnCostMultiplier)  >> config.turnCostMultiplier;
    cmdl({"HeuristicMultiplier",      "hmult"}, config.heuristicMultiplier) >> config.heuristicMultiplier;
}

enum class Direction : int8_t { N = 0, NE, E, SE, S, SW, W, NW };

double directionToDegrees(const Direction& d) {
    const double dir2angle[] = {90, 45, 0, 315, 270, 225, 180, 135};
    return dir2angle[(int)d];
}

std::string directionToString(const Direction& d) {
    const std::string dir2name[] = {"N", "NE", "E", "SE", "S", "SW", "W", "NW"};
    return dir2name[(int)d];
}


struct SlopeAtlas {
    TerrainMapFloat absolute;
    TerrainMapFloat north;
    TerrainMapFloat northeast;
    TerrainMapFloat east;
    TerrainMapFloat southeast;

    SlopeAtlas(double mapWidth, double mapHeight, double mapPitch) :
        absolute(mapWidth, mapHeight, mapPitch), 
        north(mapWidth, mapHeight, mapPitch), 
        northeast(mapWidth, mapHeight, mapPitch), 
        east(mapWidth, mapHeight, mapPitch), 
        southeast(mapWidth, mapHeight, mapPitch) {}

    double j2x(size_t j) const { return absolute.j2x(j); }
    double i2y(size_t i) const { return absolute.i2y(i); }

    size_t x2j(double x) const { return absolute.x2j(x); }
    size_t y2i(double y) const { return absolute.y2i(y); }

    double width() const { return absolute.width(); }
    double height() const { return absolute.height(); }

    size_t rows() const { return absolute.rows; }
    size_t cols() const { return absolute.cols; }
    double pitch() const { return absolute.pitch; }
};

struct Path {
    struct State {
        int i = 0;
        int j = 0;
        Direction d = Direction::N;
        bool operator==(const State& rhs) const {
            return rhs.i==i && rhs.j==j && d==rhs.d;
        }
        bool sameLocation(const State& rhs) { return rhs.i==i && rhs.j==j; }
        bool sameDirection(const State& rhs) { return d==rhs.d; }
    };

    std::vector<State> states;

    double cost = 0.0;
    double dist = 0.0;
};

namespace std {
template <>
struct hash<Path::State> {
    std::size_t operator()(const Path::State& p) const {
        // pack two int32 into one uint64
        std::size_t res;
        res  = (std::size_t)p.i << 32;
        res |= (std::size_t)p.j & 0x00000000FFFFFFFF;
        // use boost hash_combine to mix in the direction as well.
        res ^= (std::size_t)p.d + 0x9e3779b9 + (res << 6) + (res >> 2);
        return res;
    }
};
}

double octileDistance(const Path::State& a, const Path::State& b) {
    // Branchless octile distance
    int di = std::abs(a.i-b.i);
    int dj = std::abs(a.j-b.j);
    constexpr double diagonal = std::sqrt(2);
    constexpr double twoCardinalMinusDiagonal = 2-diagonal;
    return (twoCardinalMinusDiagonal*abs(di-dj) + diagonal*(di+dj)) / 2;
}

Direction directionFromTo(const Path::State& from, const Path::State& to) {
    constexpr double tan_pi_8 = std::tan(M_PI/8.0);

    const int di = to.i-from.i;
    const int dj = to.j-from.j;

    if( di <= 0 and std::abs(dj) <= tan_pi_8 * -di ) { return Direction::N; }
    if( di >= 0 and std::abs(dj) <= tan_pi_8 *  di ) { return Direction::S; }

    if( dj >= 0 and std::abs(di) <= tan_pi_8 *  dj ) { return Direction::E; }
    if( dj <= 0 and std::abs(di) <= tan_pi_8 * -dj ) { return Direction::W; }

    if( di < 0 and dj > 0 ) { return Direction::NE; }
    if( di > 0 and dj > 0 ) { return Direction::SE; }

    if( di < 0 and dj < 0 ) { return Direction::NW; }
    if( di > 0 and dj < 0 ) { return Direction::SW; }

    return to.d;
};

Path append(const Path& a, const Path& b) {
    assert(a.cost >= 0); assert(a.dist >= 0);
    assert(b.cost >= 0); assert(b.dist >= 0);
    Path c;
    c.cost = a.cost + b.cost;
    c.dist = a.dist + b.dist;
    c.states = a.states;
    c.states.insert(c.states.end(), b.states.begin(), b.states.end());
    return c;
}

Path reverse(const Path& path) {
    Path rPath = path;
    std::reverse(rPath.states.begin(), rPath.states.end());
    return rPath;
}

uint8_t setBit(uint8_t byte, uint8_t bit) {
    return byte | ((uint8_t) 1) << bit;
};
uint8_t clearBit(uint8_t byte, uint8_t bit) {
    return byte & ~(((uint8_t) 1) << bit);
};
bool checkBit(uint8_t byte, uint8_t bit) {
    return (byte >> bit) & 0x01;
};
bool checkMapBit(const TerrainMapU8& map, const Path::State& state) {
    return checkBit(map(state.i, state.j), (int)state.d);
};
void setMapBit(TerrainMapU8& map, const Path::State& state) {
    map(state.i, state.j) = setBit(map(state.i, state.j), (int)state.d);
};

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

double computeLateralSlope(const Path::State& s, const SlopeAtlas& slopeAtlas) {
    double latSlope = 90.0;
    switch (s.d) {
        case Direction::N :
            latSlope = slopeAtlas.east(s.i, s.j);
            break;
        case Direction::NE :
            latSlope = slopeAtlas.southeast(s.i, s.j);
            break;
        case Direction::E :
            latSlope = -slopeAtlas.north(s.i, s.j);
            break;
        case Direction::SE :
            latSlope = -slopeAtlas.northeast(s.i, s.j);
            break;
        case Direction::S :
            latSlope = -slopeAtlas.east(s.i, s.j);
            break;
        case Direction::SW :
            latSlope = -slopeAtlas.southeast(s.i, s.j);
            break;
        case Direction::W :
            latSlope =  slopeAtlas.north(s.i, s.j);
            break;
        case Direction::NW :
            latSlope =  slopeAtlas.northeast(s.i, s.j);
            break;
        default:
            break;
    }
    return latSlope;
}

double computeLongitudinalSlope(const Path::State& s, const SlopeAtlas& slopeAtlas) {
    double lonSlope = 90.0;
    switch (s.d) {
        case Direction::N :
            lonSlope = slopeAtlas.north(s.i, s.j);
            break;
        case Direction::NE :
            lonSlope = slopeAtlas.northeast(s.i, s.j);
            break;
        case Direction::E :
            lonSlope =  slopeAtlas.east(s.i, s.j);
            break;
        case Direction::SE :
            lonSlope =  slopeAtlas.southeast(s.i, s.j);
            break;
        case Direction::S :
            lonSlope = -slopeAtlas.north(s.i, s.j);
            break;
        case Direction::SW :
            lonSlope = -slopeAtlas.northeast(s.i, s.j);
            break;
        case Direction::W :
            lonSlope = -slopeAtlas.east(s.i, s.j);
            break;
        case Direction::NW :
            lonSlope = -slopeAtlas.southeast(s.i, s.j);
            break;
        default:
            break;
    }
    return lonSlope;
}

// Generate the successors of a node.
std::vector<Path::State> getSuccessors(const Path::State& p, const SlopeAtlas& slopeAtlas, const TerrainMapFloat& commsMap) { 

    std::vector<Path::State> succs; succs.reserve(7+2);

    // You can only point turn if the terrain is flat enough.
    if( slopeAtlas.absolute(p.i, p.j) <= config.roverPointTurnSlopeLimit ) {
        for(int dd=-1; dd<=1; ++dd) {
            Path::State s;
            s.i=p.i;
            s.j=p.j;
            s.d=static_cast<Direction>((static_cast<int>(p.d)+dd+8)%8);
            if( s.d != p.d ) {
                succs.push_back(s);
            }
        }
    }

    double lonSlope = std::abs(computeLongitudinalSlope(p, slopeAtlas));
    double latSlope = std::abs(computeLateralSlope(p, slopeAtlas));

    // If the lateral and longitudinal slopes are safe,
    // try moving one step forward or one step backward.
    if( lonSlope <= config.roverLongitudinalSlopeLimit &&
        latSlope <= config.roverLateralSlopeLimit ) {

        Path::State sf = p;
        Path::State sb = p;

        constexpr int dI[8] = {-1, -1,  0,  1,  1,  1,  0, -1};
        constexpr int dJ[8] = { 0,  1,  1,  1,  0, -1, -1, -1};

        sf.i += dI[(int)p.d];
        sf.j += dJ[(int)p.d];

        sb.i -= dI[(int)p.d];
        sb.j -= dJ[(int)p.d];

        double fwdLonSlope = std::abs(computeLongitudinalSlope(sf, slopeAtlas));
        double fwdLatSlope = std::abs(computeLateralSlope(sf, slopeAtlas));
        double bwdLonSlope = std::abs(computeLongitudinalSlope(sb, slopeAtlas));
        double bwdLatSlope = std::abs(computeLateralSlope(sb, slopeAtlas));

        auto inBounds = [&commsMap](int i, int j) -> bool {
            return (0 <= i && i<commsMap.rows && 0 <= j && j<commsMap.cols);
        };

        if( inBounds(sf.i, sf.j) && commsMap(sf.i, sf.j) &&
            fwdLonSlope <= config.roverLongitudinalSlopeLimit &&
            fwdLatSlope <= config.roverLateralSlopeLimit ) {
            succs.push_back(sf);
        }
        if( inBounds(sb.i, sb.j) && commsMap(sb.i, sb.j) &&
            bwdLonSlope <= config.roverLongitudinalSlopeLimit &&
            bwdLatSlope <= config.roverLateralSlopeLimit ) {
            succs.push_back(sb);
        }
    }
    return succs; 
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

struct Probe {
    double x,y,z;
    double priority;
};

std::pair<std::vector<Probe>, TerrainMapFloat>
generateVisibilityProbes(const TerrainMapFloat& priorityMap, const TerrainMapFloat& elevationMap) {
    std::vector<Probe> probes;
    TerrainMapFloat probeMap(priorityMap.width(), priorityMap.height(), priorityMap.pitch);
    while(probes.size() < config.numProbes) {
        int i = priorityMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        int j = priorityMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        if( priorityMap(i,j) > 0 ) {
            Probe p;
            p.x = priorityMap.j2x(j);
            p.y = priorityMap.i2y(i);
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
    Direction dir = Direction::N;
    std::vector<bool> coverage;
    double totalCoverage = 0;
};

std::pair<std::vector<Vantage>, TerrainMapFloat>
generateVantageCandidates(const TerrainMesh& tmesh,
                          const TerrainMapU8& reachMap,
                          const TerrainMapFloat& elevationMap,
                          const std::vector<Probe> probes) {
    // Generate thousands of random, reachable points at rover height above the terrain.
    // NOTE(Jordan): I'm currently duplicating each candidate vantage eight times,
    // once per rover facing direction. This is fairly wasteful. We should probably evaluate
    // all directions at each location and pick the one direction that has the best coverage.
    std::vector<Vantage> candidates;
    while(candidates.size() < config.numCandidates) {
        Path::State state;
        state.i = reachMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        state.j = reachMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

        for(int d=0; d<8; ++d) {
            state.d = static_cast<Direction>(d);

            if( checkMapBit(reachMap, state) ) {
                Vantage v;
                v.x = reachMap.j2x(state.j);
                v.y = reachMap.i2y(state.i);
                v.z = elevationMap(state.i,state.j) + config.roverHeight;
                v.dir = state.d;
                candidates.push_back(v);
                fmt::print("{},{},{} -> {:<#010b}\n", state.i, state.j, directionToString(state.d), reachMap(state.i, state.j));
            }
        }
    }

    TerrainMapFloat candidateMap(reachMap.width(), reachMap.height(), reachMap.pitch);
    for(int i=0; i<candidateMap.rows; ++i) {
        for(int j=0; j<candidateMap.cols; ++j) {
            candidateMap(i,j) = (reachMap(i,j) == 0x00) ? 0.0f : 1.0f;
        }
    }

    int progress = 0;

    // For every candidate, trace a ray to all view coverage probes.
    #pragma omp parallel for shared(progress)
    for(int ci=0; ci<candidates.size(); ++ci) {
        #pragma omp critical
        {
            fmt::print("[{}/{}] Constructing visibility map.\n", progress++, candidates.size());
        }
        auto& candidate = candidates[ci];
        candidate.coverage.resize(probes.size());

        for(int pi = 0; pi<probes.size(); ++pi) {
            auto& p = probes[pi];

            // Is this ray in front of the rover and within its field of view?
            // FIXME(Jordan): Do this calculation in 3d, not in 2d.
            {
                const double roverAngle = directionToDegrees(candidate.dir) * M_PI/180.0;
                const double rx = std::cos(roverAngle);
                const double ry = std::sin(roverAngle);

                const double cpx = p.x-candidate.x;
                const double cpy = p.y-candidate.y;
                const double cpnorm = std::sqrt(cpx*cpx+cpy*cpy);

                const double viewAngle = 180.0/M_PI * std::acos( (rx*cpx+ry*cpy) / cpnorm );
                if( std::abs(viewAngle) > config.roverFOV/2 ) {
                    continue;
                }
            }

            TerrainMesh::Ray ray;
            ray.oX = candidate.x;
            ray.oY = candidate.y;
            ray.oZ = candidate.z + config.roverHeight;
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
                if( hitAngle < config.maxVisAngle && hitDist < config.maxVisRange && std::abs(rayNorm-hitDist) < 0.05*rayNorm ) {
                    candidate.coverage[pi] = true;
                    candidate.totalCoverage += p.priority;
                }
            }
        }
        int i = candidateMap.y2i(candidate.y);
        int j = candidateMap.x2j(candidate.x);
        candidateMap(i,j) = 10+std::max<float>(candidate.totalCoverage, candidateMap(i,j));
    }
    return std::make_pair(candidates, candidateMap);
}

std::vector<Vantage> selectVantages(const std::vector<Vantage>& candidates,
                                    const std::vector<Probe>& probes) {

    std::vector<Vantage> vantages;
    std::unordered_set<int> taken;
    std::vector<unsigned char> visCounters(probes.size(), 0);

    // Make k selections. Choose the candidate that produces the greatest *new* coverage.
    for(int k = 0; k<config.numVantages; ++k) {
        fmt::print("[{}/{}] Selecting vantages.\n", k, config.numVantages);

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
                double d2 = (t.x-c.x)*(t.x-c.x)+(t.y-c.y)*(t.y-c.y)+(t.z-c.z)*(t.z-c.z);
                if( d2 < config.minVantageSeparation*config.minVantageSeparation ) {
                    tooClose = true;
                    break;
                }
            }
            if( tooClose ) { continue; }

            // Otherwise, calculate "how good" this candidate is
            // based on the visibility probes it can see.
            // There are three factors to consider here:
            //   1. Does this candidate "see" new probes that haven't been seen yet? (countWeight)
            //   2. Are the probes this candidate can see in the center of the rover's FOV? (angleWeight)
            //   3. The science priority assigned to each probe by the user. (probes[pi].priority)
            for(int pi=0; pi<probes.size(); ++pi) {
                // If this probe is visible, add to the score for this candidate.
                if( c.coverage[pi] ) {
                    // How many times has this probe been viewed?
                    int visCount = std::clamp<int>(visCounters[pi], 0, 10);
                    double countWeight = config.visMultipliers[visCount];

                    // Compute the angle between (1) the camera's normal vector
                    // and (2) the viewing vector from the candidate to the probe.
                    double visAngle = 0.0;
                    {
                        // (1)
                        Eigen::Vector3f cam;
                        double camAngle = M_PI/180.0 * 45*(int)c.dir;
                        cam << std::sin(camAngle), std::cos(camAngle), 0.0;

                        // (2)
                        Eigen::Vector3f vis;
                        vis << probes[pi].x-c.x, probes[pi].y-c.y, probes[pi].z-c.z;
                        vis /= vis.norm();

                        // How far (in deg) is the view vector from the center of the rover's FOV?
                        visAngle = 180.0/M_PI * std::acos(vis.dot(cam));
                    }
                    // Weight views in the center of the FOV more than views on the edges of the FOV.
                    double angleWeight = std::exp(-(visAngle/90.0)*(visAngle/90.0));

                    // Add to the score for this vantage candidate.
                    scores[ci] += countWeight * angleWeight * probes[pi].priority;
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

        // If the best score is zero, give up because there are no more candidates worth considering.
        if( scores[best_ci] <= 0 ) {
            fmt::print("Warning: Only {} useful vantages were found of the {} that were requested.\n",
                vantages.size(), config.numVantages);
            break;
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

TerrainMapFloat buildCoverageMap(const TerrainMesh& tmesh,
                                 const TerrainMapFloat& elevationMap,
                                 const std::vector<Vantage>& vantages) {

    TerrainMapFloat coverageMap(elevationMap.width(), elevationMap.height(), elevationMap.pitch);

    for(int vi=0; vi<vantages.size(); ++vi) {
        auto& v = vantages[vi];
        fmt::print("[{}/{}] Generating coverage maps.\n", vi, vantages.size());
        #pragma omp parallel for
        for(int ri=0; ri<coverageMap.rows; ++ri) {
            for(int rj=0; rj<coverageMap.cols; ++rj) {
                const double terrainX = elevationMap.j2x(rj);
                const double terrainY = elevationMap.i2y(ri);
                const double terrainZ = elevationMap(ri, rj);

                // Is this ray in front of the rover and within its field of view?
                // FIXME(Jordan): Do this calculation in 3d, not in 2d.
                {
                    const double roverAngle = directionToDegrees(v.dir) * M_PI/180.0;
                    const double rx = std::cos(roverAngle);
                    const double ry = std::sin(roverAngle);

                    const double vtx = terrainX-v.x;
                    const double vty = terrainY-v.y;
                    const double vtnorm = std::sqrt(vtx*vtx+vty*vty);

                    const double viewAngle = 180.0/M_PI * std::acos( (rx*vtx+ry*vty) / vtnorm );
                    if( std::abs(viewAngle) > config.roverFOV/2 ) {
                        continue;
                    }
                }

                TerrainMesh::Ray ray;
                ray.oX = v.x;
                ray.oY = v.y;
                ray.oZ = v.z+config.roverHeight;
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
                    if( hitAngle < config.maxVisAngle && hitDist < config.maxVisRange && std::abs(rayNorm-hitDist) < rayNorm*0.05 ) { 
                        coverageMap(ri, rj)++;
                    }
                } 
            }
        }
    }
    for(const auto& v : vantages) {
        const double roverAngle = directionToDegrees(v.dir) * M_PI/180.0;
        const double rx = std::cos(roverAngle);
        const double ry = std::sin(roverAngle);
        drawCircle(coverageMap, v.x+4*rx, v.y+4*ry, (float)2*vantages.size()+10, 1.0);
        drawCircle(coverageMap, v.x, v.y, (float)2*vantages.size(), 2.0);
    }
    return coverageMap;
}

std::vector<Path> multipathplan(const SlopeAtlas& slopeAtlas,
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

    return paths;
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

    // Generate all combinations of two indices into the allSites vector.
    #pragma omp parallel for
    for(int a=0; a<sites.size()-1; ++a) {
        const Path::State start = sites[a];

        for(int b=a+1; b<sites.size(); ++b) {
            std::vector<Path::State> goals;
            goals.push_back(sites[b]);
            auto paths = multipathplan(slopeAtlas, commsMap, start, goals, a);

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
    return allPaths;
}

std::vector<std::vector<Path>> planAllPairs(const std::vector<Path::State>& sites,
                                            const SlopeAtlas& slopeAtlas,
                                            const TerrainMapFloat& commsMap) {
    std::vector<std::vector<Path>> allPaths;
    allPaths.resize(sites.size());
    for(auto& pathList : allPaths) { pathList.resize(sites.size()); }

    // Generate all combinations of two indices into the allSites vector.
    #pragma omp parallel for
    for(int a=0; a<sites.size()-1; ++a) {
        const Path::State start = sites[a];

        std::vector<Path::State> goals;
        for(int b=a+1; b<sites.size(); ++b) {
            goals.push_back(sites[b]);
        }

        auto paths = multipathplan(slopeAtlas, commsMap, start, goals, a);

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

Path assembleRoute(const std::vector<int>& route,
                   const std::vector<std::vector<Path>> paths) {
    // Walk the route, lookup each path segment, and glue them all together.
    Path finalPath;
    finalPath.cost = 0;
    finalPath.dist = 0;
    for(int i=0; i<route.size()-1; ++i) {
        const Path& path = paths[route[i]][route[i+1]];
        finalPath = append(finalPath, path);
    }
    return finalPath;
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
    namespace cm = tinycolormap;
    
    // Parse the command line and populate the global config struct.
    parseCommandLine(argc, argv);

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
            savePFM(map, config.outputDir+"slope.pfm");
        }
        {
            auto map = slopeAtlas.north;
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeN.exr");
            savePFM(map, config.outputDir+"slopeN.pfm");
        }
        {
            auto map = slopeAtlas.east;
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeE.exr");
            savePFM(map, config.outputDir+"slopeE.pfm");
        }
        {
            auto map = slopeAtlas.northeast;
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeNE.exr");
            savePFM(map, config.outputDir+"slopeNE.pfm");
        }
        {
            auto map = slopeAtlas.southeast;
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeSE.exr");
            savePFM(map, config.outputDir+"slopeSE.pfm");
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
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeS.exr");
            savePFM(map, config.outputDir+"slopeS.pfm");
        }
        {
            auto map = negate(slopeAtlas.east);
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeW.exr");
            savePFM(map, config.outputDir+"slopeW.pfm");
        }
        {
            auto map = negate(slopeAtlas.northeast);
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeSW.exr");
            savePFM(map, config.outputDir+"slopeSW.pfm");
        }
        {
            auto map = negate(slopeAtlas.southeast);
            //drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"slopeNW.exr");
            savePFM(map, config.outputDir+"slopeNW.pfm");
        }


        {
            auto map = priorityMap;
            drawCircle(map, landingSiteX, landingSiteY, 10, 3.0);
            saveEXR(map, config.outputDir+"priority.exr");
            savePFM(map, config.outputDir+"priority.pfm");
        }
        {
            auto map = commsMap;
            drawCircle(map, landingSiteX, landingSiteY, 100, 3.0);
            saveEXR(map, config.outputDir+"comms.exr");
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
    //const auto paths2 = planAllPairsSLOW(allSites, slopeAtlas, commsMap);
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
        file << fmt::format("x,y,z,i,j,d\n");
        for(const auto& p : path.states) {
            Vantage v;
            v.x = elevationMap.j2x(p.j);
            v.y = elevationMap.j2x(p.i);
            v.z = elevationMap(p.i, p.j) + config.roverHeight;
            file << fmt::format("{},{},{},{},{},{}\n", v.x, v.y, v.z, p.i, p.j, (int)p.d);
        }
        file.close();
    }

    // Animate the final path.
    {
        // Create a directory for animation frames.
        std::string animDir = config.outputDir + "/anim";
        if( !std::filesystem::exists(animDir) ) {
            std::filesystem::create_directory(animDir);
        }

        int upsample = 4;
        ImageRGB baseImg(slopeAtlas.absolute, cm::ColormapType::Plasma, upsample);

        #pragma omp parallel for
        for(int f=0; f<path.states.size(); ++f) {
            fmt::print("Animating frame {}...\n", f);
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
