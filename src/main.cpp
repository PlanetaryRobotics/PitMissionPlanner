#include "terrainmesh.h"
#include "terrainmap.h"
#include "argh.h"
#include "flat_hash_map.hpp"
#include "priority_queue.h"
#include <Eigen/Dense>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <optional>
#include <unordered_set>
#include <iostream>
#include <filesystem>

#include <chrono>

struct PlannerConfiguration {
    std::string meshfile   = "../meshes/lmp_crop.ply";
    std::string outputDir  = "./";
    double landingSiteX    = 100;       // meters
    double landingSiteY    = 200;       // meters
    double landerHeight    = 2.0;       // meters

    double mapPitch             = 1.0;       // meters
    int    numProbes            = 10000;     // 
    int    numCandidates        = 10000;     // 
    int    numVantages          = 15;        // 
    double visAngle             = 55;        // degrees
    double maxVisRange          = 300;       // meters
    double minVantageSeparation = 20.0;      // meters

    double roverHeight                 =  1.0; // meters
    double roverSpeed                  = 0.01; // m/s
    double roverFOV                    =  180; // degrees
    double roverLateralSlopeLimit      = 12.0; // degrees
    double roverLongitudinalSlopeLimit = 20.0; // degrees
    double roverPointTurnSlopeLimit    =  5.0; // degrees

    const double distCostMultiplier  = 1.0;
    const double slopeCostMultiplier = 8.0;
    const double turnCostMultiplier  = 0.1;
    const double heuristicMultiplier = 1.0;
};

// Global configuration structure
PlannerConfiguration config;

enum class Direction : int8_t { N = 0, NE, E, SE, S, SW, W, NW };

double directionToDegrees(const Direction& d) {
    const double dir2angle[] = {90, 45, 0, 315, 270, 225, 180, 135};
    return dir2angle[(int)d];
}

std::string directionToString(const Direction& d) {
    const std::string dir2name[] = {"N", "NE", "E", "SE", "S", "SW", "W", "NW"};
    return dir2name[(int)d];
}


void parseCommandLine(int argc, char* argv[]) {
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
        return;
    }

    cmdl(1, config.meshfile) >> config.meshfile;
    cmdl(2, config.outputDir) >> config.outputDir;
    config.outputDir = config.outputDir + "/";

    cmdl("p", config.mapPitch) >> config.mapPitch;
    cmdl("x", config.landingSiteX) >> config.landingSiteX;
    cmdl("y", config.landingSiteY) >> config.landingSiteY;
    cmdl("lonslope", config.roverLongitudinalSlopeLimit) >> config.roverLongitudinalSlopeLimit;
    cmdl("latslope", config.roverLateralSlopeLimit) >> config.roverLateralSlopeLimit;
    cmdl("turnslope", config.roverPointTurnSlopeLimit) >> config.roverPointTurnSlopeLimit;
    cmdl("rh", config.roverHeight) >> config.roverHeight;
    cmdl("lh", config.landerHeight) >> config.landerHeight;
    cmdl("np", config.numProbes) >> config.numProbes;
    cmdl("nc", config.numCandidates) >> config.numCandidates;
    cmdl("nv", config.numVantages) >> config.numVantages;
    cmdl("va", config.visAngle) >> config.visAngle;
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

    double gridIndexToXCoord(size_t j) const { return absolute.gridIndexToXCoord(j); }
    double gridIndexToYCoord(size_t i) const { return absolute.gridIndexToYCoord(i); }

    size_t xCoordToGridIndex(double x) const { return absolute.xCoordToGridIndex(x); }
    size_t yCoordToGridIndex(double y) const { return absolute.yCoordToGridIndex(y); }

    double width() const { return absolute.width(); }
    double height() const { return absolute.height(); }

    size_t rows() const { return absolute.rows; }
    size_t cols() const { return absolute.cols; }
    size_t pitch() const { return absolute.pitch; }
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

    double cost = -1.0;
    double dist = -1.0;
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
    return checkBit(map.atIJ(state.i, state.j), (int)state.d);
};
void setMapBit(TerrainMapU8& map, const Path::State& state) {
    map.atIJ(state.i, state.j) = setBit(map.atIJ(state.i, state.j), (int)state.d);
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
            ray.oX = slopeAtlas.gridIndexToXCoord(j);
            ray.oY = slopeAtlas.gridIndexToYCoord(i);
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
    const int ROWS = commsMap.rows;
    const int COLS = commsMap.cols;

    TerrainMapU8 reachMap(commsMap.width(), commsMap.height(), commsMap.pitch);

    auto inBounds = [ROWS, COLS](int i, int j) -> bool {
        return (0 <= i && i<ROWS && 0 <= j && j<COLS);
    };

    // Generate the successors of a node.
    auto getSuccessors = [&inBounds, &slopeAtlas, &commsMap](const Path::State& p) -> std::vector<Path::State> {
        std::vector<Path::State> succs; succs.reserve(7+2);
        // You can only point turn if the terrain is flat enough.
        if( slopeAtlas.absolute(p.i, p.j) <= config.roverPointTurnSlopeLimit ) {
            for(int dd=0; dd<8; ++dd) {
                Path::State s;
                s.i=p.i;
                s.j=p.j;
                s.d=static_cast<Direction>((static_cast<int>(p.d)+dd)%8);
                if( s.d != p.d ) {
                    succs.push_back(s);
                }
            }
        }

        // Try moving one step forward or one step backward.
        Path::State sf = p;
        Path::State sb = p;

        constexpr int dX[] = {0,1,1,1,0,-1,-1,-1};
        constexpr int dY[] = {-1,-1,0,1,1,1,0,-1};

        sf.i += dY[(int)p.d];
        sf.j += dX[(int)p.d];

        sb.i -= dY[(int)p.d];
        sb.j -= dX[(int)p.d];

        // Depending on direction, figure out the lateral and longitudinal slope at the current state, p.
        // NOTE(Jordan): This is a good opportunity to optimize...
        double nSlopeFwd  = std::abs(slopeAtlas.north(sf.i, sf.j));
        double eSlopeFwd  = std::abs(slopeAtlas.east(sf.i, sf.j));
        double neSlopeFwd = std::abs(slopeAtlas.northeast(sf.i, sf.j));
        double seSlopeFwd = std::abs(slopeAtlas.southeast(sf.i, sf.j));

        double nSlopeBwd  = std::abs(slopeAtlas.north(sb.i, sb.j));
        double eSlopeBwd  = std::abs(slopeAtlas.east(sb.i, sb.j));
        double neSlopeBwd = std::abs(slopeAtlas.northeast(sb.i, sb.j));
        double seSlopeBwd = std::abs(slopeAtlas.southeast(sb.i, sb.j));

        double fwdLonSlope = 90.0;
        double fwdLatSlope = 90.0;

        double bwdLonSlope = 90.0;
        double bwdLatSlope = 90.0;

        switch (p.d) {
            case Direction::N :
                fwdLonSlope = nSlopeFwd;
                fwdLatSlope = eSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::NE :
                fwdLonSlope = neSlopeFwd;
                fwdLatSlope = seSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::E :
                fwdLonSlope = eSlopeFwd;
                fwdLatSlope = nSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::SE :
                fwdLonSlope = seSlopeFwd;
                fwdLatSlope = neSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::S :
                fwdLonSlope = nSlopeFwd;
                fwdLatSlope = eSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::SW :
                fwdLonSlope = neSlopeFwd;
                fwdLatSlope = seSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::W :
                fwdLonSlope = eSlopeFwd;
                fwdLatSlope = nSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::NW :
                fwdLonSlope = seSlopeFwd;
                fwdLatSlope = neSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            default:
                break;
        }

        // You can only move forward or backward if the next state is reachable and
        // doing so does not violate the longitudinal or lateral slope limits.
        if( inBounds(sf.i, sf.j) &&
            commsMap(sf.i, sf.j) &&
            fwdLatSlope <= config.roverLateralSlopeLimit &&
            fwdLonSlope <= config.roverLongitudinalSlopeLimit ) {
                succs.push_back(sf);
        }
        if( inBounds(sb.i, sb.j) &&
            commsMap(sb.i, sb.j) &&
            bwdLatSlope <= config.roverLateralSlopeLimit &&
            bwdLonSlope <= config.roverLongitudinalSlopeLimit ) {
                succs.push_back(sb);
        }

        return succs; 
    };

    TerrainMapU8 visited(commsMap.width(), commsMap.height(), commsMap.pitch);

    std::vector<Path::State> open;
    open.push_back(landingSite);

    int iterations = 0;
    while( !open.empty() ) {
        if( iterations++ % (1<<18) ) {
            fmt::print("[{}] Building reach map.\n", iterations);
        }

        // Pop currState from open.
        const auto currState = open.back();
        open.pop_back();

        // If currState is visited, continue (skip).
        if( checkMapBit(visited, currState) ) { continue; }

        // currState is visited.
        setMapBit(visited, currState);

        // currState is reachable.
        setMapBit(reachMap, currState);

        // Loop over successors of current state.
        const auto& successors = getSuccessors(currState);
        for(const auto& succ : successors) {
            // mark the successor as reachable
            setMapBit(reachMap, succ);
            // add the successor to the open queue.
            open.push_back(succ);
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
    const double landerX = elevationMap.gridIndexToXCoord(lj);
    const double landerY = elevationMap.gridIndexToYCoord(li);
    const double landerZ = elevationMap(li, lj) + config.landerHeight;

    for(int ri=0; ri<commsMap.rows; ++ri) {
        for(int rj=0; rj<commsMap.cols; ++rj) {
            TerrainMesh::Ray ray;
            const double roverX = elevationMap.gridIndexToXCoord(rj);
            const double roverY = elevationMap.gridIndexToYCoord(ri);
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
            
            // FIXME(Jordan): Right now this check requires a vantage
            // to be reachable facing any direction. This is too conservative,
            // but for some reason I can't reliably plan to states
            // that can't be reached in some directions.
            //
            // The condition here should be this:
            // if( checkMapBit(reachMap, state) ) {

            if( reachMap(state.i, state.j) == 0b11111111 ) {
                Vantage v;
                v.x = reachMap.gridIndexToXCoord(state.j);
                v.y = reachMap.gridIndexToYCoord(state.i);
                v.z = elevationMap(state.i,state.j) + config.roverHeight;
                v.dir = state.d;
                candidates.push_back(v);
            }
        }
    }

    TerrainMapFloat candidateMap(reachMap.width(), reachMap.height(), reachMap.pitch);
    for(int i=0; i<candidateMap.rows; ++i) {
        for(int j=0; j<candidateMap.cols; ++j) {
            candidateMap(i,j) = (reachMap(i,j) == 0) ? 0.0f : 1.0f;
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
                if( hitAngle < config.visAngle && hitDist < config.maxVisRange && std::abs(rayNorm-hitDist) < 0.05*rayNorm ) {
                    candidate.coverage[pi] = true;
                    candidate.totalCoverage += p.priority;
                }
            }
        }
        int i = candidateMap.yCoordToGridIndex(candidate.y);
        int j = candidateMap.xCoordToGridIndex(candidate.x);
        candidateMap(i,j) = std::max<float>(candidate.totalCoverage, candidateMap(i,j));
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
                double d2 = (t.x-c.x)*(t.x-c.x)+(t.y-c.y)*(t.y-c.y)+(t.z-c.z);
                if( d2 < config.minVantageSeparation*config.minVantageSeparation ) {
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

TerrainMapFloat buildCoverageMap(const TerrainMesh& tmesh, const TerrainMapFloat& elevationMap,
                            const std::vector<Vantage>& vantages) {

    TerrainMapFloat coverageMap(elevationMap.width(), elevationMap.height(), elevationMap.pitch);

    for(int vi=0; vi<vantages.size(); ++vi) {
        auto& v = vantages[vi];
        fmt::print("[{}/{}] Generating coverage maps.\n", vi, vantages.size());
        #pragma omp parallel for
        for(int ri=0; ri<coverageMap.rows; ++ri) {
            for(int rj=0; rj<coverageMap.cols; ++rj) {
                const double terrainX = elevationMap.gridIndexToXCoord(rj);
                const double terrainY = elevationMap.gridIndexToYCoord(ri);
                const double terrainZ = elevationMap(ri, rj);

                // Is this ray in front of the rover and within its field of view?
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
                    if( hitAngle < config.visAngle && hitDist < config.maxVisRange && std::abs(rayNorm-hitDist) < rayNorm*0.05 ) { 
                        coverageMap(ri, rj)++;
                    }
                } 
            }
        }
    }
    for(const auto& v : vantages) {
        drawCircle(coverageMap, v.x, v.y, (float)2*vantages.size(), 2.0);
    }
    return coverageMap;
}

std::vector<Path> multipathplan(const SlopeAtlas& slopeAtlas,
                                const TerrainMapFloat& reachMap,
                                const Path::State& start,
                                const std::vector<Path::State>& goals) {
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
        const double turn = 0.25 * std::abs(std::abs(static_cast<int>(a.d)-static_cast<int>(goalDir))-4);
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
        //fmt::print("({},{},{}) -> ({},{},{}) ... {:0.3f} {:0.3f} {:0.3f}\n", a.i, a.j, a.d, b.i, b.j, b.d, slope, dist, turn);
        return config.slopeCostMultiplier * (slope/90.0) +
               config.distCostMultiplier * dist +
               config.turnCostMultiplier * turn;
    };

    auto inBounds = [ROWS, COLS](int i, int j) -> bool {
        return (0 <= i && i<ROWS && 0 <= j && j<COLS);
    };

    // Generate the successors of a node.
    auto getSuccessors = [&inBounds, &slopeAtlas, &reachMap](const Path::State& p) -> std::vector<Path::State> {
        std::vector<Path::State> succs; succs.reserve(7+2);
        // You can only point turn if the terrain is flat enough.
        if( slopeAtlas.absolute(p.i, p.j) <= config.roverPointTurnSlopeLimit ) {
            for(int dd=0; dd<8; ++dd) {
                Path::State s;
                s.i=p.i;
                s.j=p.j;
                s.d=static_cast<Direction>((static_cast<int>(p.d)+dd)%8);
                if( s.d != p.d ) {
                    succs.push_back(s);
                }
            }
        }

        // Try moving one step forward or one step backward.
        Path::State sf = p;
        Path::State sb = p;

        constexpr int dX[] = {0,1,1,1,0,-1,-1,-1};
        constexpr int dY[] = {-1,-1,0,1,1,1,0,-1};

        sf.i += dY[(int)p.d];
        sf.j += dX[(int)p.d];

        sb.i -= dY[(int)p.d];
        sb.j -= dX[(int)p.d];

        // Depending on direction, figure out the lateral and longitudinal slope at the current state, p.
        // NOTE(Jordan): This is a good opportunity to optimize...
        double nSlopeFwd  = std::abs(slopeAtlas.north(sf.i, sf.j));
        double eSlopeFwd  = std::abs(slopeAtlas.east(sf.i, sf.j));
        double neSlopeFwd = std::abs(slopeAtlas.northeast(sf.i, sf.j));
        double seSlopeFwd = std::abs(slopeAtlas.southeast(sf.i, sf.j));

        double nSlopeBwd  = std::abs(slopeAtlas.north(sb.i, sb.j));
        double eSlopeBwd  = std::abs(slopeAtlas.east(sb.i, sb.j));
        double neSlopeBwd = std::abs(slopeAtlas.northeast(sb.i, sb.j));
        double seSlopeBwd = std::abs(slopeAtlas.southeast(sb.i, sb.j));

        double fwdLonSlope = 90.0;
        double fwdLatSlope = 90.0;

        double bwdLonSlope = 90.0;
        double bwdLatSlope = 90.0;

        switch (p.d) {
            case Direction::N :
                fwdLonSlope = nSlopeFwd;
                fwdLatSlope = eSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::NE :
                fwdLonSlope = neSlopeFwd;
                fwdLatSlope = seSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::E :
                fwdLonSlope = eSlopeFwd;
                fwdLatSlope = nSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::SE :
                fwdLonSlope = seSlopeFwd;
                fwdLatSlope = neSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::S :
                fwdLonSlope = nSlopeFwd;
                fwdLatSlope = eSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::SW :
                fwdLonSlope = neSlopeFwd;
                fwdLatSlope = seSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::W :
                fwdLonSlope = eSlopeFwd;
                fwdLatSlope = nSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            case Direction::NW :
                fwdLonSlope = seSlopeFwd;
                fwdLatSlope = neSlopeFwd;
                bwdLonSlope = nSlopeBwd;
                bwdLatSlope = eSlopeBwd;
                break;
            default:
                break;
        }

        // You can only move forward or backward if the next state is reachable and
        // doing so does not violate the longitudinal or lateral slope limits.
        if( inBounds(sf.i, sf.j) &&
            reachMap(sf.i, sf.j) &&
            fwdLatSlope <= config.roverLateralSlopeLimit &&
            fwdLonSlope <= config.roverLongitudinalSlopeLimit ) {
                succs.push_back(sf);
        }
        if( inBounds(sb.i, sb.j) &&
            reachMap(sb.i, sb.j) &&
            bwdLatSlope <= config.roverLateralSlopeLimit &&
            bwdLonSlope <= config.roverLongitudinalSlopeLimit ) {
                succs.push_back(sb);
        }

        return succs; 
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
        if( open.empty() ) { break; }

        const Path::State currState = open.top();
        open.pop();

        // If you have expanded this state before, skip it.
        if( closed.find(currState) != closed.end() ) { continue; }

        // Otherwise, add it to closed and expand it.
        closed.insert(currState);

        if( ++expansions % (1<<16) == 0 ) {
            fmt::print("Expansions: {}\n", expansions);
        }

        // If this state is a goal, remove it from the goal set!
        auto it = goalSet.find(currState);
        if( it != goalSet.end() ) {
            fmt::print("Found a path from ({},{}) to ({},{})! {} goals remaining...\n", start.i, start.j, it->i, it->j, goalSet.size()-1);
            goalSet.erase(it);
            if( goalSet.empty() ) {
                break;
            }
        }

        // Since we can't end the search yet, let's look
        // up some info on the current point and keep working.
        const double currG = nodeMap.at(currState).gscore;

        // Examine each successor of currState.
        const auto succs = getSuccessors(currState);
        for(const auto& succ : succs) {
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
            path.dist += std::sqrt((pred.i-state.i)*(pred.i-state.i)+(pred.j-state.j)*(pred.j-state.j));
            state = pred;
        }
        path.states.push_back(start);
        path.cost = nodeMap[goal].gscore;

        path = reverse(path);
        paths.push_back(path);
    }
    return paths;
}

std::vector<Path::State> routeplan(const std::vector<Path::State> allSites, const Eigen::MatrixXd& costs) {
    std::vector<int> routeIndices;
    std::vector<bool> visited(allSites.size(), false);
    int visitedCount = 1;

    int currentIdx = 0;
    visited[currentIdx] = true;
    routeIndices.push_back(currentIdx);
    double routeCost = 0;

    while(visitedCount < allSites.size()) {
        int minCostIdx = -1;
        double minCost = std::numeric_limits<double>::infinity();
        for(int i=0; i<costs.rows(); ++i) {
            const auto cost = costs(currentIdx, i);
            if( !visited[i] && cost >= 0 && cost < minCost ) {
                minCostIdx = i;
                minCost = cost;
            }
        }

        // If there are no valid routes to another city, give up.
        if( minCostIdx < 0 ) { break; }

        currentIdx = minCostIdx;
        routeIndices.push_back(minCostIdx);
        visited[minCostIdx] = true;
        routeCost += minCost;
        visitedCount++;
    }

    // Two-Opt Route Improvement
    int iterations = 0;
    bool improved = false;
    do {
        fmt::print("Improving route [{}] ...\n", iterations);
        fmt::print("Route cost: {}\n", routeCost);
        improved = false;
        for(int i=1; i<routeIndices.size()-1; ++i) {

            // Don't rewire the segment 0<->1 because that segment
            // connects the landing site to the closest vantage.
            if( i == 1 ) { continue; }

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
                // If any segment is invalid (dist < 0), forget this swap.
                double dist = 0;
                bool valid = true;
                for(int i=0; i<newIndices.size()-1; ++i) {
                    double d = costs(newIndices[i], newIndices[i+1]);
                    if( d < 0 ) { valid = false; break; }
                    dist += d;
                }
                if( !valid ) { continue; }

                // This route is shorter! Keep it.
                if( dist < routeCost ) {
                    fmt::print("[{}] Swap {:2}<->{:2} Length: {}\n", iterations, i,j, dist);
                    routeCost = dist;
                    routeIndices = newIndices;
                    improved = true;
                }
            }
        }
    } while( improved && iterations++ < 100 );

    std::vector<Path::State> route;
    for(const auto& ri : routeIndices) { route.push_back(allSites[ri]); }

    return route;
}

std::vector<std::vector<Path>> planAllPairs(const std::vector<Path::State>& sites,
                                            const SlopeAtlas& slopeAtlas,
                                            const TerrainMapFloat& reachMap) {
    std::vector<std::vector<Path>> allPaths;
    allPaths.resize(sites.size());
    for(auto& pathList : allPaths) { pathList.resize(sites.size()); }

    // Generate all combinations of two indices into the allSites vector.
    #pragma omp parallel for
    for(int a=0; a<sites.size()-1; ++a) {
        #pragma omp critical(PRINT)
        {
            fmt::print("Planning {} paths from site {:2}.\n", sites.size()-a-1, a);
        }

        const Path::State start = sites[a];

        std::vector<Path::State> goals;
        for(int b=a+1; b<sites.size(); ++b) {
            goals.push_back(sites[b]);
        }

        auto paths = multipathplan(slopeAtlas, reachMap, start, goals);

        #pragma omp critical
        {
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
    }

    // Draw all paths on their own map.
    for(int a=0; a<allPaths.size(); ++a) {
        for(int b=0; b<allPaths[0].size(); ++b) {
            if( a >= b ) { continue; }
            const auto& path = allPaths[a][b];

            TerrainMapFloat pathMap = slopeAtlas.absolute;

            // Draw the start.
            const auto& start = sites[a];
            double sx = pathMap.gridIndexToXCoord(start.j);
            double sy = pathMap.gridIndexToYCoord(start.i);
            drawCircle(pathMap, sx,sy, 100, 2.0);

            // Draw the goal.
            const auto&  goal = sites[b];
            double gx = pathMap.gridIndexToXCoord(goal.j);
            double gy = pathMap.gridIndexToYCoord(goal.i);
            drawCircle(pathMap, gx, gy, 100, 2.0);

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
            double sx = allPathsMap.gridIndexToXCoord(start.j);
            double sy = allPathsMap.gridIndexToYCoord(start.i);
            drawCircle(allPathsMap, sx,sy, 100, 2.0);

            // Draw the goal.
            const auto& goal = sites[b];
            double gx = allPathsMap.gridIndexToXCoord(goal.j);
            double gy = allPathsMap.gridIndexToYCoord(goal.i);
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

Path assembleRoute(const std::vector<Path::State>& route, const std::vector<std::vector<Path>> paths, const TerrainMapFloat& elevationMap) {
    // Stick all the paths in a hash table indexed by (start, goal) pairs.
    using StatePair = std::pair<Path::State, Path::State>;
    auto hashStatePair = [&elevationMap](const StatePair& pp) {
        std::hash<Path::State> hasher;
        auto hash = hasher(pp.first);
        hash ^= hasher(pp.second) + 0x9e3779b9 + (hash<<6) + (hash>>2);
        return hash;
    };
    ska::flat_hash_map<StatePair, Path, decltype(hashStatePair)> pathLookup(10, hashStatePair);

    for(int a=0; a<paths.size(); ++a) {
        for(int b=0; b<paths[0].size(); ++b) {
            const auto& p = paths[a][b];
            if( p.states.size() != 0 ) {
                Path::State s = p.states[0];
                Path::State g = p.states[p.states.size()-1];
                pathLookup[std::make_pair(s,g)] = p;
            }
        }
    }

    // Walk the route, lookup each path segment, and glue them all together.
    Path finalPath;
    for(int i=0; i<route.size()-1; ++i) {
        const auto& s0 = route[i]; 
        const auto& s1 = route[i+1]; 
        const auto key = std::make_pair(s0, s1);
        Path path = pathLookup.at(key);
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

    // Parse the command line and populate the global config struct.
    parseCommandLine(argc, argv);

    // Create the output directory if it doesn't exist already.
    if( !std::filesystem::exists(config.outputDir) ) {
        std::filesystem::create_directory(config.outputDir);
    }

    // Read the terrain mesh.
    TerrainMesh tmesh(config.meshfile);

    // Construct elevation, priority, and slope maps.
    auto [elevationMap, priorityMap, slopeAtlas] = buildTerrainMaps(tmesh, config.mapPitch);

    // Compute landing site coordinates.
    double landingSiteX = config.landingSiteX;
    double landingSiteY = config.landingSiteY;
    double landingSiteZ = elevationMap.atXY(landingSiteX, landingSiteY);

    int landingSiteI = elevationMap.yCoordToGridIndex(landingSiteY);
    int landingSiteJ = elevationMap.xCoordToGridIndex(landingSiteX);

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
        {
            auto map = reachMap;
            drawCircle(map, landingSiteX, landingSiteY, 100, 4.0);
            saveEXR(map, config.outputDir+"reach.exr");
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
        }
    }
    drawCircle(vantageMap, landingSiteX, landingSiteY, 100, 3.0);
    saveEXR(vantageMap, config.outputDir+"vantages.exr");

    // Draw separate coverage maps for each vantage.
    for(int vi=0; vi<vantages.size(); ++vi) {
        std::vector<Vantage> tmp;
        tmp.push_back(vantages[vi]);
        auto coverageMap = buildCoverageMap(tmesh, elevationMap, tmp);
        int j = vantageMap.xCoordToGridIndex(vantages[vi].x);
        int i = vantageMap.yCoordToGridIndex(vantages[vi].y);
        drawCircle(coverageMap, vantages[vi].x, vantages[vi].y, vantages.size()+1, 2.0);
        drawCircle(coverageMap, landingSiteX, landingSiteY, vantages.size()+10, 3.0);
        saveEXR(coverageMap, config.outputDir+fmt::format("coverage_{:02}.exr",vi));
    }

    // Construct a state for the landing site.
    // Set its direction to point toward the centroid of the vantages.
    std::vector<Path::State> allSites;
    {
        Path::State centroid;
        double cx = 0; double cy = 0;
        for(const auto& v : vantages) {
            cx += v.x;
            cy += v.y;
        }
        cx /= vantages.size(); cy /= vantages.size();
        centroid.i = slopeAtlas.yCoordToGridIndex(cy);
        centroid.j = slopeAtlas.xCoordToGridIndex(cx);
        centroid.d = Direction::N; // Don't care.

        allSites.push_back(landingState);

        for(const auto& v : vantages) {
            Path::State s;
            s.i = slopeAtlas.yCoordToGridIndex(v.y);
            s.j = slopeAtlas.xCoordToGridIndex(v.x);
            s.d = directionFromTo(s, centroid);
            allSites.push_back(s);
        }
    }

    // Save the sites to a csv file.
    {
        std::ofstream file;
        file.open(config.outputDir+"sites.csv");
        file << fmt::format("x,y,z,i,j,d\n");
        for(const auto& s : allSites) {
            Vantage v;
            v.x = elevationMap.gridIndexToXCoord(s.j);
            v.y = elevationMap.gridIndexToXCoord(s.i);
            v.z = elevationMap(s.i, s.j) + config.roverHeight;
            file << fmt::format("{:0.3f},{:0.3f},{:0.3f}, {},{},{}\n", v.x, v.y, v.z, s.i, s.j, directionToString(s.d));
        }
        file.close();
    }

    // Compute paths between all pairs of k vantages plus the landing site.
    const auto paths = planAllPairs(allSites, slopeAtlas, commsMap);

    // Organize path costs into convenient matrices for route planning.
    Eigen::MatrixXd costs(allSites.size(), allSites.size()); costs.fill(-1);
    Eigen::MatrixXd dists(allSites.size(), allSites.size()); dists.fill(-1);
    for(int a=0; a<paths.size(); ++a) {
        for(int b=0; b<paths[0].size(); ++b) {
            costs(a,b) = paths[a][b].cost;
            dists(a,b) = paths[a][b].dist;
        }
    }
    fmt::print("\nCosts: \n{}\n\n", costs);
    fmt::print("Dists: \n{}\n\n", dists);

    // Compute exploration route.
    const auto route = routeplan(allSites, costs);

    if( route.size() < allSites.size() ) {
        fmt::print("Oh no! I failed to plan a route to all vantages.\n");
    }

    // Chain paths together to create the final route.
    Path path = assembleRoute(route, paths, elevationMap);

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
            v.x = elevationMap.gridIndexToXCoord(p.j);
            v.y = elevationMap.gridIndexToXCoord(p.i);
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
            v.x = elevationMap.gridIndexToXCoord(p.j);
            v.y = elevationMap.gridIndexToXCoord(p.i);
            v.z = elevationMap(p.i, p.j) + config.roverHeight;
            file << fmt::format("{},{},{},{},{},{}\n", v.x, v.y, v.z, p.i, p.j, (int)p.d);
        }
        file.close();
    }

    return 0;
}
