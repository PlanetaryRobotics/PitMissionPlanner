#include "path.h"
#include <cassert>

Path assembleRoute(const std::vector<int> &route, const std::vector<std::vector<Path>> paths) {
    // Walk the route, lookup each path segment, and glue them all together.
    Path finalPath;
    finalPath.cost = 0;
    finalPath.dist = 0;
    for (int i = 0; i < route.size() - 1; ++i) {
        const Path &path = paths[route[i]][route[i + 1]];
        finalPath = append(finalPath, path);
    }
    if( config.returnToStart ) {
        const Path &path = paths[route[route.size()-1]][route[0]];
        finalPath = append(finalPath, path);
    }
    return finalPath;
}
Path append(const Path &a, const Path &b) {
    assert(a.cost >= 0);
    assert(a.dist >= 0);
    assert(b.cost >= 0);
    assert(b.dist >= 0);
    Path c;
    c.cost = a.cost + b.cost;
    c.dist = a.dist + b.dist;
    c.states = a.states;
    c.states.insert(c.states.end(), b.states.begin(), b.states.end());
    return c;
}

Path reverse(const Path &path) {
    Path rPath = path;
    std::reverse(rPath.states.begin(), rPath.states.end());
    return rPath;
}
uint8_t setBit(uint8_t byte, uint8_t bit) { return byte | ((uint8_t)1) << bit; };
uint8_t clearBit(uint8_t byte, uint8_t bit) { return byte & ~(((uint8_t)1) << bit); };
bool checkBit(uint8_t byte, uint8_t bit) { return (byte >> bit) & 0x01; };
bool checkMapBit(const TerrainMapU8 &map, const Path::State &state) { return checkBit(map(state.i, state.j), (int)state.d); };
void setMapBit(TerrainMapU8 &map, const Path::State &state) { map(state.i, state.j) = setBit(map(state.i, state.j), (int)state.d); };

std::vector<Path::State> getSuccessors(const Path::State &p, const SlopeAtlas &slopeAtlas, const TerrainMapFloat &commsMap) {

    std::vector<Path::State> succs;
    succs.reserve(7 + 2);

    // You can only point turn if the terrain is flat enough.
    if (slopeAtlas.absolute(p.i, p.j) <= config.roverPointTurnSlopeLimit) {
        for (int dd = -1; dd <= 1; ++dd) {
            Path::State s;
            s.i = p.i;
            s.j = p.j;
            s.d = static_cast<Direction>((static_cast<int>(p.d) + dd + 8) % 8);
            if (s.d != p.d) {
                succs.push_back(s);
            }
        }
    }

    double lonSlope = std::abs(computeLongitudinalSlope(p, slopeAtlas));
    double latSlope = std::abs(computeLateralSlope(p, slopeAtlas));

    // If the lateral and longitudinal slopes are safe,
    // try moving one step forward or one step backward.
    if (lonSlope <= config.roverLongitudinalSlopeLimit && latSlope <= config.roverLateralSlopeLimit) {

        Path::State sf = p;
        Path::State sb = p;

        constexpr int dI[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
        constexpr int dJ[8] = {0, 1, 1, 1, 0, -1, -1, -1};

        sf.i += dI[(int)p.d];
        sf.j += dJ[(int)p.d];

        sb.i -= dI[(int)p.d];
        sb.j -= dJ[(int)p.d];

        double fwdLonSlope = std::abs(computeLongitudinalSlope(sf, slopeAtlas));
        double fwdLatSlope = std::abs(computeLateralSlope(sf, slopeAtlas));
        double bwdLonSlope = std::abs(computeLongitudinalSlope(sb, slopeAtlas));
        double bwdLatSlope = std::abs(computeLateralSlope(sb, slopeAtlas));

        auto inBounds = [&commsMap](int i, int j) -> bool { return (0 <= i && i < commsMap.rows && 0 <= j && j < commsMap.cols); };

        if (inBounds(sf.i, sf.j) && commsMap(sf.i, sf.j) <= config.roverRangeFromComms && fwdLonSlope <= config.roverLongitudinalSlopeLimit &&
            fwdLatSlope <= config.roverLateralSlopeLimit) {
            succs.push_back(sf);
        }
        if (inBounds(sb.i, sb.j) && commsMap(sb.i, sb.j) <= config.roverRangeFromComms && bwdLonSlope <= config.roverLongitudinalSlopeLimit &&
            bwdLatSlope <= config.roverLateralSlopeLimit) {
            succs.push_back(sb);
        }
    }
    return succs;
}

double computeLateralSlope(const Path::State &s, const SlopeAtlas &slopeAtlas) {
    double latSlope = 90.0;
    switch (s.d) {
    case Direction::N:
        latSlope = slopeAtlas.east(s.i, s.j);
        break;
    case Direction::NE:
        latSlope = slopeAtlas.southeast(s.i, s.j);
        break;
    case Direction::E:
        latSlope = -slopeAtlas.north(s.i, s.j);
        break;
    case Direction::SE:
        latSlope = -slopeAtlas.northeast(s.i, s.j);
        break;
    case Direction::S:
        latSlope = -slopeAtlas.east(s.i, s.j);
        break;
    case Direction::SW:
        latSlope = -slopeAtlas.southeast(s.i, s.j);
        break;
    case Direction::W:
        latSlope = slopeAtlas.north(s.i, s.j);
        break;
    case Direction::NW:
        latSlope = slopeAtlas.northeast(s.i, s.j);
        break;
    default:
        break;
    }
    return latSlope;
}

double computeLongitudinalSlope(const Path::State &s, const SlopeAtlas &slopeAtlas) {
    double lonSlope = 90.0;
    switch (s.d) {
    case Direction::N:
        lonSlope = slopeAtlas.north(s.i, s.j);
        break;
    case Direction::NE:
        lonSlope = slopeAtlas.northeast(s.i, s.j);
        break;
    case Direction::E:
        lonSlope = slopeAtlas.east(s.i, s.j);
        break;
    case Direction::SE:
        lonSlope = slopeAtlas.southeast(s.i, s.j);
        break;
    case Direction::S:
        lonSlope = -slopeAtlas.north(s.i, s.j);
        break;
    case Direction::SW:
        lonSlope = -slopeAtlas.northeast(s.i, s.j);
        break;
    case Direction::W:
        lonSlope = -slopeAtlas.east(s.i, s.j);
        break;
    case Direction::NW:
        lonSlope = -slopeAtlas.southeast(s.i, s.j);
        break;
    default:
        break;
    }
    return lonSlope;
}

double directionToDegrees(const Direction &d) {
    const double dir2angle[] = {90, 45, 0, 315, 270, 225, 180, 135};
    return dir2angle[(int)d];
}

std::string directionToString(const Direction &d) {
    const std::string dir2name[] = {"N", "NE", "E", "SE", "S", "SW", "W", "NW"};
    return dir2name[(int)d];
}
double octileDistance(const Path::State &a, const Path::State &b) {
    // Branchless octile distance
    int di = std::abs(a.i - b.i);
    int dj = std::abs(a.j - b.j);
    constexpr double diagonal = std::sqrt(2);
    constexpr double twoCardinalMinusDiagonal = 2 - diagonal;
    return (twoCardinalMinusDiagonal * abs(di - dj) + diagonal * (di + dj)) / 2;
}

Direction directionFromTo(const Path::State &from, const Path::State &to) {
    constexpr double tan_pi_8 = std::tan(M_PI / 8.0);

    const int di = to.i - from.i;
    const int dj = to.j - from.j;

    if (di <= 0 and std::abs(dj) <= tan_pi_8 * -di) {
        return Direction::N;
    }
    if (di >= 0 and std::abs(dj) <= tan_pi_8 * di) {
        return Direction::S;
    }

    if (dj >= 0 and std::abs(di) <= tan_pi_8 * dj) {
        return Direction::E;
    }
    if (dj <= 0 and std::abs(di) <= tan_pi_8 * -dj) {
        return Direction::W;
    }

    if (di < 0 and dj > 0) {
        return Direction::NE;
    }
    if (di > 0 and dj > 0) {
        return Direction::SE;
    }

    if (di < 0 and dj < 0) {
        return Direction::NW;
    }
    if (di > 0 and dj < 0) {
        return Direction::SW;
    }

    return to.d;
};
