#include "maps.h"
#include <tuple>
//#include "terrainmap.h"
//#include "terrainmesh.h"
#include "priority_queue.h"
std::tuple<TerrainMapFloat, TerrainMapFloat, SlopeAtlas> buildTerrainMaps(const TerrainMesh &tmesh, const double mapPitch) {
    const double mapX = tmesh.maxX() - tmesh.minX();
    const double mapY = tmesh.maxY() - tmesh.minY();

    TerrainMapFloat elevationMap(mapX, mapY, mapPitch);
    TerrainMapFloat priorityMap(mapX, mapY, mapPitch);
    SlopeAtlas slopeAtlas(mapX, mapY, mapPitch);

    const double maxZ = tmesh.maxZ();
    const double minZ = tmesh.minZ();

    for (int i = 0; i < slopeAtlas.rows(); ++i) {
        fmt::print("[{}/{}] Building terrain maps.\n", i, slopeAtlas.rows());
#pragma omp parallel for
        for (int j = 0; j < slopeAtlas.cols(); ++j) {
            TerrainMesh::Ray ray;
            ray.oX = slopeAtlas.j2x(j);
            ray.oY = slopeAtlas.i2y(i);
            ray.oZ = maxZ + 10.0;
            ray.dX = 0.0;
            ray.dY = 0.0;
            ray.dZ = -1.0;
            const auto hit = tmesh.raytrace(ray);
            if (hit) {
                elevationMap(i, j) = hit->z;
                priorityMap(i, j) = hit->priority;

                double slope = 180.0 / M_PI * std::acos(0.0 * hit->nx + 0.0 * hit->ny + 1.0 * hit->nz);
                if (std::isnan(slope)) {
                    slope = 0.0f;
                }
                slopeAtlas.absolute(i, j) = slope;

                // Project the normal onto the YZ plane and normalize it.
                // Then compute its angle with the +Y axis.
                {
                    // Set x to zero to project onto the YZ plane.
                    double nx_on_yz = 0.0;
                    double ny_on_yz = hit->ny;
                    double nz_on_yz = hit->nz;
                    double nz_on_yz_norm = std::sqrt(nx_on_yz * nx_on_yz + ny_on_yz * ny_on_yz + nz_on_yz * nz_on_yz);
                    // Dot product with +Z axis.
                    double dot = 0 * nx_on_yz + 0 * ny_on_yz + 1 * nz_on_yz;
                    // Compute the angle with the XY plane.
                    slopeAtlas.north(i, j) = 180.0 / M_PI * std::acos(dot / nz_on_yz_norm) * ((hit->ny >= 0) ? -1.0 : 1.0);
                }

                // Project the normal onto the XZ plane and normalize it.
                // Then compute its angle with the +X axis.
                {
                    // Set y to zero to project onto the XZ plane.
                    double nx_on_xz = hit->nx;
                    double ny_on_xz = 0.0;
                    double nz_on_xz = hit->nz;
                    double nz_on_xz_norm = std::sqrt(nx_on_xz * nx_on_xz + ny_on_xz * ny_on_xz + nz_on_xz * nz_on_xz);
                    // Dot product with +Z axis.
                    double dot = 0 * nx_on_xz + 0 * ny_on_xz + 1 * nz_on_xz;
                    // Compute the angle with the XY plane.
                    slopeAtlas.east(i, j) = 180.0 / M_PI * std::acos(dot / nz_on_xz_norm) * ((hit->nx >= 0) ? -1.0 : 1.0);
                }

                // Compute the directional slope in the northeast and southeast directions.
                {
                    double eSlope = slopeAtlas.east(i, j);
                    double nSlope = slopeAtlas.north(i, j);
                    slopeAtlas.northeast(i, j) =
                        90 - 180 / M_PI * std::acos(0.5 * (std::sin(M_PI / 180 * eSlope) + std::sin(M_PI / 180 * nSlope)));
                    slopeAtlas.southeast(i, j) =
                        90 - 180 / M_PI * std::acos(0.5 * (std::sin(M_PI / 180 * eSlope) - std::sin(M_PI / 180 * nSlope)));
                }

            } else {
                elevationMap(i, j) = minZ;
                priorityMap(i, j) = 0.0;

                slopeAtlas.north(i, j) = 0.0;
                slopeAtlas.northeast(i, j) = 0.0;
                slopeAtlas.east(i, j) = 0.0;
                slopeAtlas.southeast(i, j) = 0.0;
                slopeAtlas.absolute(i, j) = 0.0;
            }
        }
    }
    return std::make_tuple(elevationMap, priorityMap, slopeAtlas);
}
TerrainMapU8 buildReachabilityMap(const TerrainMapFloat &commsMap, const SlopeAtlas &slopeAtlas, const Path::State &landingSite) {

    TerrainMapU8 reachMap(commsMap.width(), commsMap.height(), commsMap.pitch);

    std::vector<Path::State> open;
    open.push_back(landingSite);

    int iterations = 0;
    while (!open.empty()) {
        if (iterations++ % (1 << 14) == 0) {
            fmt::print("[{}] Building reach map.\n", iterations);
        }

        // Pop currState from open.
        const auto currState = open.back();
        open.pop_back();

        // Otherwise, mark it as reached.
        setMapBit(reachMap, currState);

        // Add successors of the current state to the open stack.
        const auto &successors = getSuccessors(currState, slopeAtlas, commsMap);
        for (const auto &succ : successors) {
            if (!checkMapBit(reachMap, succ)) {
                open.push_back(succ);
            }
        }
    }
    return reachMap;
}
TerrainMapFloat computeDistanceTransform(const TerrainMapFloat &commsMap) {
    TerrainMapFloat distanceMap(commsMap.width(), commsMap.height(), commsMap.pitch);
    auto neighbors = [&commsMap](std::pair<int, int> point) {
        int i = point.first;
        int j = point.second;
        std::vector<std::pair<int, int>> succs;
        succs.reserve(8);
        if (i > 0) {
            if (j > 0)
                succs.emplace_back(std::make_pair(i - 1, j - 1));
            succs.emplace_back(std::make_pair(i - 1, j));
            if (j < commsMap.cols - 1)
                succs.emplace_back(std::make_pair(i - 1, j + 1));
        }
        if (j > 0)
            succs.emplace_back(std::make_pair(i, j - 1));
        if (j < commsMap.cols - 1)
            succs.emplace_back(std::make_pair(i, j + 1));
        if (i < commsMap.rows - 1) {
            if (j > 0)
                succs.emplace_back(std::make_pair(i + 1, j - 1));
            succs.emplace_back(std::make_pair(i + 1, j));
            if (j < commsMap.cols - 1)
                succs.emplace_back(std::make_pair(i + 1, j + 1));
        }
        return succs;
    };

    PriorityQueue<std::pair<int, int>, float> open;

    for (int i = 0; i < distanceMap.rows; i++) {
        for (int j = 0; j < distanceMap.cols; j++) {
            if (commsMap(i, j) == 0)
                distanceMap(i, j) = std::numeric_limits<float>::infinity();
            else
                distanceMap(i, j) = 0.0f;
            auto curr_neighbors = neighbors(std::make_pair(i, j));
            for (const auto &pair : curr_neighbors) {
                if (commsMap(pair.first, pair.second) == 0) {
                    // on contour
                    open.insert(std::make_pair(i, j), 0.0f);
                    break;
                }
            }
        }
    }
    while (!open.empty()) {
        std::pair<int, int> curr = open.top();
        open.pop();

        auto succs = neighbors(curr);
        for (const auto &neighbor : succs) {
            float neighbordt = distanceMap(neighbor.first, neighbor.second);
            float curdt = distanceMap(curr.first, curr.second);
            if (neighbordt > curdt + 1.0f) {
                distanceMap(neighbor.first, neighbor.second) = curdt + 1.0f;
                open.insert(neighbor, curdt + 1.0f);
            }
        }
    }
    return distanceMap;
}

TerrainMapFloat buildCommsMap(const TerrainMesh &tmesh, const TerrainMapFloat &elevationMap, const Path::State &landingSite) {

    TerrainMapFloat commsMap(elevationMap.width(), elevationMap.height(), elevationMap.pitch);

    const int li = landingSite.i;
    const int lj = landingSite.j;
    const double landerX = elevationMap.j2x(lj);
    const double landerY = elevationMap.i2y(li);
    const double landerZ = elevationMap(li, lj) + config.landerHeight;

    for (int ri = 0; ri < commsMap.rows; ++ri) {
        for (int rj = 0; rj < commsMap.cols; ++rj) {
            TerrainMesh::Ray ray;
            const double roverX = elevationMap.j2x(rj);
            const double roverY = elevationMap.i2y(ri);
            const double roverZ = elevationMap(ri, rj) + config.roverHeight;
            ray.oX = landerX;
            ray.oY = landerY;
            ray.oZ = landerZ;
            ray.dX = roverX - landerX;
            ray.dY = roverY - landerY;
            ray.dZ = roverZ - landerZ;
            const auto hit = tmesh.raytrace(ray);
            if (!hit) {
                commsMap(ri, rj) = 1;
            }
            double roverDist2 =
                (landerX - roverX) * (landerX - roverX) + (landerY - roverY) * (landerY - roverY) + (landerZ - roverZ) * (landerZ - roverZ);
            double hitDist2 =
                (ray.oX - hit->x) * (ray.oX - hit->x) + (ray.oY - hit->y) * (ray.oY - hit->y) + (ray.oZ - hit->z) * (ray.oZ - hit->z);
            if (hitDist2 >= 0.98 * roverDist2) {
                commsMap(ri, rj) = 1;
            }
        }
    }

    return commsMap;
}

TerrainMapFloat buildCoverageMap(const TerrainMesh &tmesh, const TerrainMapFloat &elevationMap, const std::vector<Vantage> &vantages) {

    TerrainMapFloat coverageMap(elevationMap.width(), elevationMap.height(), elevationMap.pitch);

    for (int vi = 0; vi < vantages.size(); ++vi) {
        auto &v = vantages[vi];
        fmt::print("[{}/{}] Generating coverage maps.\n", vi, vantages.size());
#pragma omp parallel for
        for (int ri = 0; ri < coverageMap.rows; ++ri) {
            for (int rj = 0; rj < coverageMap.cols; ++rj) {
                const double terrainX = elevationMap.j2x(rj);
                const double terrainY = elevationMap.i2y(ri);
                const double terrainZ = elevationMap(ri, rj);

                // Is this ray in front of the rover and within its field of view?
                // FIXME(Jordan): Do this calculation in 3d, not in 2d.
                {
                    const double roverAngle = directionToDegrees(v.dir) * M_PI / 180.0;
                    const double rx = std::cos(roverAngle);
                    const double ry = std::sin(roverAngle);

                    const double vtx = terrainX - v.x;
                    const double vty = terrainY - v.y;
                    const double vtnorm = std::sqrt(vtx * vtx + vty * vty);

                    const double viewAngle = 180.0 / M_PI * std::acos((rx * vtx + ry * vty) / vtnorm);
                    if (std::abs(viewAngle) > config.roverFOV / 2) {
                        continue;
                    }
                }

                TerrainMesh::Ray ray;
                ray.oX = v.x;
                ray.oY = v.y;
                ray.oZ = v.z + config.roverHeight;
                ray.dX = terrainX - ray.oX;
                ray.dY = terrainY - ray.oY;
                ray.dZ = terrainZ - ray.oZ;
                double rayNorm = std::sqrt(ray.dX * ray.dX + ray.dY * ray.dY + ray.dZ * ray.dZ);
                ray.dX /= rayNorm;
                ray.dY /= rayNorm;
                ray.dZ /= rayNorm;
                const auto hit = tmesh.raytrace(ray);
                if (hit) {
                    double hitAngle = 180 / M_PI * std::acos(-ray.dX * hit->nx - ray.dY * hit->ny - ray.dZ * hit->nz);
                    double hitDist = std::sqrt((ray.oX - hit->x) * (ray.oX - hit->x) + (ray.oY - hit->y) * (ray.oY - hit->y) +
                                               (ray.oZ - hit->z) * (ray.oZ - hit->z));
                    if (hitAngle < config.maxVisAngle && hitDist < config.maxVisRange && std::abs(rayNorm - hitDist) < rayNorm * 0.05) {
                        coverageMap(ri, rj)++;
                    }
                }
            }
        }
    }
    for (const auto &v : vantages) {
        const double roverAngle = directionToDegrees(v.dir) * M_PI / 180.0;
        const double rx = std::cos(roverAngle);
        const double ry = std::sin(roverAngle);
        drawCircle(coverageMap, v.x + 4 * rx, v.y + 4 * ry, (float)2 * vantages.size() + 10, 1.0);
        drawCircle(coverageMap, v.x, v.y, (float)2 * vantages.size(), 2.0);
    }
    return coverageMap;
}
