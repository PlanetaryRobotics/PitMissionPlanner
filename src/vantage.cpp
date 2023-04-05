#include "vantage.h"
#include <Eigen/Dense>

#include <unordered_set>
std::pair<std::vector<Probe>, TerrainMapFloat> generateVisibilityProbes(const TerrainMapFloat &priorityMap,
                                                                        const TerrainMapFloat &elevationMap) {
  std::vector<Probe> probes;
  TerrainMapFloat probeMap(priorityMap.width(), priorityMap.height(), priorityMap.pitch);
  while (probes.size() < config.numProbes) {
    int i = priorityMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    int j = priorityMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    if (priorityMap(i, j) > 0) {
      Probe p;
      p.x = priorityMap.j2x(j);
      p.y = priorityMap.i2y(i);
      p.z = elevationMap(i, j);
      p.priority = priorityMap(i, j);
      probeMap(i, j) = p.priority;
      probes.push_back(p);
    }
  }
  return std::make_pair(probes, probeMap);
}
std::pair<std::vector<Vantage>, TerrainMapFloat> generateVantageCandidates(const TerrainMesh &tmesh,
                                                                           const TerrainMapU8 &reachMap,
                                                                           const TerrainMapFloat &elevationMap,
                                                                           const std::vector<Probe> probes) {
  // Generate thousands of random, reachable points at rover height above the terrain.
  // NOTE(Jordan): I'm currently duplicating each candidate vantage eight times,
  // once per rover facing direction. This is fairly wasteful. We should probably evaluate
  // all directions at each location and pick the one direction that has the best coverage.
  std::vector<Vantage> candidates;
  while (candidates.size() < config.numCandidates) {
    Path::State state;
    state.i = reachMap.rows * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    state.j = reachMap.cols * static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

    for (int d = 0; d < 8; ++d) {
      state.d = static_cast<Direction>(d);

      if (checkMapBit(reachMap, state)) {
        Vantage v;
        v.x = reachMap.j2x(state.j);
        v.y = reachMap.i2y(state.i);
        v.z = elevationMap(state.i, state.j) + config.roverHeight;
        v.dir = state.d;
        candidates.push_back(v);
        fmt::print("{},{},{} -> {:<#010b}\n", state.i, state.j, directionToString(state.d), reachMap(state.i, state.j));
      }
    }
  }

  TerrainMapFloat candidateMap(reachMap.width(), reachMap.height(), reachMap.pitch);
  for (int i = 0; i < candidateMap.rows; ++i) {
    for (int j = 0; j < candidateMap.cols; ++j) {
      candidateMap(i, j) = (reachMap(i, j) == 0x00) ? 0.0f : 1.0f;
    }
  }

  int progress = 0;

// For every candidate, trace a ray to all view coverage probes.
#pragma omp parallel for shared(progress)
  for (int ci = 0; ci < candidates.size(); ++ci) {
#pragma omp critical
    { fmt::print("[{}/{}] Constructing visibility map.\n", progress++, candidates.size()); }
    auto &candidate = candidates[ci];
    candidate.coverage.resize(probes.size());

    for (int pi = 0; pi < probes.size(); ++pi) {
      auto &p = probes[pi];

      // Is this ray in front of the rover and within its field of view?
      // FIXME(Jordan): Do this calculation in 3d, not in 2d.
      {
        const double roverAngle = directionToDegrees(candidate.dir) * M_PI / 180.0;
        const double rx = std::cos(roverAngle);
        const double ry = std::sin(roverAngle);

        const double cpx = p.x - candidate.x;
        const double cpy = p.y - candidate.y;
        const double cpnorm = std::sqrt(cpx * cpx + cpy * cpy);

        const double viewAngle = 180.0 / M_PI * std::acos((rx * cpx + ry * cpy) / cpnorm);
        if (std::abs(viewAngle) > config.roverFOV / 2) {
          continue;
        }
      }

      TerrainMesh::Ray ray;
      ray.oX = candidate.x;
      ray.oY = candidate.y;
      ray.oZ = candidate.z + config.roverHeight;
      ray.dX = p.x - candidate.x;
      ray.dY = p.y - candidate.y;
      ray.dZ = p.z - candidate.z;
      double rayNorm = std::sqrt(ray.dX * ray.dX + ray.dY * ray.dY + ray.dZ * ray.dZ);
      ray.dX /= rayNorm;
      ray.dY /= rayNorm;
      ray.dZ /= rayNorm;
      const auto hit = tmesh.raytrace(ray);
      if (hit) {
        double hitAngle = 180 / M_PI * std::acos(-ray.dX * hit->nx - ray.dY * hit->ny - ray.dZ * hit->nz);
        double hitDist = std::sqrt((ray.oX - hit->x) * (ray.oX - hit->x) + (ray.oY - hit->y) * (ray.oY - hit->y) +
                                   (ray.oZ - hit->z) * (ray.oZ - hit->z));
        if (hitAngle < config.maxVisAngle && hitDist < config.maxVisRange &&
            std::abs(rayNorm - hitDist) < 0.05 * rayNorm) {
          candidate.coverage[pi] = true;
          candidate.totalCoverage += p.priority;
        }
      }
    }
    int i = candidateMap.y2i(candidate.y);
    int j = candidateMap.x2j(candidate.x);
    candidateMap(i, j) = 10 + std::max<float>(candidate.totalCoverage, candidateMap(i, j));
  }
  return std::make_pair(candidates, candidateMap);
}

std::vector<Vantage> selectVantages(const std::vector<Vantage> &candidates, const std::vector<Probe> &probes) {

  std::vector<Vantage> vantages;
  std::unordered_set<int> taken;
  std::vector<unsigned char> visCounters(probes.size(), 0);

  // Make k selections. Choose the candidate that produces the greatest *new* coverage.
  for (int k = 0; k < config.numVantages; ++k) {
    fmt::print("[{}/{}] Selecting vantages.\n", k, config.numVantages);

    // Assign a score to every candidate.
    std::vector<float> scores(candidates.size(), 0.0f);
#pragma omp parallel for
    for (int ci = 0; ci < candidates.size(); ++ci) {
      if (taken.contains(ci)) {
        continue;
      }
      const auto &c = candidates[ci];

      // If this candidate is too close to one
      // that has already been selected, skip it.
      bool tooClose = false;
      for (const auto &ti : taken) {
        const auto &t = candidates[ti];
        double d2 = (t.x - c.x) * (t.x - c.x) + (t.y - c.y) * (t.y - c.y) + (t.z - c.z) * (t.z - c.z);
        if (d2 < config.minVantageSeparation * config.minVantageSeparation) {
          tooClose = true;
          break;
        }
      }
      if (tooClose) {
        continue;
      }

      // Otherwise, calculate "how good" this candidate is
      // based on the visibility probes it can see.
      // There are three factors to consider here:
      //   1. Does this candidate "see" new probes that haven't been seen yet? (countWeight)
      //   2. Are the probes this candidate can see in the center of the rover's FOV? (angleWeight)
      //   3. The science priority assigned to each probe by the user. (probes[pi].priority)
      for (int pi = 0; pi < probes.size(); ++pi) {
        // If this probe is visible, add to the score for this candidate.
        if (c.coverage[pi]) {
          // How many times has this probe been viewed?
          int visCount = std::clamp<int>(visCounters[pi], 0, 10);
          double countWeight = config.visMultipliers[visCount];

          // Compute the angle between (1) the camera's normal vector
          // and (2) the viewing vector from the candidate to the probe.
          double visAngle = 0.0;
          {
            // (1) A vector pointing down the boresight of the rover's camera.
            Eigen::Vector3f cam;
            double camAngle = M_PI / 180.0 * 45 * (int)c.dir;
            cam << std::sin(camAngle), std::cos(camAngle), 0.0;

            // (2) The normalized vector from the probe to the candidate.
            Eigen::Vector3f vis;
            vis << probes[pi].x - c.x, probes[pi].y - c.y, probes[pi].z - c.z;
            vis /= vis.norm();

            // How far (in deg) is the view vector from the center of the rover's FOV?
            visAngle = 180.0 / M_PI * std::acos(vis.dot(cam));
          }
          // Weight views in the center of the FOV more than views on the edges of the FOV.
          // Use a scaled gaussian w(x) = e^(-(x/90)^2). This is an arbitrary choice.
          double angleWeight = std::exp(-(visAngle / 90.0) * (visAngle / 90.0));

          // Add to the score for this vantage candidate.
          scores[ci] += countWeight * angleWeight * probes[pi].priority;
        }
      }
    }

    // Select the candidate with the highest score.
    int best_ci = 0;
    for (int ci = 0; ci < candidates.size(); ++ci) {
      if (scores[ci] > scores[best_ci]) {
        best_ci = ci;
      }
    }

    // If the best score is zero, give up because there are no more candidates worth considering.
    if (scores[best_ci] <= 0) {
      fmt::print("Warning: Only {} useful vantages were found of the {} that were requested.\n", vantages.size(),
                 config.numVantages);
      break;
    }
    taken.insert(best_ci);

    fmt::print("Score: {}\n", scores[best_ci]);

    // Add new visibility to the visCounters.
    for (int pi = 0; pi < probes.size(); ++pi) {
      if (candidates[best_ci].coverage[pi]) {
        visCounters[pi]++;
      }
    }
    vantages.push_back(candidates[best_ci]);
  }
  return vantages;
}
std::vector<Vantage> sortCCW(const std::vector<Vantage> vantages, double siteX, double siteY) {
  std::vector<Vantage> sorted = vantages;
  auto angle = [siteX, siteY](const Vantage &v0, const Vantage &v1) {
    double aX = v0.x - siteX;
    double aY = v0.y - siteY;
    double bX = v1.x - siteX;
    double bY = v1.y - siteY;
    return std::atan2(aY, aX) < std::atan2(bY, bX);
  };
  std::sort(sorted.begin(), sorted.end(), angle);
  return sorted;
}
