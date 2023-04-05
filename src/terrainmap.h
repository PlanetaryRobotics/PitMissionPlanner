#pragma once
#ifndef TERRAINMAP_H
#define TERRAINMAP_H

#include "tinyexr.h"
#include <cmath>
#include <fmt/format.h>
#include <fstream>
#include <vector>

template <typename T> class TerrainMap {
public:
  size_t rows;
  size_t cols;
  double pitch;

  TerrainMap(double mapWidth, double mapHeight, double mapPitch)
      : rows(std::ceil(mapHeight / mapPitch)), cols(std::ceil(mapWidth / mapPitch)), pitch(mapPitch),
        _data(rows * cols, 0) {}

  double width() const { return cols * pitch; }
  double height() const { return rows * pitch; }

  T &atIJ(int i, int j) { return operator()(i, j); }
  T &atXY(double x, double y) {
    int j = x2j(x);
    int i = y2i(y);
    return operator()(i, j);
  }
  const T &atIJ(int i, int j) const { return operator()(i, j); }
  const T &atXY(int x, int y) const {
    int j = x2j(x);
    int i = y2i(y);
    return operator()(i, j);
  }

  const T min() const {
    T minVal = std::numeric_limits<T>::max();
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        if (atIJ(i, j) < minVal) {
          minVal = atIJ(i, j);
        }
      }
    }
    return minVal;
  }
  const T max() const {
    T maxVal = -std::numeric_limits<T>::max();
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        if (atIJ(i, j) > maxVal) {
          maxVal = atIJ(i, j);
        }
      }
    }
    return maxVal;
  }

  double j2x(size_t j) const { return (j + 0.5) * pitch; }
  double i2y(size_t i) const { return (rows - i - 0.5) * pitch; }

  size_t x2j(double x) const { return x / pitch - 0.5; }
  size_t y2i(double y) const { return rows - 0.5 - y / pitch; }

  T &operator()(size_t i, size_t j) { return _data[i * cols + j]; }
  const T &operator()(size_t i, size_t j) const { return _data[i * cols + j]; }

  const std::vector<T> &data() const { return _data; }

private:
  std::vector<T> _data;
};

using TerrainMapFloat = TerrainMap<float>;
using TerrainMapU8 = TerrainMap<uint8_t>;

template <typename T> void drawCircle(TerrainMap<T> &map, double x, double y, double val, double rad) {
  int cj = map.x2j(x);
  int ci = map.y2i(y);
  int R = std::ceil(rad / map.pitch);

  for (int i = ci - R; i <= ci + R; ++i) {
    for (int j = cj - R; j <= cj + R; ++j) {
      if (i < 0 || i > map.rows) {
        continue;
      }
      if (j < 0 || j > map.cols) {
        continue;
      }
      if ((i - ci) * (i - ci) + (j - cj) * (j - cj) < R * R) {
        map(i, j) = static_cast<T>(val);
      }
    }
  }
}

// x and y are center of rectangle
template <typename T>
void drawTriangle(TerrainMap<T> &map, double X, double Y, int d, double val, double width = 2, double height = 4) {
  int cj = map.x2j(X);
  int ci = map.y2i(Y);

  int W = std::ceil(width / map.pitch);
  int H = std::ceil(height / map.pitch);

  double top[2] = {ci - H / 2.0, (double)cj};
  double bl[2] = {ci + H / 2.0, cj - W / 2.0};
  double br[2] = {ci + H / 2.0, cj + W / 2.0};

  auto rotCW = [ci, cj](double i, double j, double deg) {
    double di = i - ci;
    double dj = j - cj;
    double c = std::cos(deg * M_PI / 180.0);
    double s = std::sin(deg * M_PI / 180.0);
    double ri = c * di + s * dj + ci;
    double rj = -s * di + c * dj + cj;
    return std::make_pair(ri, rj);
  };
  auto det2d = [](double ai, double aj, double bi, double bj) { return ai * bj - aj * bi; };

  // Rotate the triangle vertices by the appropriate angle.
  std::tie(top[0], top[1]) = rotCW(top[0], top[1], 45 * d);
  std::tie(bl[0], bl[1]) = rotCW(bl[0], bl[1], 45 * d);
  std::tie(br[0], br[1]) = rotCW(br[0], br[1], 45 * d);

  // Compute the bounding box of the rotated triangle.
  int minI = std::min(std::min(top[0], bl[0]), br[0]);
  int maxI = std::max(std::max(top[0], bl[0]), br[0]);
  int minJ = std::min(std::min(top[1], bl[1]), br[1]);
  int maxJ = std::max(std::max(top[1], bl[1]), br[1]);

  // Loop over the bounding box of the rotated triangle.
  for (int i = minI; i <= maxI; ++i) {
    for (int j = minJ; j <= maxJ; ++j) {
      // If this point is outside the map, skip it.
      if (i < 0 || i >= map.rows) {
        continue;
      }
      if (j < 0 || j >= map.cols) {
        continue;
      }

      // Compute the barycentric coordinates of this point.
      double area = W * H;
      double beta = det2d(j - top[1], i - top[0], top[1] - br[1], top[0] - br[0]) / area;
      double gamma = det2d(j - top[1], i - top[0], bl[1] - top[1], bl[0] - top[0]) / area;
      double alpha = 1.0 - beta - gamma;

      const double eps = 1e-4;
      if (0.0 - eps <= alpha && alpha <= 1.0 + eps && 0.0 - eps <= beta && beta <= 1.0 + eps && 0.0 - eps <= gamma &&
          gamma <= 1.0 + eps) {
        map(i, j) = static_cast<T>(val);
      }
    }
  }
}

template <typename T> void saveEXR(const TerrainMap<T> &map, const std::string &filename) {
  const char *err;
  int success = SaveEXR(map.data().data(), map.cols, map.rows, 1, 0, filename.c_str(), &err);

  if (success != TINYEXR_SUCCESS) {
    fmt::print("TINYEXR ERROR: {}\n", err);
  }
}

template <typename T> void savePFM(const TerrainMap<T> &map, const std::string &filename) {
  std::ofstream os(filename, std::ios::out);
  os << fmt::format("Pf\n");
  os << fmt::format("{} {}\n", map.cols, map.rows);
  os << fmt::format("-1\n");
  for (int i = 0; i < map.rows; ++i) {
    for (int j = 0; j < map.cols; ++j) {
      float f = static_cast<float>(map.data()[i * map.cols + j]);
      os.write(reinterpret_cast<const char *>(&f), sizeof(float));
    }
  }
  os.close();
}

#endif // TERRAINMAP_H
