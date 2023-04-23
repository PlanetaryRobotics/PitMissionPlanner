#pragma once
#ifndef TERRAINMESH_H
#define TERRAINMESH_H

#include <array>
#include <embree3/rtcore.h>
#include <limits>
#include <optional>
#include <string>
#include <vector>

class TerrainMesh {
  public:
    TerrainMesh(const std::string &meshfile);
    ~TerrainMesh();

    double minX() const;
    double minY() const;
    double minZ() const;

    double maxX() const;
    double maxY() const;
    double maxZ() const;

    struct Ray {
        double oX, oY, oZ;
        double dX, dY, dZ;
    };

    struct Hit {
        double x, y, z;
        double priority;
        double nx, ny, nz;
        size_t faceIndex;
    };

    std::optional<Hit> raytrace(const Ray &ray) const;

  private:
    const std::string &filename;
    std::vector<float> vertexPositions;
    std::vector<float> vertexNormals;
    std::vector<float> vertexPriorities;
    std::vector<unsigned> faceIndices;

    RTCDevice device;
    RTCScene scene;
};

#endif // TERRAINMESH_H
