#include "terrainmesh.h"
#include "happly.h"
#include <embree3/rtcore.h>
#include <fmt/format.h>
#include <filesystem>
#include <stdexcept>
#include <optional>

TerrainMesh::TerrainMesh(const std::string& filename) : filename(filename) {
    namespace fs = std::filesystem;
    const std::string extension = fs::path(filename).extension();
    if( extension != ".ply" ) {
        throw std::runtime_error(fmt::format("Failed to load .ply mesh with unsupported extension \"{}\"", extension));
    }

    happly::PLYData plyIn(filename);

    const auto verts = plyIn.getVertexPositions();
    for(const auto& v : verts) {
        vertexPositions.push_back( static_cast<float>(v[0]) );
        vertexPositions.push_back( static_cast<float>(v[1]) );
        vertexPositions.push_back( static_cast<float>(v[2]) );
    }
    const auto faces = plyIn.getFaceIndices<size_t>();
    for(const auto& f : faces) {
        faceIndices.push_back(f[0]);
        faceIndices.push_back(f[1]);
        faceIndices.push_back(f[2]);
    }
    const auto nx = plyIn.getElement("vertex").getProperty<float>("nx");
    const auto ny = plyIn.getElement("vertex").getProperty<float>("ny");
    const auto nz = plyIn.getElement("vertex").getProperty<float>("nz");
    for(int i=0; i<nx.size(); ++i) {
        vertexNormals.push_back(nx[i]);
        vertexNormals.push_back(ny[i]);
        vertexNormals.push_back(nz[i]);
    }
    const auto vertexColors = plyIn.getVertexColors();
    for(const auto& c: vertexColors) {
        float mean = (c[0]+c[1]+c[2]) / 3.0f;
        vertexPriorities.push_back(std::pow(mean / 255.0f, 8.0));
    }

    const double minX = this->minX();
    const double minY = this->minY();
    if( std::abs(minX) > 1e-3 ) {
        throw std::runtime_error(fmt::format("Terrain mesh minimum X coordinate ({}) is not zero.\n", minX));
    }
    if( std::abs(minY) > 1e-3 ) {
        throw std::runtime_error(fmt::format("Terrain mesh minimum Y coordinate ({}) is not zero.\n", minY));
    }

    device = rtcNewDevice(NULL);
    scene = rtcNewScene(device);
    RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    float* vb = (float*)rtcSetNewGeometryBuffer(geometry, RTC_BUFFER_TYPE_VERTEX, 0,
                                                RTC_FORMAT_FLOAT3, 3*sizeof(float),
                                                vertexPositions.size()/3);
    std::memcpy(vb, vertexPositions.data(), vertexPositions.size() * sizeof(float));

    unsigned* ib = (unsigned *)rtcSetNewGeometryBuffer(geometry, RTC_BUFFER_TYPE_INDEX, 0,
                                                       RTC_FORMAT_UINT3, 3*sizeof(unsigned),
                                                       faceIndices.size()/3);
    std::memcpy(ib, faceIndices.data(), faceIndices.size() * sizeof(unsigned));

    rtcCommitGeometry(geometry);
    rtcAttachGeometry(scene, geometry);
    rtcReleaseGeometry(geometry);
    rtcCommitScene(scene);
}

TerrainMesh::~TerrainMesh() {
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);
}

double TerrainMesh::minX() const {
    double min = std::numeric_limits<double>::infinity();
    for(int i=0; i<vertexPositions.size(); i+=3) { min = std::min<double>(min, vertexPositions[i]); }
    return min;
}
double TerrainMesh::minY() const {
    double min = std::numeric_limits<double>::infinity();
    for(int i=1; i<vertexPositions.size(); i+=3) { min = std::min<double>(min, vertexPositions[i]); }
    return min;
}
double TerrainMesh::minZ() const {
    double min = std::numeric_limits<double>::infinity();
    for(int i=2; i<vertexPositions.size(); i+=3) { min = std::min<double>(min, vertexPositions[i]); }
    return min;
}
double TerrainMesh::maxX() const {
    double max = -std::numeric_limits<double>::infinity();
    for(int i=0; i<vertexPositions.size(); i+=3) { max = std::max<double>(max, vertexPositions[i]); }
    return max;
}
double TerrainMesh::maxY() const {
    double max = -std::numeric_limits<double>::infinity();
    for(int i=1; i<vertexPositions.size(); i+=3) { max = std::max<double>(max, vertexPositions[i]); }
    return max;
}
double TerrainMesh::maxZ() const {
    double max = -std::numeric_limits<double>::infinity();
    for(int i=2; i<vertexPositions.size(); i+=3) { max = std::max<double>(max, vertexPositions[i]); }
    return max;
}

std::optional<TerrainMesh::Hit>
TerrainMesh::raytrace(const TerrainMesh::Ray& ray) const {
    RTCRayHit rayhit;
    rayhit.ray.org_x = ray.oX;
    rayhit.ray.org_y = ray.oY;
    rayhit.ray.org_z = ray.oZ;
    rayhit.ray.dir_x = ray.dX;
    rayhit.ray.dir_y = ray.dY;
    rayhit.ray.dir_z = ray.dZ;
    rayhit.ray.tnear = 0.0f;
    rayhit.ray.tfar  = std::numeric_limits<float>::infinity();
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    rtcIntersect1(scene, &context, &rayhit);

    if( rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID ) {
        return {};
    }

    const float u = 1-rayhit.hit.u-rayhit.hit.v;
    const float v = rayhit.hit.u;
    const float w = rayhit.hit.v;

    const size_t tid = rayhit.hit.primID;
    const size_t Aid = faceIndices[3*tid+0];
    const size_t Bid = faceIndices[3*tid+1];
    const size_t Cid = faceIndices[3*tid+2];

    TerrainMesh::Hit hit;
    hit.x  = u*vertexPositions[3*Aid+0] + v*vertexPositions[3*Bid+0] + w*vertexPositions[3*Cid+0];
    hit.y  = u*vertexPositions[3*Aid+1] + v*vertexPositions[3*Bid+1] + w*vertexPositions[3*Cid+1];
    hit.z  = u*vertexPositions[3*Aid+2] + v*vertexPositions[3*Bid+2] + w*vertexPositions[3*Cid+2];
    hit.nx = u*vertexNormals[3*Aid+0] + v*vertexNormals[3*Bid+0] + w*vertexNormals[3*Cid+0];
    hit.ny = u*vertexNormals[3*Aid+1] + v*vertexNormals[3*Bid+1] + w*vertexNormals[3*Cid+1];
    hit.nz = u*vertexNormals[3*Aid+2] + v*vertexNormals[3*Bid+2] + w*vertexNormals[3*Cid+2];
    hit.priority = u*vertexPriorities[Aid] + v*vertexPriorities[Bid] + w*vertexPriorities[Cid];
    hit.faceIndex = tid;

    return hit;
}
