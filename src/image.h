#pragma once
#ifndef IMAGE_H
#define IMAGE_H

#include "terrainmap.h"
#include "tinycolormap.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <fmt/format.h>

class ImageRGB {
    public:
        ImageRGB(const TerrainMap<float>& map, const tinycolormap::ColormapType& cmapType, int upsample = 1, double exposure=0.0, double offset=0.0) :
            _rows(upsample*map.rows), _cols(upsample*map.cols), _exposure(exposure), _offset(offset) {
            _data.resize(_rows * _cols * _channels);

            float maxVal = map.max();
            float minVal = map.min();

            for(int i=0; i<_rows; ++i) {
                for(int j=0; j<_cols; ++j) {
                    float v = map(i/upsample, j/upsample);

                    v = std::pow<float>(2, _exposure) * v + _offset;
                    v = (v-minVal) / (maxVal-minVal);

                    const auto c = tinycolormap::GetColor(v, cmapType);

                    int pxIdx = i*_cols+j;
                    _data[_channels*pxIdx + 0] = c.ri();
                    _data[_channels*pxIdx + 1] = c.gi();
                    _data[_channels*pxIdx + 2] = c.bi();
                }
            }
        }

        int rows() const { return _rows; }
        int cols() const { return _cols; }
        int width() const { return _cols; }
        int height() const { return _rows; }
        int channels() const { return _channels; }
        const uint8_t* data() const { return _data.data(); }

        uint8_t& operator()(int i, int j, int c) { return _data[(i*_cols + j)*_channels + c]; }
        const uint8_t& operator()(int i, int j, int c) const { return _data[(i*_cols + j)*_channels + c]; }

    private:
        const int _rows = 0;
        const int _cols = 0;
        const int _channels = 3;
        const double _exposure;
        const double _offset;
        std::vector<uint8_t> _data;
};

void drawCircle(ImageRGB& image, int ci, int cj, double r, double R, const tinycolormap::Color& color) {
    for(int i=ci-R; i<=ci+R; ++i) {
        for(int j=cj-R; j<=cj+R; ++j) {
            if( i < 0 || i > image.rows() ) { continue; }
            if( j < 0 || j > image.cols() ) { continue; }
            double d2 = (i-ci)*(i-ci) + (j-cj)*(j-cj);
            if( r*r <= d2 && d2 <= R*R ) {
                image(i,j,0) = color.ri();
                image(i,j,1) = color.gi();
                image(i,j,2) = color.bi();
            }
        }
    }
}

void drawCircle(ImageRGB& image, int ci, int cj, double R, const tinycolormap::Color& color) {
    drawCircle(image, ci, cj, R, 0, color);
}

void drawTriangle(ImageRGB& image, int ci, int cj, double deg, const tinycolormap::Color& color, double W = 2, double H = 4) {
    double top[2] = {ci-H/2.0, (double)cj};
    double  bl[2] = {ci+H/2.0, cj-W/2.0};
    double  br[2] = {ci+H/2.0, cj+W/2.0};

    auto rotCW = [ci, cj](double i, double j, double deg) {
        double di = i-ci;
        double dj = j-cj;
        double c = std::cos(deg * M_PI/180.0);
        double s = std::sin(deg * M_PI/180.0);
        double ri =  c*di + s*dj + ci;
        double rj = -s*di + c*dj + cj;
        return std::make_pair(ri, rj);
    };
    auto det2d = [](double ai, double aj, double bi, double bj) {
        return ai*bj - aj*bi;
    };

    // Rotate the triangle vertices by the appropriate angle.
    std::tie(top[0], top[1]) = rotCW(top[0], top[1], deg);
    std::tie( bl[0],  bl[1]) = rotCW( bl[0],  bl[1], deg);
    std::tie( br[0],  br[1]) = rotCW( br[0],  br[1], deg);

    // Compute the bounding box of the rotated triangle.
    int minI = std::min(std::min(top[0], bl[0]), br[0]);
    int maxI = std::max(std::max(top[0], bl[0]), br[0]);
    int minJ = std::min(std::min(top[1], bl[1]), br[1]);
    int maxJ = std::max(std::max(top[1], bl[1]), br[1]);

    // Loop over the bounding box of the rotated triangle.
    for(int i=minI; i<=maxI; ++i) {
        for(int j=minJ; j<=maxJ; ++j) {
            // If this point is outside the image, skip it.
            if( i < 0 || i >= image.rows() ) { continue; }
            if( j < 0 || j >= image.cols() ) { continue; }

            // Compute the barycentric coordinates of this point.
            double area  = W*H;
            double beta  = det2d(j-top[1], i-top[0],  top[1]-br[1], top[0]-br[0]) / area;
            double gamma = det2d(j-top[1], i-top[0],  bl[1]-top[1], bl[0]-top[0]) / area;
            double alpha = 1.0 - beta - gamma;

            const double eps = 1e-4;
            if( 0.0 - eps <= alpha && alpha <= 1.0 + eps &&
                0.0 - eps <= beta  && beta  <= 1.0 + eps &&
                0.0 - eps <= gamma && gamma <= 1.0 + eps ) {
                image(i,j, 0) = color.ri();
                image(i,j, 1) = color.gi();
                image(i,j, 2) = color.bi();
            }
        }
    }
}

void drawLine(ImageRGB& image, double ai, double aj, double bi, double bj, const tinycolormap::Color& color, double W=1.0) {
    // Lazy man's line drawing algorithm...
    double length = std::sqrt((bi-ai)*(bi-ai) + (bj-aj)*(bj-aj));

    double linevec[2] = {(bi-ai)/length, (bj-aj)/length};
    double norm[2] = {-linevec[1], linevec[0]};

    double W_2 = W/2.0;
    double minI = std::min(ai, bi) - W_2;
    double maxI = std::max(ai, bi) + W_2;
    double minJ = std::min(aj, bj) - W_2;
    double maxJ = std::max(aj, bj) + W_2;

    auto dist_point_to_line = [](double  i, double  j,
                                 double ai, double aj,
                                 double bi, double bj) {
        double l2 = (bi-ai)*(bi-ai) + (bj-aj)*(bj-aj);
        if( l2 == 0.0 ) { return std::sqrt((i-ai)*(i-ai)+(j-bj)*(j-bj)); }

        double dot = (i-ai)*(bi-ai) + (j-aj)*(bj-aj);
        double t = std::clamp(dot/l2, 0.0, 1.0);
        double proj_i = ai + t*(bi-ai);
        double proj_j = aj + t*(bj-aj);
        return std::sqrt((i-proj_i)*(i-proj_i)+(j-proj_j)*(j-proj_j));
    };

    for(int i=minI; i<=maxI; ++i) {
        for(int j=minJ; j<=maxJ; ++j) {
            // If this point is outside the image, skip it.
            if( i < 0 || i >= image.rows() ) { continue; }
            if( j < 0 || j >= image.cols() ) { continue; }

            double d2l = dist_point_to_line(i,j, ai,aj, bi,bj); 
            if( d2l <= W_2 ) {
                image(i,j, 0) = color.ri();
                image(i,j, 1) = color.gi();
                image(i,j, 2) = color.bi();
            }
        }
    }
}

void savePNG(const ImageRGB& image, const std::string& path) {
    int stride = image.width() * image.channels();
    int err = stbi_write_png(path.c_str(), image.width(), image.height(), image.channels(), image.data(), stride);
}

#endif // IMAGE_H
