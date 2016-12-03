//
//  photoSet.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/21/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "photoSet.hpp"
#include <fstream>

using std::ifstream;

PhotoSet::PhotoSet() {
}

PhotoSet::~PhotoSet() {
}

void PhotoSet::init(const vector<int>& images, const string prefix, const int nimages, const int nillums, const int maxLevel, const int size, const int alloc) {
    m_images = images;
    m_nimages = nimages;
    m_nillums = nillums;
    m_prefix = prefix;
    
    for (int i = 0; i < m_nimages; ++i) {
        m_dict[images[i]] = i;
    }
    
    m_photos.resize(m_nimages);
    cerr << "Reading images: " << std::flush;
    for (int i = 0; i < m_nimages; ++i) {
        char img0[1024], img1[1024], img2[1024], img3[1024];
        sprintf(img0, "%simage/%04d%04d.jpg", prefix.c_str(), i, 0);
        sprintf(img1, "%simage/%04d%04d.ppm", prefix.c_str(), i, 0);
        sprintf(img2, "%simage/%04d%04d.png", prefix.c_str(), i, 0);
        sprintf(img3, "%simage/%04d%04d.tiff", prefix.c_str(), i, 0);
        if (ifstream(img0) || ifstream(img1)
#if defined(PMMVPS_HAVE_PNG)
            || ifstream(img2)
#endif
#if defined(PMMVPS_HAVE_TIFF)
            || ifstream(img3)
#endif
            ) {
            char iname[1024], cname[1024], mname[1024];
            sprintf(iname, "%simage/%04d", prefix.c_str(), i);
            sprintf(mname, "%smask/%08d", prefix.c_str(), i);
            sprintf(cname, "%stxt/%08d.txt", prefix.c_str(), i);
            m_photos[i].init(cname, iname, mname, nillums, maxLevel);
            if (alloc) {
                m_photos[i].alloc();
            }
            else {
                m_photos[i].alloc(1);
            }
        }
    }
    const int margin = size / 2;
    m_size = 2 * margin + 1;
}

void PhotoSet::free() {
    for (int i = 0; i < static_cast<int>(m_photos.size()); ++i) {
        // free the pyramid of images
        m_photos[i].free();
    }
}

void PhotoSet::free(const int level) {
    for (int i = 0; i < static_cast<int>(m_photos.size()); ++i) {
        // free one level of the pyramid of images
        m_photos[i].free(level);
    }
}

int PhotoSet::checkAngles(const Vector4f& coord, const vector<int>& indexes, const float minAngle, const float maxAngle, const int num) const {
    int count = 0;
    
    vector<Vector4f, Eigen::aligned_allocator<Vector4f> > rays;
    rays.resize(indexes.size());
    
    for (int i = 0; i < static_cast<int>(indexes.size()); ++i) {
        const int index = indexes[i];
        rays[i] = m_photos[index].m_center - coord;
        rays[i] /= rays[i].norm();
    }
    
    for (int i = 0; i < static_cast<int>(indexes.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(indexes.size()); ++j) {
            const float dot = std::max(-1.0f, std::min(1.0f, rays[i].dot(rays[j])));
            const float angle = acos(dot);
            if (minAngle < angle && angle < maxAngle) {
                ++count;
            }
        }
    }
    if (count < 1) {
        return -1;
    }
    else
        return 0;
}

void PhotoSet::setDistances() {
    m_distances.resize(m_nimages);
    float avedis = 0.0f;
    int denom = 0;
    
    // optical center distance
    for (int i = 0; i < m_nimages; ++i) {
        m_distances[i].resize(m_nimages);
        for (int j = 0; j < m_nimages; ++j) {
            if (i == j) {
                m_distances[i][j] = 0.0f;
            }
            else {
                m_distances[i][j] = (m_photos[i].m_center - m_photos[j].m_center).norm();
                avedis += m_distances[i][j];
                ++denom;
            }
        }
    }
    if (denom == 0) {
        return;
    }
    
    avedis /= denom;
    if (avedis == 0.0f) {
        cerr << "all the optical centers are identical." << endl;
        exit(1);
    }
    
    // angular distance
    for (int i = 0; i < m_nimages; ++i) {
        Vector4f ray0 = m_photos[i].m_oaxis;
        ray0(3) = 0.0f;
        for (int j = 0; j < m_nimages; ++j) {
            Vector4f ray1 = m_photos[j].m_oaxis;
            ray1(3) = 0.0f;
            
            m_distances[i][j] /= avedis;
            const float margin = cosf(10.0f * M_PI / 180.0f);
            const float dis = std::max(0.0f, 1.0f - ray0.dot(ray1) - margin); // angle below a specific value won't contribute to m_distances[i][j]
            m_distances[i][j] += dis;
        }
    }
}

/*
void PhotoSet::getTex(const int index, const int level, const Vector2f& icoord, const Vector2f& xaxis, const Vector2f& yaxis, vector<Vector3f>& tex, const int isNorm) const {
    m_photos[index].getTex(illum, level, icoord, xaxis, yaxis, m_size, tex, isNorm);
}

void PhotoSet::getTex(const int index, const int illum, const int level, const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, vector<Vector3f>& tex, float& weight, const int isNorm) const {
    m_photos[index].getTex(illum, level, coord, pxaxis, pyaxis, pzaxis, m_size, tex, weight, isNorm);
}
 */

Vector3f PhotoSet::getColor(const int index, const Vector4f& coord, const int level) const {
    return m_photos[index].getColor(coord, level);
}

Vector3f PhotoSet::getColor(const int index, const float fx, const float fy, const int level) const {
    return m_photos[index].getColor(fx, fy, level);
}

Vector3f PhotoSet::getColor(const int index, const int ix, const int iy, const int level) const {
    return m_photos[index].getColor(ix, iy, level);
}

Vector3f PhotoSet::getColor(const int index, const Vector4f& coord, const int level, const int illum) {
    return m_photos[index].getColor(coord, level, illum);
}

Vector3f PhotoSet::getColor(const int index, const float fx, const float fy, const int level, const int illum) {
    return m_photos[index].getColor(fx, fy, level, illum);
}

Vector3f PhotoSet::getColor(const int index, const int ix, const int iy, const int level, const int illum) {
    return m_photos[index].getColor(ix, iy, level, illum);
}

vector<Vector3f> PhotoSet::getColorSeq(const int index, const Vector4f& coord, const int level) {
    int sz = static_cast<int>(m_photos[index].m_imageSets.size());
    vector<Vector3f> colors(sz);
    for (int illum = 0; illum < sz; ++illum) {
        colors[illum] = m_photos[index].getColor(coord, illum, level);
    }
    return colors;
}

vector<Vector3f> PhotoSet::getColorSeq(const int index, const float fx, const float fy, const int level) {
    int nillums = static_cast<int>(m_photos[index].m_imageSets.size());
    vector<Vector3f> colors(nillums);
    for (int illum = 0; illum < nillums; ++illum) {
        colors[illum] = m_photos[index].getColor(fx, fy, level, illum);
    }
    return colors;
}

vector<Vector3f> PhotoSet::getColorSeq(const int index, const int ix, const int iy, const int level) {
    int sz = static_cast<int>(m_photos[index].m_imageSets.size());
    vector<Vector3f> colors(sz);
    for (int illum = 0; illum < sz; ++illum) {
        colors[illum] = m_photos[index].getColor(ix, iy, illum, level);
    }
    return colors;
}

int PhotoSet::getMask(const int index, const float fx, const float fy, const int level) const {
    return m_photos[index].Image::getMask(fx, fy, level);
}

int PhotoSet::getMask(const int index, const int ix, const int iy, const int level) const {
    return m_photos[index].Image::getMask(ix, iy, level);
}

int PhotoSet::getMask(const int index, const Vector4f& coord, const int level) const {
    return m_photos[index].getMask(coord, level);
}

int PhotoSet::getMask(const Vector4f& coord, const int level) const {
    for (int index = 0; index < m_nimages; ++index) {
        // ???
        // if a patch is legit, then it's either inside the mask, which return 255, or
        // is outside of the image plane, which return -1
        if (getMask(index, coord, level) == 0) {
            return 0;
        }
    }
    return -1;
}

int PhotoSet::getWidth(const int index, const int level) const {
    return m_photos[index].getWidth(level);
}

int PhotoSet::getHeight(const int index, const int level) const {
    return m_photos[index].getHeight(level);
}

Vector3f PhotoSet::project(const int index, const Vector4f& coord, const int level) const {
    return m_photos[index].project(coord, level);
}

float PhotoSet::computeDepth(const int index, const Vector4f& coord) const {
	return m_photos[index].computeDepth(coord);
}

int PhotoSet::image2index(const int image) const {
    map<int, int>::const_iterator pos = m_dict.find(image);
    if (pos == m_dict.end()) {
        return -1;
    }
    else {
        return pos->second;
    }
}
