//
//  photo.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/20/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "photo.hpp"

Photo::Photo() {
}

Photo::~Photo() {
}

void Photo::init(const string cname, const string iname, const string mname, const int nillums, const int maxLevel) {
    Camera::init(cname, maxLevel);
    Image::init(iname, mname, nillums, maxLevel);
}

Vector3f Photo::getColor(const Vector4f& coord, const int level) const {
    const Vector3f icoord = project(coord, level);
    return Image::getColor(icoord(0), icoord(1), level);
}

Vector3f Photo::getColor(const float fx, const float fy, const int level) const {
    return Image::getColor(fx, fy, level);
}

Vector3f Photo::getColor(const int ix, const int iy, const int level) const {
    return Image::getColor(ix, iy, level);
}

Vector3f Photo::getColor(const Vector4f& coord, const int level, const int illum) const {
    const Vector3f icoord = project(coord, level);
    return Image::getColor(icoord(0), icoord(1), level, illum);
}

Vector3f Photo::getColor(const float fx, const float fy, const int level, const int illum) const {
    return Image::getColor(fx, fy, level, illum);
}

Vector3f Photo::getColor(const int ix, const int iy, const int level, const int illum) const {
    return Image::getColor(ix, iy, level, illum);
}

int Photo::getMask(const Vector4f &coord, const int level) const {
    if (m_masks[level].empty()) {
        return -1;
    }
    
    const Vector3f icoord = project(coord, level);
    return Image::getMask(icoord(0), icoord(1), level);
}
/*
void Photo::getTex(const int level, const Vector2f& icoord, const Vector2f& xaxis, const Vector2f& yaxis, const int size, vector<Vector3f>& tex, const int isNorm) const {
    getTex(0, level, icoord, xaxis, yaxis, size, tex, 1);
}



void Photo::getTex(const int level, const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, const int size, vector<Vector3f>& tex, float& weight, const int isNorm) const {
    const int scale = 0x0001 << level;
    
    const Vector3f icoord3 = project(coord, level);
    const Vector2f icoord(icoord3(0), icoord3(1));
    
    const Vector3f xaxis3 = project(coord + pxaxis * scale, level) - icoord3;
    const Vector2f xaxis(xaxis3(0), xaxis3(1));
    
    const Vector3f yaxis3 = project(coord + pyaxis * scale, level) - icoord3;
    const Vector2f yaxis(yaxis3(0), yaxis3(1));
    
    getTex(0, level, icoord, xaxis, yaxis, size, tex, 1);
    
    Vector4f ray = m_center - coord;
    ray /= ray.norm();
    weight = std::max(0.0f, pzaxis.dot(ray));
}

void Photo::getTex(const int illum, const int level, const Vector2f& icoord, const Vector2f& xaxis, const Vector2f& yaxis, const int size, vector<Vector3f>& tex, const int isNorm) const {
    const int margin = size / 2;
    
    // check boundary condition
    const float minx = icoord[0] - size * fabs(xaxis(0)) - size * fabs(yaxis(0));
    const float maxx = icoord[0] + size * fabs(xaxis(0)) + size * fabs(yaxis(0));
    const float miny = icoord[1] - size * fabs(xaxis(1)) - size * fabs(yaxis(1));
    const float maxy = icoord[1] + size * fabs(xaxis(1)) + size * fabs(yaxis(1));
    
    tex.clear();
    if (minx < 0 || getWidth(level) - 1 <= maxx ||
        miny < 0 || getHeight(level) - 1 <= maxy) {
        return;
    }
    
    for (int y = -margin; y <= margin; ++y) {
        for (int x = -margin; x <= margin; ++x) {
            Vector2f icoordSamp = icoord + x * xaxis + y * yaxis;
            tex.push_back(getColor(icoordSamp(0), icoordSamp(1), level, illum));
        }
    }
    
    if (isNorm) {
        normalize(tex);
    }
}

void Photo::getTex(const int illum, const int level, const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, const int size, vector<Vector3f>& tex, float& weight, const int isNorm) const {
    const int scale = 0x0001 << level;
    
    const Vector3f icoord3 = project(coord, level);
    const Vector2f icoord(icoord3(0), icoord3(1));
    
    const Vector3f xaxis3 = project(coord + pxaxis * scale, level) - icoord3;
    const Vector2f xaxis(xaxis3(0), xaxis3(1));
    
    const Vector3f yaxis3 = project(coord + pyaxis * scale, level) - icoord3;
    const Vector2f yaxis(yaxis3(0), yaxis3(1));
    
    getTex(illum, level, icoord, xaxis, yaxis, size, tex, isNorm);
    
    Vector4f ray = m_center - coord;
    ray /= ray.norm();
    weight = std::max(0.0f, pzaxis.dot(ray));
}

void Photo::normalize(vector<Vector3f>& tex) {
    // compute average
    Vector3f ave;
    int sz = static_cast<int>(tex.size());
    for (int i = 0; i < sz; ++i) {
        ave += tex[i];
    }
    ave /= sz;
    
    for (int i = 0; i < sz; ++i) {
        tex[i] -= ave;
    }
    
    // compute variance
    float ssd = 0.0f;
    for (int i = 0; i < sz; ++i) {
        ssd += tex[i].squaredNorm();
    }
    ssd /= (3 * sz);
    ssd = sqrtf(ssd);
    if(ssd == 0.0f)
        ssd = 1.0f;
    
    for (int i = 0; i < sz; ++i) {
        tex[i] /= ssd;
    }
}
*/