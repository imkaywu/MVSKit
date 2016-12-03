//
//  image_set.cpp
//  PMMVPS
//
//  Created by KaiWu on Nov/1/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include <fstream>
#include "image_set.hpp"

using std::ifstream;

ImageSet::ImageSet() {
}

ImageSet::~ImageSet() {
}

void ImageSet::init(const string prefix, const int image, const string mname, const int numIllum, const int maxLevel) {
    m_nillums = numIllum;
    m_maxLevel = maxLevel;
    
    for (int n = 0; n < m_nillums; ++n) {
        char img0[1024], img1[1024], img2[1024], img3[1024];
        sprintf(img0, "%svisualize/%04d%04d.ppm", prefix.c_str(), image, n);
        sprintf(img1, "%svisualize/%04d%04d.jpg", prefix.c_str(), image, n);
        sprintf(img2, "%svisualize/%04d%04d.png", prefix.c_str(), image, n);
        sprintf(img3, "%svisualize/%04d%04d.tiff", prefix.c_str(), image, n);
        
        if (ifstream(img0) || ifstream(img1)
#if defined(PMMVPS_HAVE_PNG)
            || ifstream(img2)
#endif
#if defined(PMMVPS_HAVE_TIFF)
            || ifstream(img3)
#endif
            ) {
            char iname[1024];
            sprintf(iname, "%svisualize/%04d%04d", prefix.c_str(), image, n);
            m_images[n].init(iname, mname, maxLevel);
            // m_images[n].init(iname, maxLevel);
        }
    }
}

void ImageSet::completeName(const string& lhs, string& rhs, const int isColor) {
    // what's the point of doing this
    if (5 <= lhs.length() && lhs[lhs.length() - 4] == '.') {
        rhs = lhs;
        return;
    }
    
    // ppm jpg
    if (isColor) {
        string stmp0 = lhs + ".ppm";    string stmp1 = lhs + ".jpg";
        string stmp2 = lhs + ".png";    string stmp3 = lhs + ".tiff";
        
        if (ifstream(stmp0.c_str()))
            rhs = stmp0;
        else if (ifstream(stmp1.c_str()))
            rhs = stmp1;
#if defined(PMMVPS_HAVE_PNG)
        else if (ifstream(stmp2.c_str()))
            rhs = stmp2;
#endif
#if defined(PMMVPS_HAVE_TIFF)
        else if (ifstream(stmp3.c_str()))
            rhs = stmp3;
#endif
        else
            rhs = lhs;
    }
    // pgm pbm
    else {
        string stmp0 = lhs + ".pgm";    string stmp1 = lhs + ".pbm";
        
        if (ifstream(stmp0.c_str()))
            rhs = stmp0;
        else if (ifstream(stmp1.c_str()))
            rhs = stmp1;
        else
            rhs = lhs;
    }
}

void ImageSet::alloc(const int fast, const int filter) {
    if (m_alloc == 1) {
        return;
    }
    
    m_masks.resize(m_maxLevel);
    
    if (!m_mname.empty()) {
        if (readPGMImage(m_mname, m_masks[0], m_images[0].m_widths[0], m_images[0].m_heights[0], 0) ||
            readPBMImage(m_mname, m_masks[0], m_images[0].m_widths[0], m_images[0].m_heights[0], 0)) {
            cerr << "Read mask: " << m_mname << endl;
            for (int i = 0; i < static_cast<int>(m_masks[0].size()); ++i) {
                if (127 < static_cast<int>(m_masks[0][i])) {
                    m_masks[0][i] = (unsigned char) 255;
                }
                else {
                    m_masks[0][i] = (unsigned char) 0;
                }
            }
        }
        else {
            m_mname = "";
        }
    }
    
    buildMaskPyramid(filter);
    m_alloc = 1;
}

void ImageSet::free(const int freeLevel) {
    for (int l = 0; l < freeLevel; ++l) {
        if (!m_masks.empty()) {
            vector<unsigned char>().swap(m_masks[l]);
        }
    }
}

void ImageSet::free() {
    if(m_alloc == 1) {
        m_alloc = 0;
    }
    else {
        return;
    }
    vector<vector<unsigned char> >().swap(m_masks);
}

int ImageSet::readPBMImage(const string file, vector<unsigned char>& image, int& width, int& height, const int fast) {
    
    if (file.substr(file.length() - 3, file.length()) != "pbm") {
        return -1;
    }
    
    ifstream ifstr;
    ifstr.open(file.c_str());
    if (!ifstr.is_open()) {
        return -1;
    }
    string header;
    unsigned char byte;
    
    ifstr >> header;
    ifstr.read((char*)&byte, sizeof(unsigned char));
    
    if (header != "P4") {
        cerr << "Only accept binary pbm format: " << file << endl;
        return -1;
    }
    
    while (1) {
        ifstr.read((char*)&byte, sizeof(unsigned char));
        ifstr.putback(byte);
        if (byte == '#') {
            char buffer[1024];
            ifstr.getline(buffer, 1024);
        }
        else {
            break;
        }
    }
    ifstr >> width >> height;
    ifstr.read((char*)&byte, sizeof(unsigned char));
    
    image.clear();
    if (fast) {
        ifstr.close();
        return 0;
    }
    int reso = width * height;
    if (reso % 8 != 0) {
        reso++;
    }
    
    int count = 0;
    for (int i = 0; i < reso; ++i) {
        ifstr.read((char*)&byte, sizeof(unsigned char));
        for (int j = 0; j < 8; ++j) {
            if (byte >> 7) {
                image.push_back((unsigned char)0);
            }
            else {
                image.push_back((unsigned char)255);
            }
            count++;
            byte <<= 1;
            if (count == width * height) {
                break;
            }
        }
    }
    ifstr.close();
    return 0;
}

int ImageSet::readPGMImage(const string file, vector<unsigned char>& image, int& width, int& height, const int fast) {
    if (file.substr(file.length() - 3, file.length()) != "pgm") {
        return -1;
    }
    
    ifstream ifstr;
    ifstr.open(file.c_str());
    if (!ifstr.is_open()) {
        cerr << "Cannot open a file: " << file << endl;
        return -1;
    }
    
    string header;
    unsigned char byte;
    
    ifstr >> header;
    ifstr.read((char*)&byte, sizeof(unsigned char));
    if (header != "P5") {
        cerr << "Only accept binary pgm format" << file << " " << header << endl;
        return -1;
    }
    
    while (1) {
        ifstr.read((char*)&byte, sizeof(unsigned char));
        ifstr.putback(byte);
        if (byte == '#') {
            char buffer[1024];
            ifstr.getline(buffer, 1024);
        }
        else {
            break;
        }
    }
    int itmp;
    ifstr >> width >> height >> itmp;
    ifstr.read((char*)&byte, sizeof(unsigned char));
    
    image.clear();
    if (fast) {
        ifstr.close();
        return 0;
    }
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            ifstr.read((char*)&byte, sizeof(unsigned char));
            image.push_back(byte);
        }
    }
    ifstr.close();
    
    return 0;
}

void ImageSet::buildMaskPyramid(const int filter) {
    for (int level = 1; level < m_maxLevel; ++level) {
        const int sz = getWidth(level) * getHeight(level);
        
        m_masks[level].resize(sz);
        for (int y = 0; y < getHeight(level); ++y) {
            const int ys[2] = {2 * y, std::min(getHeight(level - 1), 2 * y + 1)};
            for (int x = 0; x < getWidth(level);  ++x) {
                const int xs[2] = {2 * x, std::min(getWidth(level - 1), 2 * x + 1)};
                int inside = 0, outside = 0;
                for (int j = 0; j < 2; ++j) {
                    for (int i = 0; i < 2; ++i) {
                        const int index = ys[j] * getWidth(level - 1) + xs[i];
                        if (m_masks[level - 1][index]) {
                            ++inside;
                        }
                        else {
                            ++outside;
                        }
                    }
                }
                const int index = y * getWidth(level) + x;
                if (0 < inside) {
                    m_masks[level][index] = (unsigned char)255;
                }
                else
                    m_masks[level][index] = (unsigned char)0;
            }
        }
    }
}

Vector3f ImageSet::getColor(const float fx, const float fy, const int level, const int illum) const {
    return m_images[illum].getColor(fx, fy, level);
}

Vector3f ImageSet::getColor(const int ix, const int iy, const int level, const int illum) const {
    return m_images[illum].getColor(ix, iy, level);
}

int ImageSet::getMask(const float fx, const float fy, const int level) const {
    if (m_alloc != 1) {
        cerr << "Image data not allocated" << endl;
        exit(1);
    }
    
    if (m_masks[level].empty()) {
        return -1;
    }
    
    const int ix = (int)floorf(fx + 0.5f);
    const int iy = (int)floorf(fy + 0.5f);
    return getMask(ix, iy, level);
}

int ImageSet::getMask(const int ix, const int iy, const int level) const {
    if (m_alloc != 1) {
        cerr << "Image data not allocated" << endl;
        exit(1);
    }
    
    if (m_masks[level].empty()) {
        return -1;
    }
    
    if (ix < 0 || getWidth(level) <= ix ||
        iy < 0 || getHeight(level) <= iy) {
        return -1;
    }
    
    const int index = iy * getWidth(level) + ix;
    return m_masks[level][index];
}

int ImageSet::getWidth(const int level) const {
    return m_images[0].m_widths[level];
}

int ImageSet::getHeight(const int level) const {
    return m_images[0].m_heights[level];
}