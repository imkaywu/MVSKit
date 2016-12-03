//
//  image_set.hpp
//  PMMVPS
//
//  Created by KaiWu on Nov/1/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef image_set_hpp
#define image_set_hpp

#include "image.hpp"

class ImageSet {
public:
    ImageSet();
    virtual ~ImageSet();
    
    virtual void init(const string prefix, const int image, const string mname, const int numIllum, const int maxLevel);
    void completeName(const string& lhs, string& rhs, const int colour);
    
    // allocate and free memories
    void alloc(const int fast = 0, const int filter = 0);
    void free();
    void free(const int freeLevel);
    
    void buildMaskPyramid(const int filter);
    
    static int readPBMImage(const string file, vector<unsigned char>& image, int& width, int& height, const int fast); // Portable BitMap
    static int readPGMImage(const string file, vector<unsigned char>& image, int& width, int& height, const int fast); // Portable GreyMap
    
//    Vector3f getTex(const )
    
    Vector3f getColor(const float fx, const float fy, const int level, const int illum = 0) const;
    Vector3f getColor(const int ix, const int iy, const int level, const int illum = 0) const;
    
    int getMask(const float fx, const float fy, const int level) const;
    int getMask(const int ix, const int iy, const int level) const;
    
    int getWidth(const int level) const;
    int getHeight(const int level) const;
    
    // number of levels
    int m_maxLevel;
    // number of illuminations
    int m_nillums;
    // a set of images from the single viewpoint, under different illuminations
    vector<Image> m_images;
    // allocate/free images
    // 0: nothing allocated, 1: memory allocated;
    int m_alloc;
    
    // name of the mask
    string m_mname;
    // a pyramid of masks
    vector<vector<unsigned char> > m_masks;
};

#endif /* image_set_hpp */
