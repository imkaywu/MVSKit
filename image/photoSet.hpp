//
//  photoSet.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/21/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef photoSet_hpp
#define photoSet_hpp

#include "photo.hpp"
#include <map>

using std::map;
using std::vector;

class PhotoSet {
public:
    PhotoSet();
    virtual ~PhotoSet();
    
    void init(const vector<int>& images, const string prefix, const int nimages, const int nillums, const int maxLevel, const int size, const int alloc);
    void free();
    void free(const int level);
    /*
    void getTex(const int index, const int level, const Vector2f& icoord, const Vector2f& xaxis, const Vector2f& yaxis, vector<Vector3f>& tex, const int isNorm = 1) const;
    void getTex(const int index, const int level, const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, vector<Vector3f>& tex, float& weight, const int isNorm = 1) const;
    void getTex(const int index, const int level, const Vector2f& icoord, const Vector2f& xaxis, const Vector2f& yaxis, vector<vector<Vector3f> >& texs, const int isNorm = 1) const;
    void getTex(const int index, const int level, const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, vector<vector<Vector3f> >& texs, float& weight, const int isNorm = 1) const;
     */
    
    Vector3f getColor(const int index, const Vector4f& coord, const int level) const;
    Vector3f getColor(const int index, const float fx, const float fy, const int level) const;
    Vector3f getColor(const int index, const int ix, const int iy, const int level) const;
    Vector3f getColor(const int index, const Vector4f& coord, const int level, const int illum);
    Vector3f getColor(const int index, const float fx, const float fy, const int level, const int illum);
    Vector3f getColor(const int index, const int ix, const int iy, const int level, const int illum);
    vector<Vector3f> getColorSeq(const int index, const Vector4f& coord, const int level);
    vector<Vector3f> getColorSeq(const int index, const float fx, const float fy, const int level);
    vector<Vector3f> getColorSeq(const int index, const int ix, const int iy, const int level);
    
    int getMask(const Vector4f& coord, const int level) const; // used only in the test
    int getMask(const int index, const Vector4f& coord, const int level) const;
    int getMask(const int index, const float fx, const float fy, const int level) const;
    int getMask(const int index, const int ix, const int iy, const int level) const;
    
    int getWidth(const int index, const int level) const;
    int getHeight(const int index, const int level) const;
    
    Vector3f project(const int index, const Vector4f& coord, const int level) const;

	float computeDepth(const int index, const Vector4f& coord) const;
    
    int image2index(const int image) const;
    int checkAngles(const Vector4f& coord, const vector<int>& indexes, const float minAngle, const float maxAngle, const int num) const;
    void setDistances();

    // image indexes
    vector<int> m_images;
    // photos from different viewpoints, under different illumination conditions
    vector<Photo, Eigen::aligned_allocator<Photo> > m_photos;
    // number of viewpoints
    int m_nimages;
    // number of illuminations
    int m_nillums;
    // window size
    int m_size;
    // root directory
    string m_prefix;
    // image-index dictionary
    map<int, int> m_dict;
    // pairwise distance based on camera center and viewing direction
    vector<vector<float> > m_distances;
};

#endif /* photoSet_hpp */
