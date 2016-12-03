//
//  photo.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/20/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef photo_hpp
#define photo_hpp

#include <vector>
#include "camera.hpp"
#include "image.hpp"

using std::vector;
using Eigen::Vector2f;

class Photo : public Camera, public Image//, public ImageSet
{
public:
    Photo();
    virtual ~Photo();
    
    // naming convention for cameras and masks
    // xxxxxxxx.FORMAT (each x is a number)
    // naming convention for images
    // aaaabbbb.FORMAT (each a, b is a number)
    // aaaa represents the index of the view, which is the same as the name of the corresponding camera and mask,
    // bbbb represents the index of the illumination
    virtual void init(const string cname, const string iname, const string mname, const int nillums, const int maxLevel = 1);
    
    Vector3f getColor(const Vector4f& coord, const int level) const;
    Vector3f getColor(const float fx, const float fy, const int level) const;
    Vector3f getColor(const int ix, const int iy, const int level) const;
    Vector3f getColor(const Vector4f& coord, const int level, const int illum) const;
    Vector3f getColor(const float fx, const float fy, const int level, const int illum) const;
    Vector3f getColor(const int ix, const int iy, const int level, const int illum) const;
    
    int getMask(const Vector4f& coord, const int level) const;
    /*
    void getTex(const int level, const Vector2f& icoord, const Vector2f& xaxis, const Vector2f& yaxis, const int size, vector<Vector3f>& tex, const int isNorm = 1) const;
    void getTex(const int level, const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, const int size, vector<Vector3f>& tex, float& weight, const int isNorm) const;
    void getTex(const int illum, const int level, const Vector2f& icoord, const Vector2f& xaxis, const Vector2f& yaxis, const int size, vector<Vector3f>& tex, const int isNorm = 1) const;
    void getTex(const int illum, const int level, const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, const int size, vector<Vector3f>& tex, float& weight, const int isNorm) const;
     */
    
    // to be implemented
    static float idot();
    static float idotC();
    static float ssd();
    static void normalize(vector<Vector3f>& tex); // it's invoked in a const function, has to be static
};

#endif /* photo_hpp */
