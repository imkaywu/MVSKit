//
//  camera.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/20/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef camera_hpp
#define camera_hpp

#include <iostream>
#include <vector>
#include "Eigen/Dense"

using std::vector;
using std::string;
using Eigen::Matrix3f;
using Eigen::Matrix3Xf;
using Eigen::Matrix4f;
using Eigen::Vector3f;
using Eigen::Vector4f;

class Camera
{
public:
    Camera();
    virtual ~Camera();
    
    virtual void init(const string cname, const int maxLevel);
    // update all camera-related parameters: projection matrices, and various axes
    void updateCamera();
    // update projection matrices from intrinsics and extrinsics
    void updateProjection();
    // set projection matrics from intrinsics and extrinsics
    void setProjection(const vector<float>& intrinsics, const vector<float>& extrinsics, Matrix3Xf& projection, const int txtType);
    // write projection to file
    void writeProjection(const string file);
    // set K
    void setK(Matrix3f& K) const;
    // set Rt
    void setRt(Matrix4f& Rt) const;
    // set R
    void setR(Matrix3f& R) const;
    // projection 2 quaternion
    static void proj2quat(Matrix4f& proj, float q[6]);
    // quaternion to projection
    static void quat2proj(const float q[6], Matrix4f& proj);
    
    Vector3f project(const Vector4f& coord, const int level) const;
    Vector4f unproject(const Vector3f& icoord, const int level) const;
	float computeDepth(const Vector4f& coord) const;
    
    float getScale(const Vector4f& coord, const int level) const;
    // get patch axes
    void getPAxes(const Vector4f& coord, const Vector4f& normal, Vector4f& pxaxis, Vector4f& pyaxis, const int level = 0) const;
    void setAxisScale(const float axisScale);
    
    // text file name of camera parameters
    string m_cname;
    // optical center
    Vector4f m_center;
    // optical axis
    Vector4f m_oaxis;
    // x-axis of the camera-centered coordinate system
    Vector3f m_xaxis;
    // y-axis of the camera-centered coordinate system
    Vector3f m_yaxis;
    // z-axis of the camera-centered coordinate system
    Vector3f m_zaxis;
    // 3x4 projection matrix for each level
    vector<Matrix3Xf> m_projections;
    // intrinsic and extrinsic camera parameters
    vector<float> m_intrinsics;
    vector<float> m_extrinsics;
    // image plane scale (fx + fy), used to compute the projection/scene-image scale
    float m_ipscale;
    // camera pamameter type
    int m_txtType;
protected:
    int m_maxLevel;
    float m_axisScale;
    
    Vector4f getCameraCenter() const;
};

#endif /* camera_hpp */
