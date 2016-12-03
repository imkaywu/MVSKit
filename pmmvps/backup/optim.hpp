//
//  optim.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef optim_hpp
#define optim_hpp

#include <vector>
#include "Eigen/Dense"
#include "patch.hpp"

using std::vector;
using Eigen::Vector3f;
using Eigen::Vector4f;

class PmMvps;

class Optim {
public:
    Optim(PmMvps& pmmvps);
    void init();
    
    //-----------------------------------------------------------------
    // Image related
    //-----------------------------------------------------------------
    // ???: fz / (fx + fy)
    // fz: the distance between the patch center and CoP
    // (fx + fy): the scale of the image plane
    float getUnit(const int index, const Vector4f& coord) const;
    // compute the x,y-axis of the patch, z-axis is the patch normal
    void getPAxes(const int index, const Vector4f& coord, const Vector4f& normal, Vector4f& pxaxis, Vector4f& pyaxis);
    void computeUnits(const Patch& patch, vector<int>& indexes, vector<float>& units, vector<Vector4f>& rays) const;
    void computeUnits(const Patch& patch, vector<float>& units) const;
    
    //-----------------------------------------------------------------
    // Optimization related
    //-----------------------------------------------------------------
    // find the visible views
    int preProcess(Patch& patch);
    // add potentially visible views based on the angle between normal and ray/line of sight
    void addImage(Patch& patch);
    // find the potentially visible views based on photo-consistency measure
    void constraintImages(Patch& patch, const float nccThreshold);
    // fix the reference image and sort the other m_tau - 1 images
    // isFixed: 0=>the first is not necessarily the reference, 1=>the first is the reference
    void sortImages(Patch& patch, const int isFixed = 1) const;
    // re-evaluate all the visible views and add the patch to the grid
    int postProcess(Patch& patch);
    // filter out invisible views based on the angle between normal and ray/line of sight
    // if the reference image is considered invisible, then this patch won't be added to the set
    void filterImagesByAngle(Patch& patch);
    // set the image with maximum photo-consistency score as the reference image
    void setRefImage(Patch& patch);
    // swap the reference image to the first position
    void swapImage(Patch& patch, const int image);
    
    //-----------------------------------------------------------------
    // Optimization
    //-----------------------------------------------------------------
    // refine patch
    void refinePatch(Patch& patch, const int time);
    int refinePatch(Patch& patch, const int time, const int ncc);
    static double cost_func(unsigned n, const double* x, double* grad, void* func_data);
    // transform from depth + normal to variables to be optimized
    void encode(const Vector4f& coord, double* const vect) const;
    void encode(const Vector4f& coord, const Vector4f& normal, double* const vect);
    // transform from optimized variables to depth + normal
    void decode(Vector4f& coord, Vector4f& normal, const double* const vect) const;
    void decode(Vector4f& coord, const double* const vect) const;
    // get texture from all visible viewpoints
    int getTex(const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, const int index, const int size, vector<Vector3f>& tex) const;
    int getTex(const Vector4f& coord, const Vector4f& pxaxis, const Vector4f& pyaxis, const Vector4f& pzaxis, const int index, const int size, vector<vector<Vector3f> >& texs) const;
    // test if it's possible to get texture
    int getTexSafe(const int index, const int size, const Vector3f& center, const Vector3f& dx, const Vector3f& dy, const int level) const;
    // level different to scale
    float myPow2(int levelDiff) const;
    // compute m_weights
    void computeWeights(const Patch& patch);
    // compute the incc score of the patch
    float computeINCC(const Vector4f& coord, const Vector4f& normal, const vector<int>& indexes, const int isRobust);
    float computeINCC(const Vector4f& coord, const Vector4f& normal, const vector<int> indexes, const Vector4f& pxaxis, const Vector4f& pyaxis, const int isRobust);
    float computeINCC2(const Vector4f& coord, const Vector4f& normal, const vector<int> indexes, const Vector4f& pxaxis, const Vector4f& pyaxis, const int isRobust);
    // compute the incc score between the reference image and all the other visible views
    void setINCCs(const Patch& patch, vector<float>& inccs, const vector<int>& indexes, const int isRobust);
    // compute the incc score between any two visible views
    void setINCCs(const Patch& patch, vector<vector<float> >& inccs, const vector<int>& indexes, const int isRobust);
    // normalize the grabbed texture
    static void normalize(vector<Vector3f>& tex);
    // normalize only scale, not used
    static void normalize(vector<vector<Vector3f> >& texs, const int sz);
    // cross correlation
    float dot(const vector<Vector3f>& tex0, const vector<Vector3f>& tex1) const;
    // sum/mean of absolute difference
    float ssd(const vector<Vector3f>& tex0, const vector<Vector3f>& tex1) const;
    // robustify the incc score
    static float robustincc(const float incc);
    // unrobustify the incc score
    static float unrobustincc(const float rincc);
    
protected:
    // compute all the camera axes
    void setAxesScales();
    
    static Optim* m_inst;
    PmMvps& m_pmmvps;
    
    // Axes of the camera-centered coordinate system
    vector<Vector3f> m_xaxes;
    vector<Vector3f> m_yaxes;
    vector<Vector3f> m_zaxes;
    // image plane scale
    // not sure what is this
    vector<float> m_ipscales;
    
    // center of projection
    Vector4f m_center;
    // ray/line of sight
    Vector4f m_ray;
    vector<int> m_indexes;
    float m_dscale;
    float m_ascale;
    
    // Grabbed texture for all visible views
    vector<vector<Vector3f> > m_texs; // m_nillums x (wsize*wsize)
    vector<vector<vector<Vector3f> > > m_textures; // m_nviews x m_nillums x (wsize*wsize)
    // Weights
    vector<float> m_weights;
    
};

#endif /* optim_hpp */
