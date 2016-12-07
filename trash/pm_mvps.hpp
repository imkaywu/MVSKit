/********************************************
difference between mex array and c/c++ array
mex array: column first
|1 4 7|
|2 5 8|
|3 6 9|
c/c++ array: row first
|1 2 3|
|4 5 6|
|7 8 9|
a pixel (i,j) has index of i * h + j in matlab, but j * w + i in c/c++
i, j is in (0, h - 1) x (0, w - 1)
*********************************************/

#ifndef pm_stereo_hpp
#define pm_stereo_hpp

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <random>

#include "Eigen/Dense"
#include "include/mvg.hpp"

using std::string;
using std::vector;

using Eigen::Matrix3Xd;
using Eigen::Vector2d;
using Eigen::Vector3d;

class BITMAP {
public:
    int w, h;
    int *data;
    BITMAP(int w_, int h_) :w(w_), h(h_) { data = new int[w*h]; }
    ~BITMAP() { delete[] data; }
    int *operator[](int y) { return &data[y*w]; } // row-first storage
};


void check_im();
BITMAP *load_bitmap(const char *filename);
void save_bitmap(BITMAP *bmp, const char *filename);
double *load_txt(const char* fname);


void patchmatch_stereo_main(char *image_base, char *mask_base, int num_img, double *matchl2r, double *matchr2l, double *nl, double *nr, vector<Matrix3Xd> &P, double *depth);
void patchmatch_stereo(vector<BITMAP *> &a, vector<BITMAP *> &b, vector<BITMAP *> &mask_ref, vector<BITMAP *> &mask_tar, double *matchl2r, double *matchr2l, double *nl, double *nr, vector<Matrix3Xd> &P, double *depth);
double dist(vector<BITMAP *> &a, vector<BITMAP *> &b, double *nl, double *nr, double ax, double ay, double bx, double by);
void improve_guess(vector<BITMAP *> &a, vector<BITMAP *> &b, double *nl, double *nr, double ax, double ay, double bx, double by, double xbest, double ybest, double dbest);
double uniform_random(double range);


#endif /* pm_stereo_h */
