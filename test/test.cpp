//
//  test.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/19/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "CImg.h"
#include "../pmmvps/pmmvps.hpp"

//using namespace std;
//using namespace Eigen;
//using namespace cimg_library;

int main()
{
    // test of removing an element from a vector
//    int dim = 6;
//    vector<Ppatch> ppatches(dim);
//    for (int i = 0; i < dim; ++i) {
//        Ppatch ppatch(new Patch());
//        ppatches[i] = ppatch;
//        cerr << ppatches[i] << endl;
//    }
//    
//    Ppatch ppatch_removed = ppatches[2];
//    ppatches.erase(std::remove(ppatches.begin(), ppatches.end(), ppatch_removed), ppatches.end());
//    for (int i = 0; i < static_cast<int>(ppatches.size()); ++i) {
//        cerr << ppatches[i] << endl;
//    }
    
    // test of reference of vector of vector
//    int dim = 6;
//    int count = 1;
//    vector<vector<int> > vecs(dim);
//    for (int i = 0; i < dim; ++i) {
//        vecs[i].resize(dim);
//        for (int j = 0; j < dim; ++j) {
//            vecs[i][j] = count++;
//        }
//    }
//    
//    for (int i = 0; i < dim; ++i) {
//        for (int j = 0; j < dim; ++j) {
//            std::cerr << vecs[i][j] << " ";
//        }
//        std::cerr << std::endl;
//    }
//    
//    vector<int>& vec = vecs[2];
//    for (int i = 0; i < dim; ++i)
//        vec[i] = 0;
//    
//    for (int i = 0; i < dim; ++i) {
//        for (int j = 0; j < dim; ++j) {
//            std::cerr << vecs[i][j] << " ";
//        }
//        std::cerr << std::endl;
//    }
    
    //----------------- test of reading jpeg
//    cimg::imagemagick_path("/opt/local/bin/convert");
//    CImg<unsigned char> image("1.jpg"), visu(500,400,1,3,0);
//    for (int x = 0; x < 10; ++x)
//    {
//        for (int y = 0; y < 10; ++y)
//        {
//            for (int c = 0; c < image.spectrum(); ++c)
//            {
//                cout << static_cast<int>(image(x, y, 0, c)) << endl;
//            }
//        }
//    }
//    const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
//    image.blur(2.5);
//    CImgDisplay main_disp(image,"Click a point"), draw_disp(visu,"Intensity profile"); // why it doesn't work?
//    while (!main_disp.is_closed() && !draw_disp.is_closed()) {
//        main_disp.wait();
//        if (main_disp.button() && main_disp.mouse_y()>=0) {
//            const int y = main_disp.mouse_y();
//            visu.fill(0).draw_graph(image.get_crop(0,y,0,0,image.width()-1,y,0,0),red,1,1,0,255,0);
//            visu.draw_graph(image.get_crop(0,y,0,1,image.width()-1,y,0,1),green,1,1,0,255,0);
//            visu.draw_graph(image.get_crop(0,y,0,2,image.width()-1,y,0,2),blue,1,1,0,255,0).display(draw_disp);
//        }
//    }
    
    //----------------- generate camera parameter files
    
//    Camera camera;
//    camera.m_intrinsics.resize(6);
//    camera.m_extrinsics.resize(6);
//    camera.m_intrinsics[0] = 765.702941895f; // fx;
//    camera.m_intrinsics[1] = 765.702941895f; // fy;
//    camera.m_intrinsics[2] = 0.0200504438012f; // skew;
//    camera.m_intrinsics[3] = 320.0f; // cx;
//    camera.m_intrinsics[4] = 240.0f; // cy;
//    camera.m_intrinsics[5] = 0.0f; // not used in m_txtType == 2
//    
//    // rotation
//    Matrix3f R;
//    R << 0.98292940855, -0.030582388863, 0.181423678994,
//    0.0355352722108, 0.999077498913, -0.0241120066494,
//    -0.180518895388, 0.0301473364234, 0.983109354973;
//    
//    // translation
//    Vector3f t;
//    t << -0.959477365017, 0.0510520711541, 0.251982748508;
//    
//    // rotation and translation
//    Matrix4f Rt;
//    Rt.block(0, 0, 3, 3) = R;
//    Rt.block(0, 3, 3, 1) = t;
//    Rt(3, 3) = 1.0f;
//    
//    // generate type 2 camera parameter files
//    float q[6];
//    camera.proj2quat(Rt, q);
//
//    for (int i = 0; i < 6; ++i) {
//        camera.m_extrinsics[i] = q[i];
//    }
//    camera.m_txtType = 2;
//    camera.writeProjection("00000001.txt");
//    
//    // generate type 1 camear parameter files
//    Matrix3f K3;
//    camera.m_txtType = 2;
//    camera.setK(K3);
//    Matrix4f K = Matrix4f::Zero();
//    K.block(0, 0, 3, 3) = K3;
//    K(3, 3) = 1.0f;
//    
//    float params[12];
//    Matrix4f proj = K * Rt;
//    for (int y = 0; y < 3; ++y) {
//        for (int x = 0; x < 4; ++x) {
//            params[4 * y + x] = proj(y, x);
//        }
//    }
//    
//    for (int i = 0; i < 6; ++i) {
//        camera.m_intrinsics[i] = params[i];
//        camera.m_extrinsics[i] = params[i + 6];
//    }
//    camera.m_txtType = 0;
//    camera.writeProjection("00000000.txt");
    
    
    //----------------- test of reading depths and normals
    
    Option option;
    option.init("C:/Users/Admin/Documents/Data/PMMVPS/", "option");
    
    PmMvps pmmvps;
    pmmvps.init(option);
    
    pmmvps.run();
    
    return 0;
}