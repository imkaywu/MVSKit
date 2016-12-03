//
//  depth_normal_init.cpp
//  PMMVPS
//
//  Created by KaiWu on Nov/2/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "depth_normal_init.hpp"
#include "pmmvps.hpp"
#include <fstream>
extern "C" {
#include "../io/io_file.h"
}

using std::ofstream;

DepthNormInit::DepthNormInit(PmMvps& pmmvps) : m_pmmvps(pmmvps) {
}

DepthNormInit::~DepthNormInit() {
}

void DepthNormInit::init(const string prefix, const int nfiles) {
    m_prefix = prefix;
    m_nplys = nfiles;
}

void DepthNormInit::createPatches() {
    const int isTest = 1;
    if (isTest) {
        m_pmmvps.m_patchManager.readPatches();
    }
    else {
    // read coords from dname
    vector<Vector3f> coords;
    readDepths(coords);

    // read normals from nname, m_views x (n_width * n_height)
    vector<vector<Vector3f> >normals;
    normals.resize(m_nplys - 1);
    readNormals(normals);
    
    // generate patches from depth and normal data
    int ncoords = 0;
    for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
        // initalize the patch with coordinate
        Ppatch ppatch(new Patch());
        ppatch->m_coord << coords[i](0), coords[i](1), coords[i](2), 1.0f;
        
        // get normal from all the potentially visible views
        vector<int> images;
        Vector3f normal3 = Vector3f::Zero();
        for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
            Vector3f icoord = m_pmmvps.m_photoSet.project(image, ppatch->m_coord, 0);
            int x = (int)(floorf(icoord(0) + 0.5f));
            int y = (int)(floorf(icoord(1) + 0.5f));
            if (m_pmmvps.m_photoSet.getMask(image, x, y, 0) <= 0) {
                continue;
            }
            // get the corresponding normal
            int index = y * m_pmmvps.m_photoSet.getWidth(image, 0) + x;
            normal3 += normals[image][index];
            images.push_back(image);
        }
        
        if (images.size() < 2 || normal3.norm() == 0.0f) {
            continue;
        }
        // set m_images and m_nimages
        ppatch->m_images = images;
        ppatch->m_nimages = static_cast<int>(ppatch->m_images.size());
        // set m_normal
        normal3 /= (float)images.size();
        normal3 /= normal3.norm();
        ppatch->m_normal << normal3(0), normal3(1), normal3(2), -ppatch->m_coord.head(3).dot(normal3);
        // organize the order of the visible views
        m_pmmvps.m_optim.sortImages(*ppatch, 0);
        // reset m_grids
        m_pmmvps.m_patchManager.setGrids(*ppatch); // m_images is fixed
        // set m_tmp
//        ppatch->m_tmp = ppatch->score2(m_pmmvps.m_nccThreshold); // a potential problem, m_ncc is not computed
        // add patch to m_pgrids and update m_dpgrids
        m_pmmvps.m_patchManager.addPatch(ppatch);
        ++ncoords;
    }
    // write patches to a PLY file
//    char pname[1024];
//    sprintf(pname, "initial_patches");
//    m_pmmvps.m_patchManager.writePatches(m_pmmvps.m_prefix, pname, true, false, false);
    
    coords.clear();
    normals.clear();
    }
}

void DepthNormInit::readDepths(vector<Vector3f>& coords) {
    char dname[1024];
    sprintf(dname, "%sply/%08d.ply", m_prefix.c_str(), 0);
    
    int ncoords;
    ply_header_read(dname, &ncoords, NULL, NULL);
    
    double *points = (double *)malloc(DIM * ncoords * sizeof(double));
    ply_read_1(dname, points, NULL, NULL, NULL, NULL);
    
    coords.resize(ncoords);
    for (int i = 0; i < ncoords; ++i) {
        coords[i] << (float)points[i * DIM + 0], (float)points[i * DIM + 1], (float)points[i * DIM + 2];
    }
    
    free(points);
}

void DepthNormInit::readNormals(vector<vector<Vector3f> >& normals) {
    for (int i = 0; i < m_nplys - 1; ++i) {
        char nname[1024];
        sprintf(nname, "%sply/%08d.ply", m_prefix.c_str(), i + 1);
        
        int ncoords;
        ply_header_read(nname, &ncoords, NULL, NULL);
        
        double *points = (double*)malloc(DIM * ncoords * sizeof(double));
        double *norms = (double*)malloc(DIM * ncoords * sizeof(double));
        ply_read_1(nname, points, norms, NULL, NULL, NULL);
        
        int width = m_pmmvps.m_photoSet.getWidth(i, 0);
        int height = m_pmmvps.m_photoSet.getHeight(i, 0);
        normals[i].resize(width * height);
        for (int n = 0; n < ncoords; ++n) {
            int x = (int)points[n * DIM + 0];
            int y = (int)points[n * DIM + 1];
            int index = y * width + x;
            Vector3f normal3;
            normal3 << (float)norms[n * DIM + 0], (float)norms[n * DIM + 1], (float)norms[n * DIM + 2];
            Matrix3f R;
            m_pmmvps.m_photoSet.m_photos[i].setR(R);
            normals[i][index] = R * normal3; // need to test
        }
        
        free(points);
        free(norms);
    }
}
// not used
void DepthNormInit::readDepthsNorms(vector<Vector3f>& coords, vector<Vector3f>& normals) {
    char pname[1024];
    sprintf(pname, "%sply/%08d.ply", m_prefix.c_str(), 0);
    
    int ncoords;
    ply_header_read(pname, &ncoords, NULL, NULL);
    
    double *points = (double*)malloc(DIM * ncoords * sizeof(double));
    double *norms = (double*)malloc(DIM * ncoords * sizeof(double));
    ply_read_1(pname, points, norms, NULL, NULL, NULL);
    
    int width = m_pmmvps.m_photoSet.getWidth(0, 0);
    for (int n = 0; n < ncoords; ++n) {
        int x = (int)points[n * DIM + 0];
        int y = (int)points[n * DIM + 1];
        int index = y * width + x;
        coords[index] << (float)points[n * DIM + 0], (float)points[n * DIM + 1], (float)points[n * DIM + 2];
        normals[index] << (float)norms[n * DIM + 0], (float)norms[n * DIM + 1], (float)norms[n * DIM + 2];
    }
    
    free(points);
    free(norms);
}

// another way to combine depth and normal together
/*
// compute the pairwise incc
float inccThreshold = INT_MAX / 2;
vector<vector<float> > inccs;
m_pmmvps.m_optim.setINCCs(*ppatch, inccs, ppatch->m_images, 1);
vector<int> visibs = ppatch->m_images;
ppatch->m_images.resize(2);

for (int i = 0; i < static_cast<int>(inccs.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(inccs[i].size()); ++j) {
        if (inccs[i][j] < inccThreshold) {
            ppatch->m_images[0] = visibs[i];
            ppatch->m_images[1] = visibs[j];
            inccThreshold = inccs[i][j];
        }
    }
}
// find the two views with maximum ncc/minimum incc score, these two are definitely the visible views
vector<Vector3f> visibNormals(2);
for (int v = 0; v < 2; ++v) {
    Vector3f icoord = m_pmmvps.m_photoSet.project(visibs[v], ppatch->m_coord, 0);
    int x = (int)(floorf(icoord(0) + 0.5f));
    int y = (int)(floorf(icoord(1) + 0.5f));
    int index = y * m_pmmvps.m_photoSet.getWidth(visibs[v], 0) + x;
    visibNormals[v] = normals[visibs[v]][index];
}

float angleThreshold = cosf(30.0f * M_PI / 180.0f);
if (visibNormals[0].dot(visibNormals[1]) <= angleThreshold) {
    continue;
}

angleThreshold = cosf(10.0f * M_PI / 180.0f);
Vector3f normal3 = (visibNormals[0] + visibNormals[1]);
Vector3f normal3Old = Vector3f::Zero();
int maxIter = 2 * m_pmmvps.m_nimages;
int iter = 0;
while (1) {
    // update patch normal
    normal3 /= static_cast<int>(ppatch->m_images.size());
    normal3 /= normal3.norm();
    ppatch->m_normal << normal3(0), normal3(1), normal3(2), -normal3.dot(ppatch->m_coord.head(3)); // need to test
    if (angleThreshold < normal3.dot(normal3Old) || maxIter <= iter++) {
        break;
    }
    normal3Old = normal3;
    
    // choose visible views based on angle between patch normal and ray of sight
    m_pmmvps.m_optim.addImage(*ppatch);
    // choose visible views based on pairwise incc score
    m_pmmvps.m_optim.constraintImages(*ppatch, m_pmmvps.m_nccThresholdBefore);
    
    normal3 = Vector3f::Zero();
    for (int v = 0; v < static_cast<int>(ppatch->m_images.size()); ++v) {
        Vector3f icoord = m_pmmvps.m_photoSet.project(ppatch->m_images[v], ppatch->m_coord, m_pmmvps.m_level);
        int x = (int)(floorf(icoord(0) + 0.5f));
        int y = (int)(floorf(icoord(1) + 0.5f));
        int index = y * m_pmmvps.m_photoSet.getWidth(visibs[v], m_pmmvps.m_level) + x;
        normal3 += normals[ppatch->m_images[v]][index];
    }
    }
    
    // check if this patch is visible
    if (m_pmmvps.m_optim.preProcess(*ppatch) == -1) {
        continue;
    }
    
    // add the patch to the correponding patch grid
    Vector3f icoord = m_pmmvps.m_photoSet.project(ppatch->m_images[0], ppatch->m_coord, 0);
    int index = icoord(1) * m_pmmvps.m_photoSet.getWidth(ppatch->m_images[0], 0) + icoord(0);
    m_pmmvps.m_patchManager.m_pgrids[ppatch->m_images[0]][index].push_back(ppatch);
*/