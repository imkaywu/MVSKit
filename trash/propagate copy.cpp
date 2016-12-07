//
//  propagate.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "propagate.hpp"
#include "pmmvps.hpp"

using std::cerr;
using std::endl;
using std::flush;
using Eigen::Vector3f;

Propagate::Propagate(PmMvps& pmmvps) : m_pmmvps(pmmvps){
}

Propagate::~Propagate() {
}

void Propagate::init() {
}

void Propagate::run(const int iter) {
    m_pmmvps.m_patchManager.collectPatches(m_queue);
    
    cerr << "Expanding patches..." << flush;
    
    propagate(iter);
}

void Propagate::propagate(const int iter) {
    while (1) {
        pPatch ppatch;
        int isEmpty = 0;
        if (m_queue.empty()) {
            isEmpty = 1;
        }
        else {
            ppatch = m_queue.top();
            m_queue.pop();
        }
        
        if (isEmpty) {
            break;
        }
        // check if this patches has been visited before
        if  (ppatch->m_iter != iter) {
            continue;
        }
//        ++ppatch->m_iter; // add 1 happens after optimization is successful
        
        int image = ppatch->m_images[0];
//        if (ppatch->m_ncc == -1.0f) {
//            m_pmmvps.m_patchManager.computeNcc(*ppatch); // compute the ncc score
//        }
        int gwidth = m_pmmvps.m_patchManager.m_gwidths[m_pmmvps.m_level];
        int gheight = m_pmmvps.m_patchManager.m_gheights[m_pmmvps.m_level];
        int start = 0, end = gwidth * gheight - 1, inc = 1;
        if (iter % 2 == 1) {
            start = gwidth * gheight - 1;
            end = 0;
            inc = -1;
        }
        
        for (int index = start; index != end; index += inc) {
            if (m_pmmvps.m_patchManager.m_pgrids[image][index].empty()) {
                continue;
            }
            
            vector<int> images;
            vector<int> indexes;
            // spatial propagation
            int lr = index + inc;
            int tb = index + gwidth * inc;
            images.push_back(image);    indexes.push_back(lr);
            images.push_back(image);    indexes.push_back(tb);
            
            // view propagation
            findViewNeighbors(ppatch, images, indexes);
            
            // propagate patch to neighboring cells
            propagatePatch(ppatch, images, indexes);
        }
    }
}

void Propagate::findViewNeighbors(pPatch& ppatch, vector<int>& images, vector<int>& indexes) {
    for (int i = 1; i < ppatch->m_nimages; ++i) {
        Vector3f icoord = m_pmmvps.m_photoSet.project(ppatch->m_images[i], ppatch->m_coord, m_pmmvps.m_level);
        int cx = ((int)floorf(icoord(0) + 0.5f)) / m_pmmvps.m_csize; // used in setGrids, setGridsImage
        int cy = ((int)floorf(icoord(1) + 0.5f)) / m_pmmvps.m_csize;
        int gwidth = m_pmmvps.m_patchManager.m_gwidths[ppatch->m_images[i]];
        images.push_back(ppatch->m_images[i]);
        indexes.push_back(cy * gwidth + cx);
    }
}

// need to be cautious that the original patch should not be changed
void Propagate::propagatePatch(pPatch& ppatch, vector<int>& images, vector<int>& indexes) {
    // number of cells to propagate to
    int ncells = static_cast<int>(images.size());
    int isSuccess = 0; // at least one successful propagation
    
    for (int i = 0; i < ncells; ++i) {
//        vector<pPatch> ppatches;
//        if (m_pmmvps.m_patchManager.m_pgrids[images[i]][indexes[i]].empty()) {
//            pPatch ppatchTmp(new Patch(*ppatch)); // ???, we need to create a new object and assign the address to ppatchTmp
//            ppatches.push_back(ppatchTmp); // only 1 ppatch
//        }
//        else {
//            ppatches = m_pmmvps.m_patchManager.m_pgrids[images[i]][indexes[i]]; // could have multiple ppatches
//        }
        
        vector<pPatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[images[i]][indexes[i]];
        int npatches = static_cast<int>(ppatches.size()); // npatches >= 0
        int niter = std::max(1, npatches);
        
        for (int n = 0; n < niter; ++n) {
            if (npatches > 0) {
                if (ppatches[n]->m_ncc == -1.0f) {
                    m_pmmvps.m_patchManager.computeNcc(*ppatches[n]);
                }
                if (ppatch->m_ncc < ppatches[n]->m_ncc) {
                    continue;
                }
                m_pmmvps.m_patchManager.removePatch(ppatches[n]);
            }
            
            Patch patch(*ppatch);
            // patch optimization
            if (m_pmmvps.m_optim.preProcess(patch) == -1) {
                ++m_fcount0;
                return;
            }
            
            m_pmmvps.m_optim.refinePatch(patch, 100);
            
            if (m_pmmvps.m_optim.postProcess(patch) == -1) {
                ++m_fcount1;
                return;
            }
            
            ++m_pcount;
            if (isSuccess == 0) {
                isSuccess = 1;
                ++ppatch->m_iter;
            }
            
            // need to do some testing here
            // if we add this new patch to the m_queue
            int add = 1;
            
            pPatch newppatch(new Patch(patch));
            m_pmmvps.m_patchManager.addPatch(newppatch);
            if (add) {
                m_queue.push(newppatch);
            }
        }
    }
}

//-------------------------------------------------
// trash code
//-------------------------------------------------
/*
void Propagate::propagate() {
    vector<int> images = m_pmmvps.m_images;
    int wsize = m_pmmvps.m_wsize;
    int hwsize = wsize / 2;
    int level = m_pmmvps.m_level;
    
    for (int index = 0; index < images.size(); ++index) {
        int cwidth = m_pmmvps.m_patchManager.m_gwidths[level];
        int cheight = m_pmmvps.m_patchManager.m_gheights[level];
        
        for (int iter = 0; iter < 3; ++iter) {
            int xstart = hwsize - 1, xend = cwidth - hwsize + 1, xinc = 1;
            int ystart = hwsize - 1, yend = cheight - hwsize + 1, yinc = 1;
            if (iter % 2 == 1) {
                xstart = cwidth - hwsize + 1; xend = hwsize - 1; xinc = -1;
                ystart = cheight - hwsize + 1; yend = hwsize - 1; yinc = -1;
            }
            for (int y = ystart; y < yend; y += yinc) {
                for (int x = xstart; x < xend; x += xinc) {
                    int ind = y * cwidth + x;
                    if (m_pmmvps.m_patchManager.m_vgrids[index][ind].empty()) {
                        continue;
                    }
                    
                    vector<pPatch> ppatches = m_pmmvps.m_patchManager.m_vgrids[index][ind];
                    int nPatches = static_cast<int>(ppatches.size());
                    
                    for (int n = 0; n < nPatches; ++n) {
                        pPatch ppatch = ppatches[n];
                        // ---------spatial propagation
                        // propagation along the x-axis
                        ind = y * cwidth + x + xinc;
                        vector<pPatch>& ppatchesNeighbor = m_pmmvps.m_patchManager.m_vgrids[index][ind];
                        int nPatchesNeighbor = static_cast<int>(ppatchesNeighbor.size());
                        if (ppatchesNeighbor.empty()) {
                            ppatchesNeighbor.push_back(ppatch);
                        }
                        else {
                            for (int m = 0; m < nPatchesNeighbor; ++m) {
                                pPatch ppatchNeighbor = ppatchesNeighbor[m];
                                // calculate the photo-consistency score of ppatch
                                float score0 = ppatch->score(m_pmmvps.m_nccThreshold);
                                // calculate the ncc score of ppatchNeighbor
                                float score1 = ppatchNeighbor->score(m_pmmvps.m_nccThreshold);
                                
                                if (score0 < score1) {
                                    ppatchesNeighbor[m] = ppatch;
                                }
//                                vector<vector<Vector3f> > texs;
//                                m_pmmvps.m_photoSet.m_photos[index].getTex(ppatch->m_vgrids[ind], m_pmmvps.m_level, Vector2f(1.0f, 0.0f), Vector2f(0.0f, 1.0f), wsize, texs[0], 1);
//                                m_pmmvps.m_photoSet.m_photos[index].getTex(ppatchNeighbor->m_vgrids[ind], m_pmmvps.m_level, Vector2f(1.0f, 0.0f), Vector2f(0.0f, 1.0f), wsize, texs[1], 1);
                            }
                        }
                        
                        // propagation along the y-axis
                        ind = (y + yinc) & cwidth + x;
                        ppatchesNeighbor = m_pmmvps.m_patchManager.m_vgrids[index][ind];
                        if (ppatchesNeighbor.empty()) {
                            ppatchesNeighbor.push_back(ppatch);
                        }
                        else {
                            for (int m = 0; m < nPatchesNeighbor; ++m) {
                                pPatch ppatchNeighbor = ppatchesNeighbor[m];
                                // calculate the photo-consistency score of ppatch
                                float score0 = ppatch->score(m_pmmvps.m_nccThreshold);
                                // calculate the photo-consistency score of ppatchNeighbor
                                float score1 = ppatchNeighbor->score(m_pmmvps.m_nccThreshold);
                                
                                if (score0 < score1) {
                                    ppatchesNeighbor[m] = ppatch;
                                }
                            }
                        }
                        
                        // ---------view propagation
                        
                    }
                }
            }
        }
    }
}

void Propagate::propagate() {
    while (1) {
        pPatch ppatch;
        int isEmpty = 0;
        if (m_queue.empty()) {
            isEmpty = 1;
        }
        else {
            ppatch = m_queue.top();
            m_queue.pop();
        }
        
        if (isEmpty) {
            break;
        }
        int wsize = m_pmmvps.m_wsize;
        int hwsize = wsize / 2;
        int image = ppatch->m_images[0];
        int cwidth = m_pmmvps.m_patchManager.m_gwidths[image];
        int cheight = m_pmmvps.m_patchManager.m_gheights[image];
        
        for (int iter = 0; iter < 3; ++iter) {
            int xstart = hwsize - 1, xend = cwidth - hwsize + 1, xinc = 1;
            int ystart = hwsize - 1, yend = cheight - hwsize + 1, yinc = 1;
            if (iter % 2 == 1) {
                xstart = cwidth - hwsize + 1; xend = hwsize - 1; xinc = -1;
                ystart = cheight - hwsize + 1; yend = hwsize - 1; yinc = -1;
            }
            for (int y = ystart; y < yend; y += yinc) {
                for (int x = xstart; x < xend; x += xinc) {
                    propagatePatch(image, x, y, xinc, yinc);
                }
            }
        }
    }
}

int Propagate::propagatePatch(const int image, const int x, const int y, const int xinc, const int yinc) {
    int width = m_pmmvps.m_patchManager.m_gwidths[image];
    int index = y * width + x;
    if (m_pmmvps.m_patchManager.m_pgrids[image][index].empty()) {
        cerr << "No patches in the current grid cell" << endl;
        return -1;
    }
    
    vector<pPatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[image][index];
    int npatches = static_cast<int>(ppatches.size());
    
    for (int n = 0; n < npatches; ++n) {
        vector<int> images;
        vector<int> indexes;
        findNeighbors(ppatches[n], image, x, y, xinc, yinc, images, indexes);
        propagatePatch(ppatches[n], images, indexes);
    }
    
    return 0;
}

void Propagate::findNeighbors(pPatch& ppatch, const int image, const int x, const int y, const int xinc, const int yinc, vector<int>& images, vector<int>& indexes) {
    const int width = m_pmmvps.m_patchManager.m_gwidths[image];
    
    // spatial propagation
    images.clear();
    indexes.clear();
    images.push_back(image);
    indexes.push_back(y * width + x + xinc);
    indexes.push_back((y + yinc) * width + x);
    
    // view propagation
    for (int i = 1; i < static_cast<int>(ppatch->m_images.size()); ++i) {
        Vector3f icoord = m_pmmvps.m_photoSet.project(ppatch->m_images[i], ppatch->m_coord, m_pmmvps.m_level);
        // !!!transform to the cell coordinates, need to test
        int cx = (icoord(0) + m_pmmvps.m_csize - 1) / m_pmmvps.m_csize;
        int cy = (icoord(1) + m_pmmvps.m_csize - 1) / m_pmmvps.m_csize;
        images.push_back(ppatch->m_images[i]);
        indexes.push_back(cy * width + cx);
    }
    
}
*/