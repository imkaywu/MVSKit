//
//  propagate.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include <random>
#include "propagate.hpp"
#include "pmmvps.hpp"

using std::cerr;
using std::endl;
using std::flush;

Propagate::Propagate(PmMvps& pmmvps) : m_pmmvps(pmmvps){
}

Propagate::~Propagate() {
}

void Propagate::init() {
    MAX_NUM_OF_PROPAG = 2;
    MAX_NUM_OF_PATCHES = MAX_NUM_OF_PROPAG * m_pmmvps.m_csize * m_pmmvps.m_csize;
}

void Propagate::run(const int iter) {
    m_ecount = 0;
    m_fcount0 = 0;
    m_fcount1 = 0;
    m_pcount = 0;
    
    time_t starttime = time(NULL);
    
    m_pmmvps.m_patchManager.clearCounts();
    m_pmmvps.m_patchManager.clearFlags();
    
    if (!m_queue.empty()) {
        cerr << "queue is not empty in propagate" << endl;
        exit(1);
    }
    m_pmmvps.m_patchManager.collectPatches(m_queue);
    
    cerr << "Expanding patches..." << flush;
    
	// PM-MVS in image space
    propagatePmImage(iter);
	// PM-MVS in scene space
	//propagatePmScene(iter);
	// PMVS
    //propagatePmvs(iter);
    
    cerr << endl
         << "---- EXPANSION: " << (time(NULL) - starttime) << " secs ----" << endl;
    cerr << "total pass fail0 fail1 refinepatch: "
        << m_ecount << " " << m_pcount << " " << m_fcount0 << " " << m_fcount1 << " " << m_pcount + m_fcount1 << endl;
    cerr << "total pass fail0 fail1 refinepatch (%): "
        << 100 * m_ecount / (float)m_ecount << " "
        << 100 * m_pcount / (float)m_ecount << " "
        << 100 * m_fcount0 / (float)m_ecount << " "
        << 100 * m_fcount1 / (float)m_ecount << " "
        << 100 * (m_pcount + m_fcount1) / (float)m_ecount << endl;
}

//-------------------------------------------
// implementation of PM-MVS in image space
//-------------------------------------------

// for each patch, we propagate MAX_NUM_OF_PROPAG to the spatial/view neighbours if the cell is not full.
// if the cell is full, we refine the worst MAX_NUM_OF_PROPAG patches
void Propagate::propagatePmImage(const int iter) {
    for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
		cerr << "---------------------" << endl
			 << "propagate on image " << image << endl
			 << "---------------------" << endl;

        int gwidth = m_pmmvps.m_patchManager.m_gwidths[image];
        int gheight = m_pmmvps.m_patchManager.m_gheights[image];
        int start = 0, end = gwidth * gheight - 1, inc = 1;
        if (iter % 2 == 1) {
            start = gwidth * gheight - 1;
            end = 0;
            inc = -1;
        }

        for (int index = start; index != end; index += inc) {
            vector<Ppatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[image][index];
            m_pmmvps.m_patchManager.sortPatches(ppatches, 0); // descending order
            int npatches = static_cast<int>(ppatches.size());
            if (npatches == 0) {
                continue;
            }
            else if (npatches > MAX_NUM_OF_PATCHES) {
                for (int i = npatches - 1; i >= MAX_NUM_OF_PATCHES; --i) {
                    m_pmmvps.m_patchManager.removePatch(ppatches[i]);
                }
                npatches = MAX_NUM_OF_PATCHES;
            }
            cerr << "# source patches: " << npatches << endl;
            
            for (int n = 0; n < npatches; ++n) {
                // spatial propagation
                if (ppatches[n]->m_images[0] == image) {
                    // propagate left-right
                    propagatePatch(ppatches[n], image, index + inc);
                    // propagate up-down
                    propagatePatch(ppatches[n], image, index + inc * gwidth);
                }
                // view propagation, temporarily commented
                //else {
                //    const float fx = ppatches[n]->m_grids[0](0);
                //    const float fy = ppatches[n]->m_grids[0](1);
                //    const int cx = (int)(floorf(fx + 0.5f)) / m_pmmvps.m_csize;
                //    const int cy = (int)(floorf(fy + 0.5f)) / m_pmmvps.m_csize;
                //    // left-right
                //    propagatePatch(ppatches[n], ppatches[n]->m_images[0], cy * gwidth + cx);
                //    // up-down
                //    propagatePatch(ppatches[n], ppatches[n]->m_images[0], (cy + inc) * gwidth + cx);
                //}
            }
        }
    }
}

int Propagate::propagatePatch(const Ppatch& ppatch, const int image, const int index) {
    vector<Ppatch>& ppatches = m_pmmvps.m_patchManager.m_pgrids[image][index];
    m_pmmvps.m_patchManager.sortPatches(ppatches, 0); // descending order
    int npatches = static_cast<int>(ppatches.size());
    if (npatches > MAX_NUM_OF_PATCHES) {
        for (int i = npatches - 1; i >= MAX_NUM_OF_PATCHES; --i) {
            m_pmmvps.m_patchManager.removePatch(ppatches[i]);
        }
        npatches = MAX_NUM_OF_PATCHES;
    }
    cerr << "# dest patches: " << npatches << endl;
    
    // random number generator
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(-0.5, 0.5); // [-1, 1)
    
    // the center coordinate of the cell
    const int gwidth = m_pmmvps.m_patchManager.m_gwidths[image];
    const int cx = index % gwidth;
    const int cy = index / gwidth;
    Vector3f icoord;
    icoord << (m_pmmvps.m_csize * (2 * cx + 1) - 1) / 2.0f,
              (m_pmmvps.m_csize * (2 * cy + 1) - 1) / 2.0f,
              1.0f;
    
    // propagation and refinement
    int pcount = 0, rcount = 0;
    for (int iter = 0; iter < MAX_NUM_OF_PROPAG; ++iter) {
        npatches = static_cast<int>(ppatches.size());
        Ppatch newppatch = NULL;
        if (npatches < MAX_NUM_OF_PATCHES) {
            // determine the increment of the icoord;
            Vector3f newicoord;
            newicoord << distribution(generator) * m_pmmvps.m_csize, distribution(generator) * m_pmmvps.m_csize, 0.0f;
            newicoord = icoord + newicoord;
            newppatch = generatePatch(ppatch, newicoord);
			if (!newppatch) { // newppatch == NULL
				continue;
			}
        }
        else {
            m_pmmvps.m_patchManager.sortPatches(ppatches, 0); // descending order;
            Vector3f icoord = m_pmmvps.m_photoSet.project(image, ppatches[MAX_NUM_OF_PATCHES - 1]->m_coord, m_pmmvps.m_level);
            newppatch = generatePatch(ppatch, icoord);
			if (!newppatch || newppatch->m_ncc < ppatches[MAX_NUM_OF_PATCHES - 1]->m_ncc) { // newppatch != NULL
				continue;
			}
        }

		//const int flag = checkCounts(*newppatch);
		//if (flag) {
		//	return -1;
		//}
		//++m_ecount;
        
        //----------patch optimization
        if (m_pmmvps.m_optim.preProcess(*newppatch) == -1) {
            ++m_fcount0;
            //cerr << "failed before optimization" << endl;
            continue;
        }
        
        //cerr << "before refinement: " << newppatch->m_coord << endl << newppatch->m_normal << endl;
        m_pmmvps.m_optim.refinePatch(*newppatch, 100);
        //cerr << "after refinement: " << newppatch->m_coord << endl << newppatch->m_normal << endl;
        
        if (m_pmmvps.m_optim.postProcess(*newppatch) == -1) {
            ++m_fcount1;
            //cerr << "failed after optimization" << endl;
            continue;
        }
        
        if (npatches == MAX_NUM_OF_PATCHES) {
            m_pmmvps.m_patchManager.removePatch(ppatches[MAX_NUM_OF_PATCHES - 1]);
            ++rcount;
        }
		else {
			++pcount;
		}

		//const int add = updateCounts(*newppatch);

		m_pmmvps.m_patchManager.addPatch(newppatch);

		//if (add) {
		//	m_queue.push(ppatch);
		//}
    }
    cerr << "# propag: " << pcount << endl;
    cerr << "# refine: " << rcount << endl;
    
    return 0;
}

Ppatch Propagate::generatePatch(const Ppatch& ppatch, const Vector3f& icoord) const {
    Patch newpatch;
    // update the coord of the new patch
	const int image = ppatch->m_images[0];
    const float depth = m_pmmvps.m_photoSet.m_photos[image].m_oaxis.dot(ppatch->m_coord);
    Vector3f newicoord = depth * icoord;
    newpatch.m_coord = m_pmmvps.m_photoSet.m_photos[image].unproject(newicoord, m_pmmvps.m_level);
	// update the normal of the new patch
	newpatch.m_normal = ppatch->m_normal;
    // set grids and images [???the grid of the reference SHOULD be fixed since only the depth is updated]
    m_pmmvps.m_patchManager.setGridsImages(newpatch, ppatch->m_images);
	if (newpatch.m_images.empty()) {
		return NULL;
	}
    // compute the ncc score
    m_pmmvps.m_patchManager.computeNcc(newpatch);
    return Ppatch(new Patch(newpatch));
}

//-------------------------------------------
// implementation of PM-MVS in scene space
//-------------------------------------------

void Propagate::propagatePmScene(const int iter) {
	while (1) {
		Ppatch ppatch;
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

		// 
	}
}

//-------------------------------------------
// implementation of PM-MVS (old code)
//-------------------------------------------

// need to be cautious that the original patch should not be changed
// use the view we're propagating to as the reference view
int Propagate::propagatePatch(Ppatch& ppatch, vector<int>& indexes) {
    int refImage = ppatch->m_images[0];
    int ncells = static_cast<int>(indexes.size()); // number of neighboring cells to propagate to
    int isSuccess = 0; // at least one successful propagation
    
    for (int i = 0; i < ncells; ++i) {
        vector<Ppatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[refImage][indexes[i]]; // patches in the neighboring cell
        m_pmmvps.m_patchManager.sortPatches(ppatches, 0); // descending order
        int npatches = std::min(static_cast<int>(ppatches.size()), MAX_NUM_OF_PATCHES); // number of patches to get refined in the neighboring cell
        int niter = std::max(1, npatches);
        
        vector<Ppatch>& ppatchesRef = m_pmmvps.m_patchManager.m_pgrids[refImage][indexes[i]]; // used to remove patches from grid
        
        //----------initialize a new patch for each cell
        Patch patch(*ppatch);
        // estimate depth of the current patch
        const float depth = m_pmmvps.m_photoSet.m_photos[refImage].m_oaxis.dot(ppatch->m_coord);
        // backproject to get the coord for the neighbouring pixel
        const int cx = indexes[i] % m_pmmvps.m_patchManager.m_gwidths[refImage];
        const int cy = indexes[i] / m_pmmvps.m_patchManager.m_gwidths[refImage];
        const float fx = (m_pmmvps.m_csize * (2 * cx + 1) - 1) / 2.0f;
        const float fy = (m_pmmvps.m_csize * (2 * cy + 1) - 1) / 2.0f;
        Vector3f icoord;
        icoord << std::max(0.0f, std::min((float)m_pmmvps.m_photoSet.getWidth(refImage, m_pmmvps.m_level), fx)),
                  std::max(0.0f, std::min((float)m_pmmvps.m_photoSet.getWidth(refImage, m_pmmvps.m_level), fy)),
                  1.0f; // probably won't be outside of the image, but just in case
        icoord = depth * icoord;
        patch.m_coord = m_pmmvps.m_photoSet.m_photos[refImage].unproject(icoord, m_pmmvps.m_level);
        // set grids, the grid of the reference SHOULD be fixed since only the depth is updated
        m_pmmvps.m_patchManager.setGrids(patch);
        // recompute the ncc score
        m_pmmvps.m_patchManager.computeNcc(patch);
        // check if the grid coord to image coord is correct
        if (patch.m_grids[0](0) != cx || patch.m_grids[0](1) != cy) {
            cerr << "!!!pixel coordinate wrong" << endl;
        }
        
        for (int n = 0; n < niter; ++n) {
            cerr << "target: cell index: " << indexes[i] << ", total patches: " << std::max(1, npatches) << ", patch index: " << n + 1 << endl;
            // a NON-empty cell
            if (npatches > 0) {
                // they're the SAME patch, and the ncc score is not higher
                if (ppatch == ppatches[n] || ppatch->m_ncc <= ppatches[n]->m_ncc) {
                    cerr << "patch ncc is too lower to get propagated any more" << endl;
                    return -1;
                }
            }
            
            // generate a new patch for each patch in this cell
            Patch newpatch(patch);
            
            //----------patch optimization
            if (m_pmmvps.m_optim.preProcess(newpatch) == -1) {
                ++m_fcount0;
                cerr << "failed before optimization" << endl;
                break; // this patch won't work in this cell
            }
            
            m_pmmvps.m_optim.refinePatch(newpatch, 100);
            
            // the optimized patch is lower than the original one
            if (npatches > 0 && newpatch.m_ncc <= ppatches[n]->m_ncc) {
                cerr << "the optimized patch has lower ncc score" << endl;
                continue;
            }
            
            if (m_pmmvps.m_optim.postProcess(newpatch) == -1) {
                ++m_fcount1;
                cerr << "failed after optimization" << endl;
                break; // this patch won't work in this cell
            }
            
            ++m_pcount;
            isSuccess = 1;
            
            // patches that is currently being refined need to be removed from the grid
//            cerr << "# patches before removing: " << npatches << endl;
            ppatchesRef.erase(remove(ppatchesRef.begin(), ppatchesRef.end(), ppatches[n]), ppatchesRef.end()); // ppatches is ordered
//            cerr << "# patches after removing: " << static_cast<int>(ppatchesRef.size()) << endl;
//            cerr << "# patches after removing (double checking): " << static_cast<int>(m_pmmvps.m_patchManager.m_pgrids[refImage][indexes[i]].size()) << endl;

            Ppatch newppatch(new Patch(newpatch));
            
            m_pmmvps.m_patchManager.addPatch(newppatch);

        }
        // to make sure each cell has most 2 * MAX_NUM_OF_PATCHES patches
        m_pmmvps.m_patchManager.sortPatches(ppatchesRef, 0); // descending order
        for (int n = static_cast<int>(ppatchesRef.size()) - 1; n >= 2 * MAX_NUM_OF_PATCHES; --n) {
            Ppatch ppatch = ppatches[n];
            m_pmmvps.m_patchManager.removePatch(ppatch);
        }
    }
    if (isSuccess == 0) {
        return 0;
    }
//    ++ppatch->m_iter;
    return -1; // current patch successfully propagated, the lower ncc patch is less likely to succeed
}

// not used
void Propagate::findViewNeighbors(Ppatch& ppatch, vector<int>& images, vector<int>& indexes) {
    for (int i = 1; i < ppatch->m_nimages; ++i) {
        Vector3f icoord = m_pmmvps.m_photoSet.project(ppatch->m_images[i], ppatch->m_coord, m_pmmvps.m_level);
        int cx = ((int)floorf(icoord(0) + 0.5f)) / m_pmmvps.m_csize; // used in setGrids, setGridsImages
        int cy = ((int)floorf(icoord(1) + 0.5f)) / m_pmmvps.m_csize;
        int gwidth = m_pmmvps.m_patchManager.m_gwidths[ppatch->m_images[i]];
        images.push_back(ppatch->m_images[i]);
        indexes.push_back(cy * gwidth + cx);
    }
}

//-------------------------------------------
// implementation of PMVS
//-------------------------------------------
void Propagate::propagatePmvs(const int iter) {
	while (1) {
		Ppatch ppatch;
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

		vector<vector<Vector4f, Eigen::aligned_allocator<Vector4f> > > canCoords;
		findEmptyBlocks(ppatch, canCoords);

		for (int i = 0; i < static_cast<int>(canCoords.size()); ++i) {
			for (int j = 0; j < static_cast<int>(canCoords[i].size()); ++j) {
				const int flag = expandSub(ppatch, canCoords[i][j]);
				if (flag == -1) {
					ppatch->m_dflag |= (0x0001) << i;
				}
			}
		}
	}
}

void Propagate::findEmptyBlocks(const Ppatch& ppatch, vector<vector<Vector4f, Eigen::aligned_allocator<Vector4f> > >& canCoords){
    const int dnum = 6; // dnum must be at most 8, because m_dflag is char???
    const Patch& patch = *ppatch;
    
    // six empty directions
    Vector4f xdir, ydir;
    ortho(ppatch->m_normal, xdir, ydir);
    
    vector<float> fill(dnum);
    std::fill(fill.begin(), fill.end(), 0.0f);
    
    const float radius = computeRadius(patch);
    const float radiuslow = radius / 6.0f;
    const float radiushigh = radius * 2.5f;
    
    vector<Ppatch> neighbors;
    m_pmmvps.m_patchManager.findNeighbors(patch, neighbors, 4.0, 1);
	//cerr << "# of neighbors: " << neighbors.size() << endl;
    
    vector<Ppatch>::iterator bpatch = neighbors.begin();
    vector<Ppatch>::iterator epatch = neighbors.end();
    while (bpatch != epatch) {
        const Vector4f diff = (*bpatch)->m_coord - ppatch->m_coord;
        Eigen::Vector2f f2;
        f2 << diff.dot(xdir), diff.dot(ydir);
        const float len = f2.norm();
        if (len < radiuslow || radiushigh < len) {
            ++bpatch;
            continue;
        }
        f2 = f2 / len;
        
        float angle = atan2f(f2(1), f2(0));
        if (angle < 0.0f) {
            angle += 2 * M_PI;
        }
        
        const float findex = angle / (2 * M_PI / dnum);
        const int lindex = (int)floorf(findex);
        const int hindex = lindex + 1;
        
        fill[lindex % dnum] += hindex - findex;
        fill[hindex % dnum] += findex - lindex;
        ++bpatch;
    }
    
    canCoords.resize(dnum);
    for (int i = 0; i < dnum; ++i) {
        if (0.0f < fill[i]) {
            continue;
        }
        if (ppatch->m_dflag & (0x0001 << i)) {
            continue;
        }
        const float angle = 2 * M_PI * i / dnum;
        Vector4f canCoord = ppatch->m_coord + cosf(angle) * radius * xdir + sinf(angle) * radius * ydir;
        canCoords[i].push_back(canCoord);
    }
}

float Propagate::computeRadius(const Patch &patch) {
    const int minnum = 2;
    vector<float> units;
    m_pmmvps.m_optim.computeUnits(patch, units);
    vector<float> vftmp = units;
    nth_element(vftmp.begin(), vftmp.begin() + minnum - 1, vftmp.end()); // why select the second smallest number as the raidus
    return (*(vftmp.begin() + minnum - 1)) * m_pmmvps.m_csize;
}

void Propagate::ortho(const Vector4f& z, Vector4f& x, Vector4f& y) {
    if (fabsf(z(0)) > 0.5f) {
        x << z(1), -z(0), 0, 0;
    }
    else if(fabsf(z(1)) > 0.5f) {
        x << 0, z(2), -z(1), 0;
    }
    else {
        x << -z(2), 0, z(0), 0;
    }
    x /= x.norm();
    y << z(1) * x(2) - z(2) * x(1),
         z(2) * x(0) - z(0) * x(2),
         z(0) * x(1) - z(1) * x(0),
         0;
}

int Propagate::expandSub(const Ppatch& orgppatch, const Vector4f& canCoord) {
    // Choose the closest one
    Patch patch;
    patch.m_coord = canCoord;
    patch.m_normal = orgppatch->m_normal;
    patch.m_flag = 1;
    
    m_pmmvps.m_patchManager.setGridsImages(patch, orgppatch->m_images);
    if (patch.m_images.empty())
        return -1;
    
    //-----------------------------------------------------------------
    // Check bimages and mask. Then, initialize possible visible images
    // if mask doesn't exist, which returns -1, or
    // if mask return 255 if the patch does exist
    if (m_pmmvps.m_photoSet.getMask(patch.m_coord, m_pmmvps.m_level) == 0) {
        return -1;
    }
    
    // Check m_counts and maybe m_pgrids
    const int flag = checkCounts(patch);
    if (flag)
        return -1;
    
    ++m_ecount;
    //-----------------------------------------------------------------
    // Preprocess
    if (m_pmmvps.m_optim.preProcess(patch) == -1) {
        ++m_fcount0;
        return -1;
    }
    
    //-----------------------------------------------------------------
    m_pmmvps.m_optim.refinePatch(patch, 100);
    
    //-----------------------------------------------------------------
    if (m_pmmvps.m_optim.postProcess(patch) == -1) {
        ++m_fcount1;
        return -1;
    }
    ++m_pcount;
    
    //-----------------------------------------------------------------
    // Finally
    Ppatch ppatch(new Patch(patch));
    
    //patch.m_images = orgppatch->m_images;
    const int add = updateCounts(patch);
    
    m_pmmvps.m_patchManager.addPatch(ppatch);
    
    if (add) {
        m_queue.push(ppatch);
    }    
    
    return 0;
}

int Propagate::checkCounts(Patch& patch) {
    int full = 0;  int empty = 0;
    
    vector<int>::iterator begin = patch.m_images.begin();
    vector<int>::iterator end = patch.m_images.end();
    vector<Vector2i>::iterator begin2 = patch.m_grids.begin();
    
    while (begin != end) {
        const int index = *begin;
        
        const int ix = (*begin2)(0);
        const int iy = (*begin2)(1);
        if (ix < 0 || m_pmmvps.m_patchManager.m_gwidths[index] <= ix ||
            iy < 0 || m_pmmvps.m_patchManager.m_gheights[index] <= iy) {
            ++begin;
            ++begin2;
            continue;
        }
        
        const int index2 = iy * m_pmmvps.m_patchManager.m_gwidths[index] + ix;
        
        int flag = 0;
        if (!m_pmmvps.m_patchManager.m_pgrids[index][index2].empty()) {
            flag = 1;
        }
        if (flag) {
            ++full;
            ++begin;
            ++begin2;
            continue;
        }
        
        if (m_pmmvps.m_countThreshold1 <= m_pmmvps.m_patchManager.m_counts[index][index2]) {
            ++full;
        }
        else {
            ++empty;
        }
        ++begin;
        ++begin2;
    }
    
    //First expansion is expensive and make the condition strict
    if (m_pmmvps.m_depth <= 1) {
        if (empty < m_pmmvps.m_minImageNumThreshold && full != 0) {
            return 1;
        }
        else {
            return 0;
        }
    }
    else {
        if (empty < m_pmmvps.m_minImageNumThreshold - 1 && full != 0) {
            return 1;
        }
        else {
            return 0;
        }
    }
}

int Propagate::updateCounts(const Patch& patch) {
    int full = 0;  int empty = 0;
    
    {
        vector<int>::const_iterator begin = patch.m_images.begin();
        vector<int>::const_iterator end = patch.m_images.end();
        vector<Vector2i>::const_iterator begin2 = patch.m_grids.begin();
        
        while (begin != end) {
            const int index = *begin;
            
            const int ix = (*begin2)(0);
            const int iy = (*begin2)(1);
            if (ix < 0 || m_pmmvps.m_patchManager.m_gwidths[index] <= ix ||
                iy < 0 || m_pmmvps.m_patchManager.m_gheights[index] <= iy) {
                ++begin;
                ++begin2;
                continue;
            }
            
            const int index2 = iy * m_pmmvps.m_patchManager.m_gwidths[index] + ix;
            
            if (m_pmmvps.m_countThreshold1 <= m_pmmvps.m_patchManager.m_counts[index][index2]) {
                ++full;
            }
            else {
                ++empty;
            }
            
            ++m_pmmvps.m_patchManager.m_counts[index][index2];
            ++begin;
            ++begin2;
        }
    }
    
    {
        vector<int>::const_iterator begin = patch.m_vimages.begin();
        vector<int>::const_iterator end = patch.m_vimages.end();
        vector<Vector2i>::const_iterator begin2 = patch.m_vgrids.begin();
        
        while (begin != end) {
            const int index = *begin;
            
            const int ix = (*begin2)(0);
            const int iy = (*begin2)(1);
            if (ix < 0 || m_pmmvps.m_patchManager.m_gwidths[index] <= ix ||
                iy < 0 || m_pmmvps.m_patchManager.m_gheights[index] <= iy) {
                ++begin;
                ++begin2;
                continue;
            }
            
            const int index2 = iy * m_pmmvps.m_patchManager.m_gwidths[index] + ix;
            
            if (m_pmmvps.m_countThreshold1 <= m_pmmvps.m_patchManager.m_counts[index][index2]) {
                ++full;
            }
            else {
                ++empty;
            }
            ++m_pmmvps.m_patchManager.m_counts[index][index2];
            ++begin;
            ++begin2;
        }
    }
    
    if (empty != 0) {
        return 1;
    }
    else {
        return 0;
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
                    if (m_pmmvps.m_patchManager.m_vpgrids[index][ind].empty()) {
                        continue;
                    }
                    
                    vector<Ppatch> ppatches = m_pmmvps.m_patchManager.m_vpgrids[index][ind];
                    int nPatches = static_cast<int>(ppatches.size());
                    
                    for (int n = 0; n < nPatches; ++n) {
                        Ppatch ppatch = ppatches[n];
                        // ---------spatial propagation
                        // propagation along the x-axis
                        ind = y * cwidth + x + xinc;
                        vector<Ppatch>& ppatchesNeighbor = m_pmmvps.m_patchManager.m_vpgrids[index][ind];
                        int nPatchesNeighbor = static_cast<int>(ppatchesNeighbor.size());
                        if (ppatchesNeighbor.empty()) {
                            ppatchesNeighbor.push_back(ppatch);
                        }
                        else {
                            for (int m = 0; m < nPatchesNeighbor; ++m) {
                                Ppatch ppatchNeighbor = ppatchesNeighbor[m];
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
                        ppatchesNeighbor = m_pmmvps.m_patchManager.m_vpgrids[index][ind];
                        if (ppatchesNeighbor.empty()) {
                            ppatchesNeighbor.push_back(ppatch);
                        }
                        else {
                            for (int m = 0; m < nPatchesNeighbor; ++m) {
                                Ppatch ppatchNeighbor = ppatchesNeighbor[m];
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
        Ppatch ppatch;
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
    
    vector<Ppatch> ppatches = m_pmmvps.m_patchManager.m_pgrids[image][index];
    int npatches = static_cast<int>(ppatches.size());
    
    for (int n = 0; n < npatches; ++n) {
        vector<int> images;
        vector<int> indexes;
        findNeighbors(ppatches[n], image, x, y, xinc, yinc, images, indexes);
        propagatePatch(ppatches[n], images, indexes);
    }
    
    return 0;
}

void Propagate::findNeighbors(Ppatch& ppatch, const int image, const int x, const int y, const int xinc, const int yinc, vector<int>& images, vector<int>& indexes) {
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