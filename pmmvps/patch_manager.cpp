//
//  patch_manager.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/19/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include <fstream>
#include "patch_manager.hpp"
#include "pmmvps.hpp"

using std::ifstream;
using std::ofstream;
using Eigen::Vector2i;
using Eigen::Vector3f;

Ppatch PatchManager::m_MAXDEPTH(new Patch());
Ppatch PatchManager::m_BACKGROUND(new Patch());

PatchManager::PatchManager(PmMvps& pmmvps) : m_pmmvps(pmmvps) {
}

void PatchManager::init()
{
    m_nimages = m_pmmvps.m_nimages;
    m_pgrids.clear();      m_pgrids.resize(m_nimages);
    m_vpgrids.clear();     m_vpgrids.resize(m_nimages);
    m_dpgrids.clear();     m_dpgrids.resize(m_nimages);
    m_gheights.clear();    m_gheights.resize(m_nimages);
    m_gwidths.clear();     m_gwidths.resize(m_nimages);
    m_counts.clear();      m_counts.resize(m_nimages);
    
    for (int index = 0; index < m_nimages; ++index)
    {
        const int gheight = (m_pmmvps.m_photoSet.getHeight(index, m_pmmvps.m_level) + m_pmmvps.m_csize - 1) / m_pmmvps.m_csize;
        const int gwidth = (m_pmmvps.m_photoSet.getWidth(index, m_pmmvps.m_level) + m_pmmvps.m_csize - 1) / m_pmmvps.m_csize;
        const int ngrids = gheight * gwidth;
        
        m_gheights[index] = gheight;
        m_gwidths[index] = gwidth;
        m_pgrids[index].resize(ngrids);
        m_vpgrids[index].resize(ngrids);
        m_dpgrids[index].resize(ngrids);
        m_counts[index].resize(ngrids);
        
        for (int i = 0; i < ngrids; ++i) {
            m_dpgrids[index][i] = m_MAXDEPTH;
        }
    }
}

void PatchManager::image2index(Patch& patch) {
    vector<int> newimages;
    for (int i = 0; i < static_cast<int>(patch.m_images.size()); ++i) {
        const int index = m_pmmvps.m_photoSet.image2index(patch.m_images[i]);
        if (index != -1) {
            newimages.push_back(index);
        }
    }
    
    patch.m_images.swap(newimages);
}

void PatchManager::index2image(Patch& patch) {
    for (int i = 0; i < static_cast<int>(patch.m_images.size()); ++i) {
        patch.m_images[i] = m_pmmvps.m_photoSet.m_images[patch.m_images[i]];
    }
    
    for (int i = 0; i < static_cast<int>(patch.m_vimages.size()); ++i) {
        patch.m_vimages[i] = m_pmmvps.m_photoSet.m_images[patch.m_vimages[i]];
    }
}

void PatchManager::collectPatches(const int target) {
    m_ppatches.clear();
    
    for (int index = 0; index < m_pmmvps.m_nimages; ++index) {
        for (int i = 0; i < static_cast<int>(m_pgrids[index].size()); ++i) {
            vector<Ppatch>::iterator bpatch = m_pgrids[index][i].begin();
            vector<Ppatch>::iterator epatch = m_pgrids[index][i].end();
            while (bpatch != epatch) {
                (*bpatch)->m_id = -1;
                ++bpatch;
            }
        }
    }
    
    int count = 0;
    for (int index = 0; index < m_pmmvps.m_nimages; ++index) {
        for (int i = 0; i < static_cast<int>(m_pgrids[index].size()); ++i) {
            vector<Ppatch>::iterator bpatch = m_pgrids[index][i].begin();
            vector<Ppatch>::iterator epatch = m_pgrids[index][i].end();
            while (bpatch != epatch) {
                if ((*bpatch)->m_id == -1) {
                    (*bpatch)->m_id = count++;
                    
                    if (target == 0 || (*bpatch)->m_fix == 0)
                        m_ppatches.push_back(*bpatch);
                }
                ++bpatch;
            }
        }
    }
}

void PatchManager::collectPatches(priority_queue<Ppatch, vector<Ppatch>, PatchCmp>& pqPatches) {
    for (int index = 0; index < m_pmmvps.m_nimages; ++index) {
        for (int i = 0; i < static_cast<int>(m_pgrids[index].size()); ++i) {
            vector<Ppatch>::iterator bpatch = m_pgrids[index][i].begin();
            vector<Ppatch>::iterator epatch = m_pgrids[index][i].end();
            while (bpatch != epatch) {
                if ((*bpatch)->m_flag == 0) {
                    (*bpatch)->m_flag = 1;
                    pqPatches.push(*bpatch);
                }
                ++bpatch;
            }
        }
    }
}

void PatchManager::collectPatches(const int index, priority_queue<Ppatch, vector<Ppatch>, PatchCmp>& pqPatches) {
    for (int i = 0; i < static_cast<int>(m_pgrids[index].size()); ++i) {
        vector<Ppatch>::iterator bpatch = m_pgrids[index][i].begin();
        vector<Ppatch>::iterator epatch = m_pgrids[index][i].end();
        
        while (bpatch != epatch) {
            if ((*bpatch)->m_images[0] == index && (*bpatch)->m_flag == 0) {
                (*bpatch)->m_flag = 1;
                pqPatches.push(*bpatch);
            }
            ++bpatch;
        }
    }
}

void PatchManager::clearCounts() {
    for (int index = 0; index < m_pmmvps.m_nimages; ++index) {
        vector<unsigned char>::iterator begin = m_counts[index].begin();
        vector<unsigned char>::iterator end = m_counts[index].end();
        while (begin != end) {
            *begin = (unsigned char)0;
            ++begin;
        }
    }
}

void PatchManager::clearFlags() {
    vector<Ppatch>::iterator bppatch = m_ppatches.begin();
    vector<Ppatch>::iterator eppatch = m_ppatches.end();
    while (bppatch != eppatch) {
        (*bppatch)->m_flag = 0;
        ++bppatch;
    }
}

void PatchManager::addPatch(Ppatch& ppatch)
{
    vector<int>::iterator bimage = ppatch->m_images.begin();
    vector<int>::iterator eimage = ppatch->m_images.end();
    vector<Vector2i>::iterator bgrid = ppatch->m_grids.begin();
    
    while (bimage != eimage) {
        const int image = *bimage;
        const int index = (*bgrid)(1) * m_gwidths[image] + (*bgrid)(0);
        m_pgrids[image][index].push_back(ppatch);
        ++bimage;
        ++bgrid;
    }
    
    if (m_pmmvps.m_depth == 0) {
        return;
    }
    
    bimage = ppatch->m_vimages.begin();
    eimage = ppatch->m_vimages.end();
    bgrid = ppatch->m_vgrids.begin();
    while(bimage != eimage)
    {
        const int image = *bimage;
        const int index = (*bgrid)(1) * m_gwidths[image] + (*bgrid)(0); // y * width + x;
        m_vpgrids[image][index].push_back(ppatch);
        ++bimage;
        ++bgrid;
    }
    
    updateDepthMaps(ppatch);
}

void PatchManager::updateDepthMaps(Ppatch &ppatch)
{
    for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
        const Vector3f icoord = m_pmmvps.m_photoSet.project(image, ppatch->m_coord, m_pmmvps.m_level);
        
        const float fx = icoord(0) / m_pmmvps.m_csize;
        const int x[2] = {(int)floor(fx), (int)ceil(fx)};
        const float fy = icoord(1) / m_pmmvps.m_csize;
        const int y[2] = {(int)floor(fy), (int)ceil(fy)};
        
        const float depth = m_pmmvps.m_photoSet.m_photos[image].m_oaxis.dot(ppatch->m_coord);
        
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                if(x[i] < 0 || m_gwidths[image] <= x[i] || y[j] < 0 || m_gheights[image] <= y[j]){
                    continue;
                }
                const int index = y[j] * m_gwidths[image] + x[i];
                if(m_dpgrids[image][index] == m_MAXDEPTH) {
                    m_dpgrids[image][index] = ppatch;
                }
                else {
                    const float d = m_pmmvps.m_photoSet.m_photos[image].m_oaxis.dot(m_dpgrids[image][index]->m_coord);
                    if (depth < d) {
                        m_dpgrids[image][index] = ppatch;
                    }
                }
            }
        }
    }
}

void PatchManager::setGridsImages(Patch& patch, vector<int>& images) const {
    patch.m_images.clear();
    patch.m_grids.clear();
    vector<int>::const_iterator bimage = images.begin();
    vector<int>::const_iterator eimage = images.end();
    while (bimage != eimage) {
        const Vector3f icoord = m_pmmvps.m_photoSet.project(*bimage, patch.m_coord, m_pmmvps.m_level);
        const int ix = ((int)floorf(icoord(0) + 0.5f)) / m_pmmvps.m_csize;
        const int iy = ((int)floorf(icoord(1) + 0.5f)) / m_pmmvps.m_csize;
        if (0 <= ix && ix < m_gwidths[*bimage] &&
            0 <= iy && iy < m_gheights[*bimage]) {
            patch.m_images.push_back(*bimage);
            patch.m_grids.push_back(Vector2i(ix, iy));
        }
        ++bimage;
    }
}

void PatchManager::setGrids(Patch& patch) {
    patch.m_grids.clear();
    for (int i = 0; i < static_cast<int>(patch.m_images.size()); ++i) {
        const int image = patch.m_images[i]; // those are index, image2index() called
        Vector3f icoord = m_pmmvps.m_photoSet.project(image, patch.m_coord, m_pmmvps.m_level);
        const int ix = ((int)floorf(icoord(0) + 0.5f)) / m_pmmvps.m_csize;
        const int iy = ((int)floorf(icoord(1) + 0.5f)) / m_pmmvps.m_csize;
        patch.m_grids.push_back(Vector2i(ix, iy));
    }
}

void PatchManager::setVGrids(Patch& patch) {
	patch.m_vgrids.clear();
	for (int i = 0; i < static_cast<int>(patch.m_vimages.size()); i++) {
		const int image = patch.m_vimages[i];
		Vector3f icoord = m_pmmvps.m_photoSet.project(image, patch.m_coord, m_pmmvps.m_level);
		const int ix = ((int)floorf(icoord(0) + 0.5f)) / m_pmmvps.m_csize;
		const int iy = ((int)floorf(icoord(1) + 0.5f)) / m_pmmvps.m_csize;
		patch.m_vgrids.push_back(Vector2i(ix, iy));
	}
}

void PatchManager::setVImagesVGrids(Ppatch& ppatch) {
    setVImagesVGrids(*ppatch);
}

void PatchManager::setVImagesVGrids(Patch& patch) {
    vector<int> visib;
    visib.resize(m_pmmvps.m_nimages);
    fill(visib.begin(), visib.end(), 0);
    
    // collect all the visible images from m_images and m_vimages
    vector<int>::iterator bimage = patch.m_images.begin();
    vector<int>::iterator eimage = patch.m_images.end();
    while (bimage != eimage) {
        visib[*bimage] = 1;
        ++bimage;
    }
    
    bimage = patch.m_vimages.begin();
    eimage = patch.m_vimages.end();
    while (bimage != eimage) {
        visib[*bimage] = 1;
        ++bimage;
    }
    
    for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
        if (visib[image]) {
            continue;
        }
        int ix, iy;
        if (isVisible0(patch, image, ix, iy, m_pmmvps.m_neighborThreshold) == 0) {
            continue;
        }
        
        // check edge, see the original source code, not implemented yet
        
        patch.m_vimages.push_back(image);
        patch.m_vgrids.push_back(Vector2i(ix, iy));
    }
}

void PatchManager::removePatch(const Ppatch& ppatch) {
    for (int i = 0; i < static_cast<int>(ppatch->m_images.size()); ++i) {
        const int image = ppatch->m_images[i];
        const int& ix = ppatch->m_grids[i](0);
        const int& iy = ppatch->m_grids[i](1);
        const int index = iy * m_gwidths[image] + ix;
        m_pgrids[image][index].erase(remove(m_pgrids[image][index].begin(),
                                            m_pgrids[image][index].end(),
                                            ppatch),
                                     m_pgrids[image][index].end());
    }
    
    for (int i = 0; i < static_cast<int>(ppatch->m_vimages.size()); ++i) {
        const int image = ppatch->m_vimages[i];
        const int& ix = ppatch->m_vgrids[i](0);
        const int& iy = ppatch->m_vgrids[i](1);
        const int index = iy * m_gwidths[image] + ix;
        m_vpgrids[image][index].erase(remove(m_vpgrids[image][index].begin(),
                                             m_vpgrids[image][index].end(),
                                             ppatch),
                                      m_vpgrids[image][index].end());
    }
}

int PatchManager::isVisible0(const Patch& patch, const int image, int& ix, int& iy, const float strict) {
    const Vector3f icoord = m_pmmvps.m_photoSet.project(image, patch.m_coord, m_pmmvps.m_level);
    ix = ((int)floorf(icoord(0) + 0.5f)) / m_pmmvps.m_csize;
    iy = ((int)floorf(icoord(1) + 0.5f)) / m_pmmvps.m_csize;
    
    return isVisible(patch, image, ix, iy, strict);
}

int PatchManager::isVisible(const Patch& patch, const int image, const int& ix, const int& iy, const float strict) {
    
    const int gheight = m_gheights[image];
    const int gwidth = m_gwidths[image];
    
    if (ix < 0 || gwidth <= ix ||
        iy < 0 || gheight <= iy) {
        return 0;
    }
    
    if (m_pmmvps.m_depth == 0) {
        return 1;
    }
    
    int isEmpty = 0;
    Ppatch dppatch = m_MAXDEPTH;
    const int index = iy * gwidth + ix;
    
    if (m_dpgrids[image][index] == m_MAXDEPTH) {
        isEmpty = 1;
    }
    else {
        dppatch = m_dpgrids[image][index];
    }
    
    if (isEmpty == 1) {
        return 1;
    }
    
    Vector4f ray = patch.m_coord - m_pmmvps.m_photoSet.m_photos[image].m_center;
    ray /= ray.norm();
    const float diff = ray.dot(patch.m_coord - dppatch->m_coord);
    const float factor = std::min(2.0, 2.0 + ray.dot(patch.m_normal));
    
    // ???reasons
    if (diff < m_pmmvps.m_optim.getUnit(image, patch.m_coord) * m_pmmvps.m_csize * strict * factor) {
        return 1;
    }
    else {
        return 0;
    }
}

void PatchManager::setScales(Patch& patch) const {
    const float unit = m_pmmvps.m_optim.getUnit(patch.m_images[0], patch.m_coord); // unit = 2 * p_z / (f_x + f_y) = p_x / x = p_y / y
    const float unit2 = 2.0f * unit; // ∆p_x or ∆p_y cooresponding to 2 pixels replacement. From my understanding, unit2 could be randomly initialized a value
    Vector4f ray = patch.m_coord - m_pmmvps.m_photoSet.m_photos[patch.m_images[0]].m_center;
    ray /= ray.norm();
    
    const int num = std::min(m_pmmvps.m_tau, static_cast<int>(patch.m_images.size()));
    
    // how many pixel difference per unit along the the light of sight
    // if the patch is moved along the line of sight by unit2, m_dscale records the pixel movement in all the visible views except for the reference
    for (int i = 1; i < num; ++i) {
        Vector3f diff = m_pmmvps.m_photoSet.project(patch.m_images[i], patch.m_coord, m_pmmvps.m_level) -
        m_pmmvps.m_photoSet.project(patch.m_images[i], patch.m_coord - unit2 * ray, m_pmmvps.m_level);
        patch.m_dscale += diff.norm();
    }
    
    // set m_dscale to the vertical distance where average pixel move is half pixel
    patch.m_dscale /= num - 1;
    patch.m_dscale = unit2 / patch.m_dscale;
    
    patch.m_ascale = atan(patch.m_dscale / (unit * m_pmmvps.m_wsize / 2.0f));
}

void PatchManager::computeNcc(Patch& patch) const {
    m_pmmvps.m_optim.computeWeights(patch);
    patch.m_ncc = 1.0f - m_pmmvps.m_optim.unrobustincc(m_pmmvps.m_optim.computeINCC(patch.m_coord, patch.m_normal, patch.m_images, 1));
}

void PatchManager::sortPatches(vector<Ppatch>& ppatches, const int ascend) const {
    int npatches = static_cast<int>(ppatches.size());
    if (npatches == 0) {
        return;
    }
    for (int n = 0; n < npatches; ++n) {
        if (ppatches[n]->m_ncc < 0.0f) {
            m_pmmvps.m_patchManager.computeNcc(*ppatches[n]);
        }
    }
    if (npatches == 1) {
        return;
    }
    for (int i = 0; i < npatches; ++i) {
        for (int j = i + 1; j < npatches; ++j) {
            if (ascend) {
                if (ppatches[i]->m_ncc > ppatches[j]->m_ncc) {
                    std::swap(ppatches[i], ppatches[j]);
                }
            }
            else {
                if (ppatches[i]->m_ncc < ppatches[j]->m_ncc) {
                    std::swap(ppatches[i], ppatches[j]);
                }
            }
        }
    }
}

void PatchManager::readPatches() {
//    for (int i = 0; i < m_pmmvps.m_nimages; ++i) {
//    const int image = m_pmmvps.m_images[i];
    char buffer[1024];
    sprintf(buffer, "%sply/%08d.patch", m_pmmvps.m_prefix.c_str(), 0);
    ifstream ifstr;
    ifstr.open(buffer);
    if (!ifstr.is_open()) {
        return;
    }
    
    string header;
    int pnum;
    ifstr >> header >> pnum;
//    cerr << image << " " << pnum << " patches" << endl;
    for (int p = 0; p < pnum; ++p) {
        Ppatch ppatch(new Patch());
        ifstr >> *ppatch;
        ppatch->m_fix = 0;
        ppatch->m_tmp = ppatch->score2(m_pmmvps.m_nccThreshold);
        ppatch->m_vimages.clear();
        image2index(*ppatch);
        if (ppatch->m_images.empty()) {
            return;
        }
        
        setGrids(*ppatch);
        addPatch(ppatch);
    }
//    }
    ifstr.close();
}

void PatchManager::readPatches(const int iter) {
	char buffer[1024];
	sprintf(buffer, "%sply/%08d.patch", m_pmmvps.m_prefix.c_str(), iter);
	ifstream ifstr;
	ifstr.open(buffer);
	if (!ifstr.is_open()) {
		return;
	}

	string header;
	int pnum;
	ifstr >> header >> pnum;
	//    cerr << image << " " << pnum << " patches" << endl;
	for (int p = 0; p < pnum; ++p) {
		Ppatch ppatch(new Patch());
		ifstr >> *ppatch;
		ppatch->m_fix = 0;
		ppatch->m_tmp = ppatch->score2(m_pmmvps.m_nccThreshold);
		ppatch->m_vimages.clear();
		image2index(*ppatch);
		if (ppatch->m_images.empty()) {
			return;
		}

		setGrids(*ppatch);
		//setVGrids(*ppatch);
		addPatch(ppatch);
	}
	ifstr.close();
}

void PatchManager::writePatches(const string prefix, bool bExportPLY, bool bExportPatch, bool bExportPSet) {
    collectPatches(1);
    
    if (bExportPLY)
    {
        char buffer[1024];
        sprintf(buffer, "%s.ply", prefix.c_str());
        writePly(m_ppatches, buffer);
    }
    
    if (bExportPatch)
    {
        char buffer[1024];
        sprintf(buffer, "%s.patch", prefix.c_str());
        ofstream ofstr;
        ofstr.open(buffer);
        ofstr << "PATCHES" << endl
			  << static_cast<int>(m_ppatches.size()) << endl;
        for (int p = 0; p < static_cast<int>(m_ppatches.size()); ++p) {
            Patch patch = *m_ppatches[p];
            index2image(patch);
            ofstr << patch << "\n";
        }
        ofstr.close();
    }
    
    if (bExportPSet)
    {
//        char buffer[1024];
//        sprintf(buffer, "%s.pset", prefix.c_str());
//        ofstream ofstr;
//        ofstr.open(buffer);
//        for (int p = 0; p < (int)m_ppatches.size(); ++p)
//            ofstr << m_ppatches[p]->m_coord[0] << ' '
//            << m_ppatches[p]->m_coord[1] << ' '
//            << m_ppatches[p]->m_coord[2] << ' '
//            << m_ppatches[p]->m_normal[0] << ' '
//            << m_ppatches[p]->m_normal[1] << ' '
//            << m_ppatches[p]->m_normal[2] << "\n";
//        ofstr.close();
    }
}

void PatchManager::writePly(const vector<Ppatch>& patches, const string filename) {
    ofstream ofstr;
    ofstr.open(filename.c_str());
    ofstr << "ply" << '\n'
    << "format ascii 1.0" << '\n'
    << "element vertex " << (int)patches.size() << '\n'
    << "property float x" << '\n'
    << "property float y" << '\n'
    << "property float z" << '\n'
    << "property float nx" << '\n'
    << "property float ny" << '\n'
    << "property float nz" << '\n'
    << "property uchar diffuse_red" << '\n'
    << "property uchar diffuse_green" << '\n'
    << "property uchar diffuse_blue" << '\n'
    << "end_header" << '\n';
    
    vector<Ppatch>::const_iterator bpatch = patches.begin();
    vector<Ppatch>::const_iterator bend = patches.end();
    
    while (bpatch != bend) {
        // Get color
        Vector3i color;
        
        const int mode = 0;
        // 0: color from images
        // 1: fix
        // 2: angle
        if (mode == 0) {
            int denom = 0;
            Vector3f colorf = Vector3f::Zero();;
            for (int i = 0; i < (int)(*bpatch)->m_images.size(); ++i) {
                const int image = (*bpatch)->m_images[i];
                if (m_pmmvps.m_nillums == 1) {
                    colorf += m_pmmvps.m_photoSet.getColor(image, (*bpatch)->m_coord, m_pmmvps.m_level);
                }
                else {
                    colorf += m_pmmvps.m_photoSet.getColor(image, (*bpatch)->m_coord, m_pmmvps.m_level, 0);
                }
                denom++;
            }
            colorf /= denom;
            color(0) = std::min(255,(int)floor(colorf(0) + 0.5f));
            color(1) = std::min(255,(int)floor(colorf(1) + 0.5f));
            color(2) = std::min(255,(int)floor(colorf(2) + 0.5f));
        }
        else if (mode == 1) {
//            if ((*bpatch)->m_tmp == 1.0f) {
//                color[0] = 255;
//                color[1] = 0;
//                color[2] = 0;
//            }
//            else {
//                color[0] = 255;
//                color[1] = 255;
//                color[2] = 255;
//            }
        }
        else if (mode == 2) {
//            float angle = 0.0f;
//            vector<int>::iterator bimage = (*bpatch)->m_images.begin();
//            vector<int>::iterator eimage = (*bpatch)->m_images.end();
//            
//            while (bimage != eimage) {
//                const int index = *bimage;
//                Vector4f ray = m_pmmvps.m_photoSet.m_photos[index].m_center - (*bpatch)->m_coord;
//                ray[3] = 0.0f;
//                unitize(ray);
//                
//                angle += acos(ray * (*bpatch)->m_normal);
//                ++bimage;
//            }
//            
//            angle = angle / (M_PI / 2.0f);
//            float r, g, b;
//            Image::gray2rgb(angle, r, g, b);
//            color[0] = (int)(r * 255.0f);
//            color[1] = (int)(g * 255.0f);
//            color[2] = (int)(b * 255.0f);
        }
        
        ofstr << (*bpatch)->m_coord(0) << ' '
        << (*bpatch)->m_coord(1) << ' '
        << (*bpatch)->m_coord(2) << ' '
        << (*bpatch)->m_normal(0) << ' '
        << (*bpatch)->m_normal(1) << ' '
        << (*bpatch)->m_normal(2) << ' '
        << color(0) << ' ' << color(1) << ' ' << color(2) << '\n';
        ++bpatch;
    }
    ofstr.close();
}

void PatchManager::writePly(const vector<Ppatch>& patches, const string filename, const vector<Vector3i>& colors) {
    ofstream ofstr;
    ofstr.open(filename.c_str());
    ofstr << "ply" << '\n'
    << "format ascii 1.0" << '\n'
    << "element vertex " << (int)patches.size() << '\n'
    << "property float x" << '\n'
    << "property float y" << '\n'
    << "property float z" << '\n'
    << "property float nx" << '\n'
    << "property float ny" << '\n'
    << "property float nz" << '\n'
    << "property uchar diffuse_red" << '\n'
    << "property uchar diffuse_green" << '\n'
    << "property uchar diffuse_blue" << '\n'
    << "end_header" << '\n';
    
    vector<Ppatch>::const_iterator bpatch = patches.begin();
    vector<Ppatch>::const_iterator bend = patches.end();
    vector<Vector3i>::const_iterator colorb = colors.begin();
    
    while (bpatch != bend) {
        ofstr << (*bpatch)->m_coord(0) << ' '
        << (*bpatch)->m_coord(1) << ' '
        << (*bpatch)->m_coord(2) << ' '
        << (*bpatch)->m_normal(0) << ' '
        << (*bpatch)->m_normal(1) << ' '
        << (*bpatch)->m_normal(2) << ' '
        << *colorb << '\n';
        ++bpatch;
        ++colorb;
    }
    ofstr.close();
}

// used in the test
void PatchManager::findNeighbors(const Patch& patch, vector<Ppatch>& neighbors, const float scale, const int margin, const int skipvis) {
    const float radius = 1.5 * margin * m_pmmvps.m_propagate.computeRadius(patch);
	//cerr << "radius: " << radius << endl;
    
    float unit = 0.0f;
    for (int i = 0; i < static_cast<int>(patch.m_images.size()); ++i) {
        unit += m_pmmvps.m_optim.getUnit(patch.m_images[i], patch.m_coord);
    }
    unit /= static_cast<int>(patch.m_images.size());
    unit *= m_pmmvps.m_csize;
    
    vector<int>::const_iterator bimage = patch.m_images.begin();
    vector<int>::const_iterator eimage = patch.m_images.end();
    vector<Vector2i>::const_iterator bgrid = patch.m_grids.begin();
    
    while (bimage != eimage) {
        const int image = *bimage;
        const int& ix = (*bgrid)(0);
        const int& iy = (*bgrid)(1);
        
        for (int i = -margin; i <= margin; ++i) {
            const int ytmp = iy + i;
            if (ytmp < 0 || m_pmmvps.m_patchManager.m_gheights[image] <= ytmp) {
                continue;
            }
            for (int j = -margin; j <= margin; ++j) {
                const int xtmp = ix + j;
                if (xtmp < 0 || m_pmmvps.m_patchManager.m_gwidths[image] <= xtmp) {
                    continue;
                }
                const int index = ytmp * m_pmmvps.m_patchManager.m_gwidths[image] + xtmp;
                vector<Ppatch>::const_iterator bpatch = m_pmmvps.m_patchManager.m_pgrids[image][index].begin();
                vector<Ppatch>::const_iterator epatch = m_pmmvps.m_patchManager.m_pgrids[image][index].end();
                
                while (bpatch != epatch) {
                    if (m_pmmvps.isNeighborRadius(patch, **bpatch, unit, m_pmmvps.m_neighborThreshold * scale, radius)) {
                        neighbors.push_back(*bpatch);
                    }
                    ++bpatch;
                }
                
                bpatch = m_pmmvps.m_patchManager.m_vpgrids[image][index].begin();
                epatch = m_pmmvps.m_patchManager.m_vpgrids[image][index].end();
                while (bpatch != epatch) {
                    if (m_pmmvps.isNeighborRadius(patch, **bpatch, unit, m_pmmvps.m_neighborThreshold * scale, radius)) {
                        neighbors.push_back(*bpatch);
                    }
                    ++bpatch;
                }
            }
        }
        ++bimage;
        ++bgrid;
    }
    
    sort(neighbors.begin(), neighbors.end());
    neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
}

// not used
float PatchManager::computeUnit(const Patch& patch) const {
    float unit = 0.0f;
    for (int i = 0; i < static_cast<int>(patch.m_images.size()); ++i) {
        unit += m_pmmvps.m_optim.getUnit(patch.m_images[i], patch.m_coord);
    }
    unit /= static_cast<int>(patch.m_images.size());
    unit *= m_pmmvps.m_csize;
    return unit;
}
