//
//  filter.cpp
//  PMMVPS
//
//  Created by KaiWu on Nov/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "filter.hpp"
#include "pmmvps.hpp"
#include <ctime>

using std::cerr;
using std::endl;
using std::flush;
using std::time_t;
using std::time;

Filter::Filter(PmMvps& pmmvps) : m_pmmvps(pmmvps) {
}

void Filter::init() {
}

void Filter::run() {
    setDepthMapsVGridsVPGridsAddPatchV(0);
	//string name = m_pmmvps.m_prefix + "ply/before_refine";
	//m_pmmvps.m_patchManager.writePatches(name, false, true, false);
    
    filterOutside();
    setDepthMapsVGridsVPGridsAddPatchV(1);
	//name = m_pmmvps.m_prefix + "ply/refined_patches_before_refine_f1";
	//m_pmmvps.m_patchManager.writePatches(name, true, false, false);
    
    filterExact();
    setDepthMapsVGridsVPGridsAddPatchV(1);
	//name = m_pmmvps.m_prefix + "ply/refined_patches_before_refine_f2";
	//m_pmmvps.m_patchManager.writePatches(name, true, false, false);

    filterNeighbor(1);
    setDepthMapsVGridsVPGridsAddPatchV(1);
	//name = m_pmmvps.m_prefix + "ply/refined_patches_before_refine_f3";
	//m_pmmvps.m_patchManager.writePatches(name, true, false, false);
    
	filterSmallGroups();
	setDepthMapsVGridsVPGridsAddPatchV(1);
	//name = m_pmmvps.m_prefix + "ply/refined_patches_before_refine_f4";
	//m_pmmvps.m_patchManager.writePatches(name, true, false, false);
}

void Filter::filterOutside() {
    time_t tv;
    time(&tv);
    time_t curtime = tv;
    cerr << "FilterOutside" << endl;
    // notice (1) to avoid removing patch of which m_fix = 1
    m_pmmvps.m_patchManager.collectPatches(1);
    
    const int psize = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
    m_gains.resize(psize);
    
    cerr << "main body" << flush;
    
	filterOutsideSub();
    
    // delete patches with positive m_gains
    int count = 0;
    float ave = 0.0f;
    float ave2 = 0.0f;
    int denom = 0;
    
    for (int p = 0; p < psize; ++p) {
        ave += m_gains[p];
        ave2 += m_gains[p] * m_gains[p];
        ++denom;
        
        if (m_gains[p] < 0.0f) {
            m_pmmvps.m_patchManager.removePatch(m_pmmvps.m_patchManager.m_ppatches[p]);
            ++count;
        }
    }
    
    if (denom == 0) {
        denom = 1;
    }
    
    ave /= denom;
    ave2 /= denom;
    ave2 = sqrtf(std::max(0.0f, ave2 - ave * ave));
    cerr << "Gain (ave/var): " << ave << " " << ave2 << endl;
    
    time(&tv);
    cerr << static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size()) << " -> "
        << static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size()) - count << " ("
        << 100 * static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size() - count) / static_cast<float>(m_pmmvps.m_patchManager.m_ppatches.size())
        << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << endl;
}

void Filter::filterOutsideSub() {
    const int sz = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
    
    for (int p = 0; p < sz; ++p) {
        Ppatch& ppatch = m_pmmvps.m_patchManager.m_ppatches[p];
        m_gains[p] = computeGain(*ppatch);
    }
}

float Filter::computeGain(const Patch& patch) {
    float gain = patch.score2(m_pmmvps.m_nccThreshold);
    
    const int sz = static_cast<int>(patch.m_images.size());
    for (int i = 0; i < sz; ++i) {
        const int& index = patch.m_images[i];
        const int& ix = patch.m_grids[i](0);
        const int& iy = patch.m_grids[i](1);
        const int index2 = iy * m_pmmvps.m_patchManager.m_gwidths[index] + ix;
        
        float maxpressure = 0.0f;
        for (int j = 0; j < static_cast<int>(m_pmmvps.m_patchManager.m_pgrids[index][index2].size()); ++j) {
            if (!m_pmmvps.isNeighbor(patch, *m_pmmvps.m_patchManager.m_pgrids[index][index2][j], m_pmmvps.m_neighborThreshold1)) {
                maxpressure = std::max(maxpressure, m_pmmvps.m_patchManager.m_pgrids[index][index2][j]->m_ncc - m_pmmvps.m_nccThreshold);
            }
        }
        gain -= maxpressure;
    }
    
	const int vsz = static_cast<int>(patch.m_vimages.size());
	for (int i = 0; i < vsz; ++i) {
		const int& index = patch.m_vimages[i];
		const float pdepth = m_pmmvps.m_photoSet.computeDepth(index, patch.m_coord);
		const int& ix = patch.m_vgrids[i](0);
		const int& iy = patch.m_vgrids[i](1);
		const int& index2 = iy * m_pmmvps.m_patchManager.m_gwidths[index] + ix;
		float maxpressure = 0.0f;

		for (int j = 0; j < static_cast<int>(m_pmmvps.m_patchManager.m_pgrids[index][index2].size()); ++j) {
			const float bdepth = m_pmmvps.m_photoSet.computeDepth(index, m_pmmvps.m_patchManager.m_pgrids[index][index2][j]->m_coord);
			if (pdepth < bdepth && !m_pmmvps.isNeighbor(patch, *m_pmmvps.m_patchManager.m_pgrids[index][index2][j], m_pmmvps.m_neighborThreshold1)) {
				maxpressure = std::max(maxpressure, m_pmmvps.m_patchManager.m_pgrids[index][index2][j]->m_ncc - m_pmmvps.m_nccThreshold);
			}
		}
		gain -= maxpressure;
	}
    
    return gain;
}

void Filter::filterExact() {
    time_t tv;
    time(&tv);
    time_t curtime = tv;
    cerr << "Filter Exact: " << flush;
    
    // cannot use (1) because we use patch.m_id to set newimages
    m_pmmvps.m_patchManager.collectPatches();
    const int psize = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
    
    m_newimages.clear();        m_newimages.resize(psize);
    m_removeimages.clear();     m_removeimages.resize(psize);
    m_newgrids.clear();         m_newgrids.resize(psize);
    m_removegrids.clear();      m_removegrids.resize(psize);
    
    filterExactSub();
    
    for (int p = 0; p < psize; ++p) {
        if (m_pmmvps.m_patchManager.m_ppatches[p]->m_fix) {
            continue;
        }
        
        for (int i = 0; i < static_cast<int>(m_removeimages[p].size()); ++i) {
            const int index = m_removeimages[p][i];
            const int ix = m_removegrids[p][i](0);
            const int iy = m_removegrids[p][i](1);
            const int index2 = iy * m_pmmvps.m_patchManager.m_gwidths[index] + ix;
            m_pmmvps.m_patchManager.m_pgrids[index][index2].
            erase(remove(m_pmmvps.m_patchManager.m_pgrids[index][index2].begin(),
                         m_pmmvps.m_patchManager.m_pgrids[index][index2].end(),
                         m_pmmvps.m_patchManager.m_ppatches[p]),
                  m_pmmvps.m_patchManager.m_pgrids[index][index2].end()); // we cannot use removePatch
        }
    }
    
    int count = 0;
    for (int p = 0; p < psize; ++p) {
        if (m_pmmvps.m_patchManager.m_ppatches[p]->m_fix) {
            continue;
        }
        
        Patch& patch = *m_pmmvps.m_patchManager.m_ppatches[p];
        patch.m_nimages = static_cast<int>(m_newimages[p].size());
        patch.m_images.swap(m_newimages[p]);
        patch.m_grids.swap(m_newgrids[p]);
        
        if (m_pmmvps.m_minImageNumThreshold <= static_cast<int>(patch.m_images.size())) {
            m_pmmvps.m_optim.setRefImage(patch);
            m_pmmvps.m_patchManager.setGrids(patch);
        }
        
		if (static_cast<int>(patch.m_images.size()) < m_pmmvps.m_minImageNumThreshold) {
            m_pmmvps.m_patchManager.removePatch(m_pmmvps.m_patchManager.m_ppatches[p]);
            ++count;
        }
    }
    time(&tv);
    cerr << (int)m_pmmvps.m_patchManager.m_ppatches.size() << " -> "
        << static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size()) - count << " ("
        << 100 * (static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size()) - count) / static_cast<float>(m_pmmvps.m_patchManager.m_ppatches.size())
        << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << endl;
}

void Filter::filterExactSub() {
//    const int psize = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
//    vector<vector<int> > newimages, removeimages;
//    vector<vector<Vector2i> > newgrids, removegrids;
//    newimages.resize(psize);    removeimages.resize(psize);
//    newgrids.resize(psize);     removegrids.resize(psize);
    
    for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
        cerr << "*" << flush;
        
        const int& w = m_pmmvps.m_patchManager.m_gwidths[image];
        const int& h = m_pmmvps.m_patchManager.m_gheights[image];
        int index = -1;
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                ++index;
                for (int i = 0; i < static_cast<int>(m_pmmvps.m_patchManager.m_pgrids[image][index].size()); ++i) {
                    const Patch& patch = *m_pmmvps.m_patchManager.m_pgrids[image][index][i];
                    if (patch.m_fix) {
                        continue;
                    }
                    int safe = 0;
                    
                    if (m_pmmvps.m_patchManager.isVisible(patch, image, x, y, m_pmmvps.m_neighborThreshold1)) {
                        safe = 1;
                    }
                    // use 4 neighbors
                    else if (0 < x && m_pmmvps.m_patchManager.isVisible(patch, image, x - 1, y, m_pmmvps.m_neighborThreshold1)) {
                        safe = 1;
                    }
                    else if (x < w - 1 && m_pmmvps.m_patchManager.isVisible(patch, image, x + 1, y, m_pmmvps.m_neighborThreshold1)) {
                        safe = 1;
                    }
                    else if (0 < y && m_pmmvps.m_patchManager.isVisible(patch, image, x, y - 1, m_pmmvps.m_neighborThreshold1)) {
                        safe = 1;
                    }
                    else if (y < h - 1 && m_pmmvps.m_patchManager.isVisible(patch, image, x, y + 1, m_pmmvps.m_neighborThreshold1)) {
                        safe = 1;
                    }
                    
                    if (safe) {
                        m_newimages[patch.m_id].push_back(image);
                        m_newgrids[patch.m_id].push_back(Vector2i(x, y));
                    }
                    else {
                        m_removeimages[patch.m_id].push_back(image);
                        m_removegrids[patch.m_id].push_back(Vector2i(x, y));
                    }
                }
            }
        }
    }
}

void Filter::filterNeighbor(const int times) {
    time_t tv;
    time(&tv);
    time_t curtime = tv;
    cerr << "FilterNeighbor:\t" << flush;
    
    // (1) to avoid removing patch of which m_fix = 1
    m_pmmvps.m_patchManager.collectPatches(1);
    if (m_pmmvps.m_patchManager.m_ppatches.empty()) {
        return;
    }
    
    const int psize = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
    m_rejects.resize(psize);
    for (int p = 0; p < psize; ++p) {
        m_rejects[p] = 0;
    }
    
    int count = 0;
    for (m_time = 0; m_time < times; ++m_time) {
        filterNeighborSub();

		vector<Ppatch>::iterator bpatch = m_pmmvps.m_patchManager.m_ppatches.begin();
		vector<Ppatch>::iterator epatch = m_pmmvps.m_patchManager.m_ppatches.end();
		vector<int>::iterator breject = m_rejects.begin();
    
		while (bpatch != epatch) {
			if (*breject == m_time + 1) {
				count++;
				m_pmmvps.m_patchManager.removePatch(*bpatch);
			}
			++bpatch;
			++breject;
		}
	}
    
    time(&tv);
    cerr << static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size()) << "-> "
        << static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size()) - count << " ("
        << 100 * (static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size()) - count) / static_cast<float>(m_pmmvps.m_patchManager.m_ppatches.size())
        << "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << endl;
}

void Filter::filterNeighborSub() {
    const int psize = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
    
    for (int p = 0; p < psize; ++p) {
        Ppatch& ppatch = m_pmmvps.m_patchManager.m_ppatches[p];
        if (m_rejects[p]) {
            continue;
        }
        
        vector<Ppatch> neighbors;
        m_pmmvps.m_patchManager.findNeighbors(*ppatch, neighbors, 4, 2, 1); // scale, margin, skipvis
        
        if (static_cast<int>(neighbors.size()) < 6) {
            m_rejects[p] = m_time + 1;
        }
        else if (filterQuad(*ppatch, neighbors)) {
            m_rejects[p] = m_time + 1;
        }
    }
}

int Filter::filterQuad(const Patch& patch, const vector<Ppatch>& neighbors) const {
    vector<vector<float> > A;
    vector<float> b, x;
    
    Vector4f xdir, ydir;
    ortho(patch.m_normal, xdir, ydir);
    
    const int nsize = static_cast<int>(neighbors.size());
    float h = 0.0f;
    for (int n = 0; n < nsize; ++n) {
        h += (neighbors[n]->m_coord - patch.m_coord).norm();
    }
    h /= nsize;
    
    A.resize(nsize);
    b.resize(nsize);
    
    vector<float> fxs, fys, fzs;
    fxs.resize(nsize);
    fys.resize(nsize);
    fzs.resize(nsize);
    for (int n = 0; n < nsize; ++n) {
        A[n].resize(5);
        Vector4f diff = neighbors[n]->m_coord - patch.m_coord;
        fxs[n] = diff.dot(xdir) / h;
        fys[n] = diff.dot(ydir) / h;
        fzs[n] = diff.dot(patch.m_normal);
        
        A[n][0] = fxs[n] * fxs[n];
        A[n][1] = fys[n] * fys[n];
        A[n][2] = fxs[n] * fys[n];
        A[n][3] = fxs[n];
        A[n][4] = fys[n];
        b[n] = fzs[n];
    }
    x.resize(5);
    lls(A, b, x);
    
    // compute residual divided by m_dscael;
    const int inum= std::min(m_pmmvps.m_tau, static_cast<int>(patch.m_images.size()));
    float unit = 0.0f;
    for (int i = 0; i < inum; ++i) {
        unit += m_pmmvps.m_optim.getUnit(patch.m_images[i], patch.m_coord);
    }
    unit /= inum;
    
    float residual = 0.0f;
    for (int n = 0; n < nsize; ++n) {
        const float res = x[0] * (fxs[n] * fxs[n]) +
                          x[1] * (fys[n] * fys[n]) +
                          x[2] * (fxs[n] * fys[n]) +
                          x[3] * fxs[n] +
                          x[4] * fys[n] - fzs[n];
        residual += fabsf(res) / unit;
    }
    residual /= (nsize - 5);
    
    if (residual < m_pmmvps.m_quadThreshold) {
        return 0;
    }
    else {
        return 1;
    }
}

void Filter::ortho(const Vector4f& z, Vector4f& x, Vector4f& y) const {
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

void Filter::lls(const vector<vector<float> >& A, const vector<float>& b, vector<float>& x) const {
    const int m = static_cast<int>(A.size());
    const int n = static_cast<int>(A[0].size());
    
    Eigen::MatrixXf matA(m, n);
    Eigen::VectorXf vecb(m);
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            matA(y, x) = A[y][x];
        }
    }
    for (int i = 0; i < m; ++i) {
        vecb(i) = b[i];
    }
    Eigen::VectorXf vecx = matA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vecb);
    
    for (int i = 0; i < n; ++i) {
        x[i] = vecx(i);
    }
}

void Filter::filterSmallGroups() {
	time_t tv;
	time(&tv);
	time_t curtime = tv;
	cerr << "FilterGroups:\t" << flush;
	m_pmmvps.m_patchManager.collectPatches();
	if (m_pmmvps.m_patchManager.m_ppatches.empty()) {
		return;
	}

	const int psize = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
	vector<int> label;
	label.resize(psize);
	fill(label.begin(), label.end(), -1);

	list<int> untouch;
	vector<Ppatch>::iterator bpatch = m_pmmvps.m_patchManager.m_ppatches.begin();
	for (int p = 0; p < psize; ++p) {
		untouch.push_back(p);
		(*bpatch)->m_flag = p;
		++bpatch;
	}

	int id = -1;
	while (!untouch.empty()) {
		const int pid = untouch.front();
		untouch.pop_front();

		if (label[pid] != -1) {
			continue;
		}

		label[pid] = ++id;
		list<int> ltmp;
		ltmp.push_back(pid);

		while (!ltmp.empty()) {
			const int ptmp = ltmp.front();
			ltmp.pop_front();
			filterSmallGroupsSub(ptmp, id, label, ltmp);
		}
	}
	++id;

	vector<int> size;
	size.resize(id);
	vector<int>::iterator biter = label.begin();
	vector<int>::iterator eiter = label.end();
	while (biter != eiter) {
		++size[*biter];
		++biter;
	}

	const int threshold = std::max(20, psize / 10000);
	cerr << threshold << endl;

	biter = size.begin();
	eiter = size.end();
	while (biter != eiter) {
		if (*biter < threshold) {
			*biter = 0;
		}
		else {
			*biter = 1;
		}
		++biter;
	}

	int count = 0;

	biter = label.begin();
	eiter = label.end();
	bpatch = m_pmmvps.m_patchManager.m_ppatches.begin();
	while (biter != eiter) {
		if ((*bpatch)->m_fix) {
			++biter;
			++bpatch;
			continue;
		}

		if (size[*biter] == 0) {
			m_pmmvps.m_patchManager.removePatch(*bpatch);
			++count;
		}
		++biter;
		++bpatch;
	}
	time(&tv);
	cerr << (int)m_pmmvps.m_patchManager.m_ppatches.size() << " -> "
		<< (int)m_pmmvps.m_patchManager.m_ppatches.size() - count << " ("
		<< 100 * ((int)m_pmmvps.m_patchManager.m_ppatches.size() - count) / (float)m_pmmvps.m_patchManager.m_ppatches.size()
		<< "%)\t" << (tv - curtime) / CLOCKS_PER_SEC << " secs" << endl;
}

void Filter::filterSmallGroupsSub(const int pid, const int id, vector<int>& label, list<int>& ltmp) const {
	// find neighbors of ptmp and set their ids
	const Patch& patch = *m_pmmvps.m_patchManager.m_ppatches[pid];
	const int index = patch.m_images[0];
	const int ix = patch.m_grids[0](0);
	const int iy = patch.m_grids[0](1);
	const int gwidth = m_pmmvps.m_patchManager.m_gwidths[index];
	const int gheight = m_pmmvps.m_patchManager.m_gheights[index];

	for (int y = -1; y <= 1; ++y) {
		const int iytmp = iy + y;
		if (iytmp < 0 || gheight <= iytmp) {
			continue;
		}
		for (int x = -1; x <= 1; ++x) {
			const int ixtmp = ix + x;
			if (ixtmp < 0 || gwidth <= ixtmp) {
				continue;
			}
			const int index2 = iytmp * gwidth + ixtmp;
			vector<Ppatch>::iterator bgrid = m_pmmvps.m_patchManager.m_pgrids[index][index2].begin();
			vector<Ppatch>::iterator egrid = m_pmmvps.m_patchManager.m_pgrids[index][index2].end();
			while (bgrid != egrid) {
				const int itmp = (*bgrid)->m_flag;
				if (label[itmp] != -1) {
					++bgrid;
					continue;
				}

				if (m_pmmvps.isNeighbor(patch, **bgrid, m_pmmvps.m_neighborThreshold2)) {
					label[itmp] = id;
					ltmp.push_back(itmp);
				}
				++bgrid;
			}
			bgrid = m_pmmvps.m_patchManager.m_vpgrids[index][index2].begin();
			egrid = m_pmmvps.m_patchManager.m_vpgrids[index][index2].end();
			while (bgrid != egrid) {
				const int itmp = (*bgrid)->m_flag;
				if (label[itmp] != -1) {
					++bgrid;
					continue;
				}

				if (m_pmmvps.isNeighbor(patch, **bgrid, m_pmmvps.m_neighborThreshold2)) {
					label[itmp] = id;
					ltmp.push_back(itmp);
				}
				++bgrid;
			}
		}
	}
}

void Filter::setDepthMaps() {
    for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
		fill(m_pmmvps.m_patchManager.m_dpgrids[image].begin(), m_pmmvps.m_patchManager.m_dpgrids[image].end(), m_pmmvps.m_patchManager.m_MAXDEPTH);
    }
    setDepthMapsSub();
}

void Filter::setDepthMapsSub() {
    for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
        const int gwidth = m_pmmvps.m_patchManager.m_gwidths[image];
        const int gheight = m_pmmvps.m_patchManager.m_gheights[image];
        
        vector<Ppatch>::iterator bpatch = m_pmmvps.m_patchManager.m_ppatches.begin();
        vector<Ppatch>::iterator epatch = m_pmmvps.m_patchManager.m_ppatches.end();
        
        while (bpatch != epatch) {
            Ppatch& ppatch = *bpatch;
            const Vector3f icoord = m_pmmvps.m_photoSet.project(image, ppatch->m_coord, m_pmmvps.m_level);
            
            const float fx = icoord(0) / m_pmmvps.m_csize;
            const int xs[2] = {(int)floorf(fx), (int)ceilf(fx)};
            const float fy = icoord(1) / m_pmmvps.m_csize;
            const int ys[2] = {(int)floorf(fy), (int)ceilf(fy)};
            const float depth = m_pmmvps.m_photoSet.m_photos[image].m_oaxis.dot(ppatch->m_coord);
            
            for (int j = 0; j < 2; ++j) {
                for (int i = 0; i < 2; ++i) {
                    if (xs[i] < 0 || gwidth <= xs[i] || ys[j] < 0 || gheight <= ys[j]) {
                        continue;
                    }
                    const int index = ys[j] * gwidth + xs[i];
                    
                    if (m_pmmvps.m_patchManager.m_dpgrids[image][index] == m_pmmvps.m_patchManager.m_MAXDEPTH) {
                        m_pmmvps.m_patchManager.m_dpgrids[image][index] = ppatch;
                    }
                    else {
                        const float dtmp = m_pmmvps.m_photoSet.m_photos[image].m_oaxis.dot(m_pmmvps.m_patchManager.m_dpgrids[image][index]->m_coord);
                        if (depth < dtmp) {
                            m_pmmvps.m_patchManager.m_dpgrids[image][index] = ppatch;
                        }
                    }
                }
            }
            ++bpatch;
        }
    }
}

void Filter::setDepthMapsVGridsVPGridsAddPatchV(const int additive) {
    m_pmmvps.m_patchManager.collectPatches();
    setDepthMaps();
    
    // clear m_vpgrids
    for (int image = 0; image < m_pmmvps.m_nimages; ++image) {
        vector<vector<Ppatch> >::iterator bvpatch = m_pmmvps.m_patchManager.m_vpgrids[image].begin();
        vector<vector<Ppatch> >::iterator evpatch = m_pmmvps.m_patchManager.m_vpgrids[image].end();
        while (bvpatch != evpatch) {
            (*bvpatch).clear();
            ++bvpatch;
        }
    }
    
    if (additive == 0) {
        vector<Ppatch>::iterator bpatch = m_pmmvps.m_patchManager.m_ppatches.begin();
        vector<Ppatch>::iterator epatch = m_pmmvps.m_patchManager.m_ppatches.end();
        while (bpatch != epatch) {
            (*bpatch)->m_vimages.clear();
            (*bpatch)->m_vgrids.clear();
            ++bpatch;
        }
    }
    
	setVGridsVPGrids();

	addPatchV();
}

void Filter::setVGridsVPGrids() {
    // add patches to m_vpgrids
    const int psize = static_cast<int>(m_pmmvps.m_patchManager.m_ppatches.size());
    for (int p = 0; p < psize; ++p) {
        Ppatch& ppatch = m_pmmvps.m_patchManager.m_ppatches[p];
        m_pmmvps.m_patchManager.setVImagesVGrids(ppatch);
    }
}

void Filter::addPatchV() {
    vector<Ppatch>::iterator bpatch = m_pmmvps.m_patchManager.m_ppatches.begin();
    vector<Ppatch>::iterator epatch = m_pmmvps.m_patchManager.m_ppatches.end();
    
    while (bpatch != epatch) {
        Ppatch& ppatch = *bpatch;
        vector<int>::iterator bimage = ppatch->m_vimages.begin();
        vector<int>::iterator eimage = ppatch->m_vimages.end();
        vector<Vector2i>::iterator bgrid = ppatch->m_vgrids.begin();
        
        while (bimage != eimage) {
            const int image = *bimage;
            const int& ix = (*bgrid)(0);
            const int& iy = (*bgrid)(1);
            const int index = iy * m_pmmvps.m_patchManager.m_gwidths[image] + ix;
            m_pmmvps.m_patchManager.m_vpgrids[image][index].push_back(ppatch);
            
            ++bimage;
            ++bgrid;
        }
        ++bpatch;
    }
}