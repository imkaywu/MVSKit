//
//  patch.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/19/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "patch.hpp"

Patch::Patch () {
    m_ncc = -1.0f;
    m_nimages = 0;
    m_iter = 0;
    m_collected = 0;
    m_flag = 0;
    m_dflag = 0;
    m_fix = 0;
    m_dscale = 0.0f; // doesn't have to, initialized implicitly to 0.0f
    m_ascale = 0.0f;
}

//float Patch::score (const float threshold) const {
//    return std::max(0.0f, m_ncc - threshold) * static_cast<int>(m_images.size());
//}

float Patch::score2(const float threshold) const {
    return std::max(0.0f, m_ncc - threshold) * static_cast<int>(m_images.size());
}

std::istream& operator >>(std::istream& istr, Patch& rhs) {
    std::string header;
    int itmp;
    istr >> header >> rhs.m_coord >> rhs.m_normal >> rhs.m_ncc
        >> rhs.m_dscale >> rhs.m_ascale;
    
    if (header == "PATCHA") {
        int type;
        Vector4f dir;
        istr >> type >> dir;
    }
    
    istr >> itmp;
    rhs.m_images.resize(itmp);
    for (int i = 0; i < itmp; ++i) {
        istr >> rhs.m_images[i];
    }
    
    istr >> itmp;
    rhs.m_vimages.resize(itmp);
    for (int i = 0; i < itmp; ++i) {
        istr >> rhs.m_vimages[i];
    }
    
    return istr;
}

std::ostream& operator <<(std::ostream& ostr, const Patch& rhs) {
    ostr << "PATCHES" << endl
        << rhs.m_coord << endl
        << rhs.m_normal << endl
        << rhs.m_ncc << ' '
        << rhs.m_dscale << ' '
        << rhs.m_ascale << endl
        << static_cast<int>(rhs.m_images.size()) << endl;
    
    for (int i = 0; i < static_cast<int>(rhs.m_images.size()); ++i) {
        ostr << rhs.m_images[i] << ' ';
    }
    ostr << endl;
    
    for (int i = 0; i < static_cast<int>(rhs.m_vimages.size()); ++i) {
        ostr << rhs.m_vimages[i] << ' ';
    }
    ostr << endl;
    
    return ostr;
}

std::istream& operator >>(std::istream& istr, Vector4f& v) {
    return istr >> v(0) >> v(1) >> v(2) >> v(3);
}

std::ostream& operator <<(std::ostream& ostr, const Vector4f& v) {
    return ostr << v(0) << " " << v(1) << " " << v(2) << " " << v(3);
}

