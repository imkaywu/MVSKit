//
//  option.cpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#include "option.hpp"
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

using std::cerr;
using std::endl;
using std::ifstream;

Option::Option() {
    m_level = 1;
    m_csize = 2;
    m_wsize = 7;
    m_nccThreshold = 0.7f;
    m_minImageNum = 3;
    m_cpu = 4;
    m_setEdge = 0;
    m_useBound = 0;
    m_useVisData = 0;
    m_sequence = -1;
    m_flag = -10;
    m_maxAngleThreshold = 10.0f * M_PI / 180.0f;
    m_quadThreshold = 2.5f;
}

void Option::init(const string prefix, const string option) {
    m_prefix = prefix;
    m_option = option;
    ifstream ifstr;
    string optionFile = prefix + option;
    ifstr.open(optionFile.c_str());
    while (1) {
        string name;
        ifstr >> name;
        if(ifstr.eof()) {
            break;
        }
        if (name[0] == '#') {
            char buffer[1024];
            ifstr.putback('#');
            ifstr.getline(buffer, 1024);
            continue;
        }
        if (name == "image") {
            ifstr >> m_nimages;
        }
        else if (name == "illum") {
            ifstr >> m_nillums;
        }
        else if (name == "level") {
            ifstr >> m_level;
        }
        else if (name == "csize") {
            ifstr >> m_csize;
        }
        else if (name == "threshold") {
            ifstr >> m_nccThreshold;
        }
        else if (name == "wsize") {
            ifstr >> m_wsize;
        }
        else if (name == "minImageNum") {
            ifstr >> m_minImageNum;
        }
        else if (name == "CPU") {
            ifstr >> m_cpu;
        }
        else if (name == "setEdge") {
            ifstr >> m_setEdge;
        }
        else if (name == "useBound") {
            ifstr >> m_useBound;
        }
        else if (name == "useVisData") {
            ifstr >> m_useVisData;
        }
        else if (name == "sequence") {
            ifstr >> m_sequence;
        }
        else if (name == "maxAngle") {
            ifstr >> m_maxAngleThreshold;
            m_maxAngleThreshold *= M_PI / 180.0f;
        }
        else if (name == "quad") {
            ifstr >> m_quadThreshold;
        }
        else if (name == "images") {
            ifstr >> m_flag;
            if (m_flag == -1) {
                int firstImage, lastImage;
                ifstr >> firstImage >> lastImage;
                for (int i = firstImage; i < lastImage; ++i) {
                    m_images.push_back(i);
                }
            }
            else if(0 < m_flag) {
                for (int i = 0; i < m_flag; ++i) {
                    int index;
                    ifstr >> index;
                    m_images.push_back(index);
                }
            }
            else {
                cerr << "flag is not valid: " << m_flag << endl;
                exit(1);
            }
        }
        else {
            cerr << "Unrecognizable option: " << name << endl;
            exit(1);
        }
    }
    ifstr.close();
    
    // m_flag is initialized to -10
    if (m_flag == -10) {
        cerr << "m_flag not specified: " << m_flag << endl;
        exit(1);
    }
    
    for (int i = 0; i < static_cast<int>(m_images.size()); ++i) {
        m_dict[m_images[i]] = i;
    }
    
    initVisdata();
    
    cerr << "--------------------------------------------------" << endl
    << "--- Summary of specified options ---" << endl;
    cerr << "# of images: " << (int)m_images.size();
    if (m_flag == -1)
        cerr << " (range specification)" << endl;
    else
        cerr << " (enumeration)" << endl;
    
    cerr << "level: " << m_level << "  csize: " << m_csize << endl
        << "nccThreshold: " << m_nccThreshold << "  wsize: " << m_wsize << endl
        << "minImageNum: " << m_minImageNum << "  CPU: " << m_cpu << endl
        << "useVisData: " << m_useVisData << "  sequence: " << m_sequence << endl;
    cerr << "--------------------------------------------------" << endl;
}

void Option::initVisdata() {
    if (m_useVisData == 0) {
        const int nimages = static_cast<int>(m_images.size());
//        m_visdata.resize(nimages);
        m_visdata2.resize(nimages);
        for (int y = 0; y < nimages; ++y) {
            for (int x = 0; x < nimages; ++x) {
                if (x != y) {
                    m_visdata2[y].push_back(x);
                }
                else {
                    // do something
                }
            }
        }
    }
    else {
        // use vis.dat
    }
}