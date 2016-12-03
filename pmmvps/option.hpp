//
//  option.hpp
//  PMMVPS
//
//  Created by KaiWu on Oct/24/16.
//  Copyright © 2016 KaiWu. All rights reserved.
//

#ifndef option_hpp
#define option_hpp

#include <string>
#include <vector>
#include <map>

using std::string;
using std::vector;
using std::map;

struct Option {
public:
    Option();
    void init(const string prefix, const string option);
    
    // number of views
    int m_nimages;
    // number of illuminations
    int m_nillums;
    // number of levels
    int m_level;
    // cell size
    int m_csize;
    // ncc threshold
    float m_nccThreshold;
    // window size
    int m_wsize;
    // minimum number of images that a patch should be visible
    int m_minImageNum;
    // number of CPU
    int m_cpu;
    //
    int m_setEdge;
    //
    int m_useBound;
    // use data of visible views
    int m_useVisData;
    // use a sequence of images
    int m_sequence;
    // minimum baseline angle
    float m_maxAngleThreshold;
    // quad threshold
    float m_quadThreshold;
    
    // root directory
    string m_prefix;
    // option file
    string m_option;
    // type of image indexes
    // -1: firstimage lastimage
    // nimages: a list of image indexes
    int m_flag;
    // images
    vector<int> m_images;
    // image-index pair
    map<int, int> m_dict;
    // visible views, from SfM
    vector<vector<int> > m_visdata;
    // visible views. index, not image
    vector<vector<int> > m_visdata2;
    
protected:
    void initVisdata();
};

#endif /* option_hpp */
