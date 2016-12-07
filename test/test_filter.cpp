#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "CImg.h"
#include "pmmvps/pmmvps.hpp"

int main(int argc, char* argv[]) {
	Option option;
	option.init("C:/Users/Admin/Documents/Data/PMMVPS/", "option");

	PmMvps pmmvps;
	pmmvps.init(option);

	pmmvps.m_depth = 1;
	pmmvps.m_patchManager.readPatches(1);

	pmmvps.m_filter.run();
}