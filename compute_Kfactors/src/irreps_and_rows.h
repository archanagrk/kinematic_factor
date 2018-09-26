#ifndef __IRREPANDROWS__
#define __IRREPANDROWS__

//std lib
#include <cmath>
#include "math.h"
#include <string>
#include <vector>
#include <iostream>



using namespace std;


namespace IrrepName {
  std::vector<std::string>  getIrrep(int& twoJ, int& P, string& lg);
  int irrepRows(string& irrep);
}


#endif
