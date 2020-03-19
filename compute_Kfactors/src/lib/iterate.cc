#include "iterate.h"

  //**********************************************************************************************************************



 std::vector<Vector3d> iter::itermom(double max_mom)
  {

    std::vector<Vector3d> moms;
    Vector3d tmp(3,1);

    for(int i=-int(sqrt(max_mom)); i <= int(sqrt(max_mom)); i++){
      for(int j = -int(sqrt(max_mom)); j <= int(sqrt(max_mom)); j++){
        for(int k = -int(sqrt(max_mom)); k <= int(sqrt(max_mom)); k++){


          double mom_sq = (pow(i,2)+pow(j,2)+pow(k,2));
        
          // Cut-off on mom1
          if(max_mom >= mom_sq){
            tmp << -i,-j,-k;
            moms.push_back(tmp);
            
          }
        }
      }
    }

    return moms;
  }

  //**********************************************************************************************************************