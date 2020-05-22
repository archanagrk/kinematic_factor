#include "iterate.h"

  //**********************************************************************************************************************

    /* function to iterate over each index the three-mometna  */

  //**********************************************************************************************************************

 std::vector<Vector3d> iter::itermom(double max_mom, double min_mom, ADAT::Array1dO<ADAT::Array1dO<int>> omit_mom, bool canonical)
  {

    std::vector<Vector3d> moms;
    Vector3d mom(3,1);
    Vector3d remove(3,1);

    for(int i=-int(sqrt(max_mom)); i <= int(sqrt(max_mom)); i++){
      for(int j = -int(sqrt(max_mom)); j <= int(sqrt(max_mom)); j++){
        for(int k = -int(sqrt(max_mom)); k <= int(sqrt(max_mom)); k++){


          double mom_sq = (pow(i,2)+pow(j,2)+pow(k,2));
        
          // Cut-off on mom1
          if(max_mom >= mom_sq && mom_sq >= min_mom){
            mom << -i,-j,-k;

            if(canonical){
              XMLArray::Array<int> mom_tmp(3);
              mom_tmp[0] = mom(0); mom_tmp[1] = mom(1); mom_tmp[2] = mom(2);
              XMLArray::Array<int> mom_can = Hadron::canonicalOrder(mom_tmp);

              if(mom_tmp == mom_can){
                if(omit_mom.size() > 0){
                  bool ispresent = false;

                  for(int omit = 1; omit <= omit_mom.size(); omit++){
                    remove << omit_mom[omit][1], omit_mom[omit][2], omit_mom[omit][3];
                    
                    if(mom == remove){ispresent = true;}
                  }

                  if(!ispresent){moms.push_back(mom);}
                  
                }

                else{moms.push_back(mom);}                
              }
            }

            else{
              if(omit_mom.size() > 0){
                bool ispresent = false;

                for(int omit = 1; omit <= omit_mom.size(); omit++){
                  remove << omit_mom[omit][1], omit_mom[omit][2], omit_mom[omit][3];

                  if(mom == remove){ispresent = true;}
                }

                if(!ispresent) moms.push_back(mom);
                
              }

              else{moms.push_back(mom);}
              //std::cout <<  mom.transpose() << std::endl;
            }
          }
        }
      }
    }

    if(moms.size() == 0){std::cout << "No moms to iterate over" << std::endl; exit(1);}

    return moms;
  }

  //**********************************************************************************************************************