#include "naming.h"

  //**********************************************************************************************************************



string naming::name(int npt, Ph::tripKey two_abs_lam, Vector3d mom1 , Vector3d mom_curr, Vector3d mom3, irrep_label rep1, irrep_label rep_curr, irrep_label rep3, string LG1, string LG_curr, string LG3, string lev1, string lev3)
    {
        string name, tmp_nm; 

        for(int pts = 0; pts < npt; pts++){
            if(pts == 2){
            string mst = "_p" + std::to_string(int(mom1(0))) + std::to_string(int(mom1(1)))+ std::to_string(int(mom1(2)));

            int mom1_sq = mom1.squaredNorm();

            if(mom1_sq != 0){
                tmp_nm = "xH"+ std::to_string(get<0>(two_abs_lam)/2) + LG1 + rep1.irrep + "r" + std::to_string(rep1.row) + mst;
            }
            else{
                tmp_nm =  "x" + rep1.irrep + "r" + std::to_string(rep1.row) + mst;
            }

            tmp_nm += "_" + lev1;
            
            }

            if(pts == 1){
            string mst = "_p" + std::to_string(int(mom_curr(0))) + std::to_string(int(mom_curr(1))) + std::to_string(int(mom_curr(2)));

            int mom_curr_sq = mom_curr.squaredNorm();

            if(mom_curr_sq != 0){
                tmp_nm = "xH"+ std::to_string(get<2>(two_abs_lam)/2) + LG_curr + rep_curr.irrep + "r" + std::to_string(rep_curr.row) + mst;
            }
            else{
                tmp_nm = "x" + rep_curr.irrep + "r" + std::to_string(rep_curr.row) + mst;
            }

            }
            if(pts == 0){
            string mst = "_p" + std::to_string(int(mom3(0))) + std::to_string(int(mom3(1))) + std::to_string(int(mom3(2)));

            int mom2_sq = mom3.squaredNorm();

            if(mom2_sq != 0){
                tmp_nm = "H"+ std::to_string(get<1>(two_abs_lam)/2) + LG3 + rep3.irrep + "r" + std::to_string(rep3.row) + mst;
            }
            else{
                tmp_nm = rep3.irrep + "r" + std::to_string(rep3.row) + mst;
            }

            tmp_nm += "_" + lev3;

            }

            name += tmp_nm;
        }

        return name;
    }

  //**********************************************************************************************************************