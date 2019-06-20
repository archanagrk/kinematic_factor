
#include "phase.h"

double Round(double x) {return std::round(1000.0*x)/1000.0;}


  //**********************************************************************************************************************

    /* Composition of two Wigner-D matrices R'R, R'^-1R, R'R^-1 */

  //**********************************************************************************************************************

    std::complex <double> Ph::comp_Wigner_d(int twoJ, int twolam1, int twolam2,  double a1, double b1, double c1, double a2, double b2, double c2, int n){
        std::complex <double> D = 0;
            switch(n)
            {
                case 1: //R'^-1R
                    for(int i = -twoJ; i <= twoJ; i = i+2){
                        D += conj(Hadron::Wigner_D(twoJ, i, twolam1, a1, b1, c1))*
                            Hadron::Wigner_D(twoJ, i, twolam2, a2, b2, c2);
                    }
                    break;
                    
                case 2: //R'R^-1
                    for(int i = -twoJ; i <= twoJ; i = i+2){
                        D += Hadron::Wigner_D(twoJ, twolam1, i, a1, b1, c1)*
                        conj(Hadron::Wigner_D(twoJ, twolam2, i, a2, b2, c2));

                    }
                    break;
                    
                default: //R'R
                    for(int i = -twoJ; i <= twoJ; i = i+2){
                        D += Hadron::Wigner_D(twoJ, twolam1, i, a1, b1, c1)*
                            Hadron::Wigner_D(twoJ, i, twolam2, a2, b2, c2);
                    }
                    break;
                    
                    
            }
            
            return D;
    }

  //**********************************************************************************************************************

    /* Compute the phase factor after rotation of a helicity operator or state from Appendix A of Shulz paper */

  //**********************************************************************************************************************


Ph::phChars Ph::phaseFactor(int twoJ1, int twoJ2, int twoJCurr, Eigen::Vector3d mom1, Eigen::Vector3d mom2, bool compute){

    Ph::phChars out;
    map < Ph::tripKey , complex<double>> phase_hel;

    //if you want to do the phase correction
    //**********************************************************************************************************************
    if(compute) { 

        XMLArray::Array<int> mom2_r(3); mom2_r[0] = mom2(0); mom2_r[1] = mom2(1); mom2_r[2] = mom2(2);
        mom2_r = Hadron::canonicalOrder(mom2_r);   //ref mom is wrt mom2 - the vector particle - change it to canon and find R
        Eigen::Vector3d mom2_ref;
        mom2_ref  << mom2_r[0],mom2_r[1],mom2_r[2];
        
        
        std::vector<double> r2 = LittleGrp::refAngles(mom2_ref);
        std::vector<double> r_mom2 = LittleGrp::refAngles(mom2);
        std::vector<double> r_mom1 = LittleGrp::refAngles(mom1);


        Eigen::MatrixXd R_mom2 = Rot::eulerRotMat(r_mom2[0],r_mom2[1],r_mom2[2]);
        Eigen::MatrixXd R_ref_mom2 = Rot::eulerRotMat(r2[0],r2[1],r2[2]);
        Eigen::MatrixXd R_without_rounding = R_ref_mom2*R_mom2.inverse();

        Eigen::MatrixXd R = R_without_rounding.unaryExpr(&Round);


        Eigen::VectorXd  mom1_f(4); mom1_f << 0.0 , mom1(0), mom1(1), mom1(2);
        Eigen::VectorXd n_mom1_f = R*mom1_f;


        Eigen::Vector3d n_mom1;
        n_mom1 << round(n_mom1_f(1,0)), round(n_mom1_f(2,0)), round(n_mom1_f(3,0));


        Eigen::Vector3d n_mom_curr; n_mom_curr =  mom2_ref - n_mom1;
        Eigen::Vector3d mom_curr; mom_curr = mom2 - mom1;


        Eigen::VectorXd  mom_curr_f(4); mom_curr_f << 0.0 , mom_curr(0), mom_curr(1), mom_curr(2);
        Eigen::VectorXd n_mom_curr_f = R*mom_curr_f;


        Eigen::Vector3d n_mom_curr_check;
        n_mom_curr_check << round(n_mom_curr_f(1,0)), round(n_mom_curr_f(2,0)), round(n_mom_curr_f(3,0));




        std::vector<double> r_n_mom1 = LittleGrp::refAngles(n_mom1);


        /* Specialization for canonical mom */

    
        double mom1_sq = mom1.squaredNorm();
        double mom2_sq = mom2.squaredNorm();
        double mom_curr_sq = mom_curr.squaredNorm();


        if((n_mom1 == mom1 && mom2_ref == mom2)) {

            phase_hel = Ph::cnst_phase(twoJ1, twoJ2, twoJCurr); 

            out.mom1 = mom1; //mom1
            out.mom2 = mom2; // mom2
            out.lam_phase = phase_hel; //phase
            out.r = MatrixXd::Identity(4, 4); // the rotation is identity

            return out;

        }

        /* specialization when photon is at rest */



        std::vector<double> r_mom_curr = LittleGrp::refAngles(mom_curr);
        std::vector<double> r_n_mom_curr = LittleGrp::refAngles(n_mom_curr);

    

        Eigen::MatrixXd R_mom_curr = Rot::eulerRotMat(r_mom_curr[0],r_mom_curr[1],r_mom_curr[2]);
        Eigen::MatrixXd R_n_mom_curr = Rot::eulerRotMat(r_n_mom_curr[0],r_n_mom_curr[1],r_n_mom_curr[2]);

        Eigen::MatrixXd R_mom1 = Rot::eulerRotMat(r_mom1[0],r_mom1[1],r_mom1[2]);
        Eigen::MatrixXd R_n_mom1 = Rot::eulerRotMat(r_n_mom1[0],r_n_mom1[1],r_n_mom1[2]);

        //cout << "R" << endl << R << endl;
        
        phase_hel = Ph::calc_phase(twoJ1, twoJ2, twoJCurr, mom1_sq, mom2_sq, mom_curr_sq, r_mom1,r_n_mom1, r_mom2,r2, r_mom_curr, r_n_mom_curr);     
                                        
        
        out.mom1 = n_mom1; //new mom1
        out.mom2 = mom2_ref; // new mom2 - the canonical
        out.lam_phase = phase_hel; //the phase
        out.r = R; // the rotation performed


    }
    //dont factor in the phases
    //**********************************************************************************************************************
    else{

        phase_hel = Ph::cnst_phase(twoJ1, twoJ2, twoJCurr); 

        out.mom1 = mom1; //mom1
        out.mom2 = mom2; // mom2
        out.lam_phase = phase_hel; //phase
        out.r = MatrixXd::Identity(4, 4);; // the rotation is identity

    }

    return out;


}


//**********************************************************************************************************************

// Calculate - D1,D2 and D_curr - the phases of the source, sink and current

//**********************************************************************************************************************


map < Ph::tripKey , complex<double>> Ph::calc_phase(int twoJ1, int twoJ2, int twoJCurr, double mom1_sq, 
            double mom2_sq, double mom_curr_sq, vector<double> r_mom1, vector<double> r_n_mom1, vector<double> r_mom2, vector<double> r2, vector<double> r_mom_curr, vector<double> r_n_mom_curr)
{

    map < Ph::tripKey , complex<double>> phase_hel; 
    Ph::tripKey lambd; //helicities for the 3 particles
    complex<double> ph; //phase

    for(int two_abs_lam1 = -twoJ1; two_abs_lam1 <= twoJ1; two_abs_lam1 = two_abs_lam1 +2 ){
        for(int two_abs_lam2 = -twoJ2; two_abs_lam2 <= twoJ2; two_abs_lam2 = two_abs_lam2 +2 ){
            for(int two_abs_lamCurr = -twoJCurr; two_abs_lamCurr <= twoJCurr; two_abs_lamCurr = two_abs_lamCurr +2 ){
        
            complex<double> D2 = 0;
            complex<double> D1 = 0;
            complex<double> DCurr = 0;
        
        
    
            if(mom2_sq){
                for(int i = -twoJ2; i <= twoJ2; i = i+2){
                    for(int j = -twoJ2; j <= twoJ2; j = j+2){
                
                        D2 += conj(Hadron::Wigner_D(twoJ2, i, two_abs_lam2, r2[0], r2[1], r2[2])) * Ph::comp_Wigner_d(twoJ2, i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJ2,  j, two_abs_lam2, r_mom2[0], r_mom2[1], r_mom2[2]);
                    }
                }
            }
            else{D2 = 1;}
            
            if(mom1_sq){
                for(int i = -twoJ1; i <= twoJ1; i = i+2){
                    for(int j = -twoJ1; j <= twoJ1; j = j+2){
                    
                        D1 += conj(Hadron::Wigner_D(twoJ1, i, two_abs_lam1, r_n_mom1[0], r_n_mom1[1], r_n_mom1[2])) *  Ph::comp_Wigner_d(twoJ1, i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJ1,  j, two_abs_lam1, r_mom1[0], r_mom1[1], r_mom1[2]);
                    }
                }
            }
                
            else{D1 = 1;}
                
                
            if(mom_curr_sq){
                for(int i = -twoJCurr; i <= twoJCurr; i = i+2){
                    for(int j = -twoJCurr; j <= twoJCurr; j = j+2){

                        DCurr += conj(Hadron::Wigner_D(twoJCurr, i, two_abs_lamCurr, r_n_mom_curr[0], r_n_mom_curr[1], r_n_mom_curr[2])) *  Ph::comp_Wigner_d(twoJCurr, i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJCurr,  j, two_abs_lamCurr, r_mom_curr[0], r_mom_curr[1], r_mom_curr[2]);
                    }
                    
                }
            }
            else{DCurr = 1;}
            

            lambd = std::make_tuple(two_abs_lam1, two_abs_lam2 , two_abs_lamCurr);
            //ph =  D1* conj(D2)*DCurr;
            
            ph.real(KfUt::truncate(std::real(D1* conj(D2)), 6)); // the current is a four vector dont add the phase but use the rotation matrix
            ph.imag(KfUt::truncate(std::imag(D1* conj(D2)), 6));
            //if(!twoJ1){cout << "lam" << two_abs_lamCurr << "D" << DCurr << "D" << D1 << "D" << D2 << "ph" << ph << endl;}
            //ph = 1;
            phase_hel.insert( std::make_pair( lambd, ph  ) );

            }
        }
    }
    return phase_hel; 
}; 


map < Ph::tripKey , complex<double>> Ph::cnst_phase(int twoJ1, int twoJ2, int twoJCurr)
{

    map < Ph::tripKey , complex<double>> phase_hel; 
    Ph::tripKey lambd; //helicities for the 3 particles
    complex<double> ph; //phase

    for(int two_abs_lam1 = -twoJ1; two_abs_lam1 <= twoJ1; two_abs_lam1 = two_abs_lam1 +2 ){
        for(int two_abs_lam2 = -twoJ2; two_abs_lam2 <= twoJ2; two_abs_lam2 = two_abs_lam2 +2 ){
            for(int two_abs_lamCurr = -twoJCurr; two_abs_lamCurr <= twoJCurr; two_abs_lamCurr = two_abs_lamCurr +2 ){

                lambd = std::make_tuple(two_abs_lam1, two_abs_lam2 , two_abs_lamCurr);
                ph = 1;
                phase_hel.insert( std::make_pair( lambd, ph  ) );

            }
        }
    }

    return phase_hel;

};



