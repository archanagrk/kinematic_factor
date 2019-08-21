
#include "subd_pol_vec.h"

//**********************************************************************************************************************

//**********************************************************************************************************************

  /* Perform the subductions by adding over the J_z and abs_lams  */

//**********************************************************************************************************************


//**********************************************************************************************************************

/* The product of the subduction and polarization vectors and finally a sum over +-lambda */


//In the case of scalar there is no polarization tensor so multiply by identity
//In the case of vector multiple by a polarization 4-vector
//mo support of tensors yet

//**********************************************************************************************************************

map< int, Eigen::MatrixXcd > SubdPol::Subduce_with_pol(double& mom_sq, double& mass_sq,  int& twoJ, const irrep_label& irrep,
                                   const string& little_group,double R1_phi, double R1_theta, double R1_psi, bool curr){
    
    /* Variables */
    map< int, complex<double> > Sub;
    map< int, Eigen::MatrixXcd > Sub_with_pol;
    typedef std::complex<double> cd;
    complex<double> z_i(0.,1.);
    double zero = 0.0;
    


    //At rest subductions    
    if(mom_sq == 0.0){

            Sub = Subd::subduce_oct(irrep);

            if(twoJ == 2)
            {
                for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++)
                {   //cout << "subinside: " << (it->second) << " " << (it->first) << endl;
                    //cout << " " << PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi, curr) << endl;
                    Sub_with_pol.insert(make_pair( (it->first) , (it->second) * KfUt::Gmunu() * PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi, curr)));
                }
                
            
            }

            else if(twoJ == 0){
                Eigen::MatrixXcd unit_vec(4,1);
                unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
                for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    
                    Sub_with_pol.insert(make_pair( (it->first) , (it->second) * unit_vec ));
                    
                }
            
            }

            else{
                    cerr << "Tensors and fermions not coded" << endl; exit(1);
                
            }

        }
        
    //In flight    
    else{

        if(twoJ%2){cerr << "Fermions not coded" << endl; exit(1);} //fermions
        
        else{
            Sub = Subd::subduce_lg_boson(irrep, little_group);

            if(twoJ == 2){
                for(map< int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    //cout << "subinside: " << (it->second) << " " << (it->first) << endl;
                    //cout << "pols " << PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi, curr) << endl;
                    Sub_with_pol.insert(make_pair( (it->first) , (it->second) * KfUt::Gmunu() * PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi, curr) ));}}

            else if(twoJ == 0){
                Eigen::MatrixXcd unit_vec(4,1);
                unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
                for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    Sub_with_pol.insert(make_pair( (it->first) ,(it->second) * unit_vec )) ;}}


            else{cerr << "Tensors not coded" << endl; exit(1);}     

        }  //bosons
        
        }
            
            

    return Sub_with_pol;
    
};


