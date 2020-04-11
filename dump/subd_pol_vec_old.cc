
#include "subd_pol_vec.h"

//**********************************************************************************************************************

//**********************************************************************************************************************

  /* Perform the subductions by adding over the J_z and abs_lams  */

//**********************************************************************************************************************


sub_hel SubdPol::Subduce_all(double& mom_sq, double& mass_sq,  int& twoJ, const irrep_label& irrep,
                                                   const string& little_group,double R1_phi, double R1_theta, double R1_psi){


  map< int, complex<double> > Sub;
  Eigen::MatrixXcd sum = Eigen::MatrixXcd::Zero(4,1);
  typedef std::complex<double> cd;
  sub_hel out;
    int two_hel;

  switch(int(mom_sq)){

  case 0:{
        Sub = Subd::subduce_oct(irrep);
        if(twoJ == 2){
          for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
            //cout << "lam" <<  (it->first) << "lam" << "\n";
            //cout << irrep.irrep << "\n";
            //cout << PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi) << "\n";
            sum +=  (it->second) * PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi); two_hel = it->first;}}
        else{
        Eigen::MatrixXcd unit_vec(4,1);
        unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
        for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
          sum +=  (it->second) * unit_vec; two_hel = it->first;}}
        break;}


  default:{
        if(twoJ%2){
          Eigen::MatrixXcd unit_vec(4,1);
          unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
          Sub  = Subd::subduce_lg_fermion(irrep,little_group);
          for(map< int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
            sum +=  (it->second) * unit_vec; two_hel = it->first;}
   } //fermions

        else{
          Sub = Subd::subduce_lg_boson(irrep, little_group);
          if(twoJ == 2){
            for(map< int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                //cout << irrep.irrep << "\n";
                //cout << PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi) << "\n";
              //cout << "lam" <<  (it->first) << "lam" << "\n";
              //cout << "Polstart" <<getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi)<< "Polend";
              sum +=  (it->second) * PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi);  two_hel = it->first;}}
          else{
            Eigen::MatrixXcd unit_vec(4,1);
            unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
            for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
              sum +=  (it->second) * unit_vec;  two_hel = it->first;}}
  }  //bosons
        break;}


}
  out.sum = sum;
  out.two_hel = two_hel;
  return out;

};

//**********************************************************************************************************************

/* The product of the subduction and polarization vectors and finally a sum over +-lambda */


//In the case of scalar there is no polarization tensor so multiply by identity
//In the case of vector multiple by a polarization 4-vector
//mo support of tensors yet

//**********************************************************************************************************************

map< int, Eigen::MatrixXcd > SubdPol::Subduce_with_pol(double& mom_sq, double& mass_sq,  int& twoJ, const irrep_label& irrep,
                                   const string& little_group,double R1_phi, double R1_theta, double R1_psi){
    
    /* Variables */
    map< int, complex<double> > Sub;
    map< int, Eigen::MatrixXcd > Sub_with_pol;
    typedef std::complex<double> cd;
    

    switch(int(mom_sq)){

        //At rest subductions    
        case 0:{

            Sub = Subd::subduce_oct(irrep);

            if(twoJ == 2){
                for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    Sub_with_pol.insert(make_pair( (it->first) , (it->second) * KfUt::Gmunu() * PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi)));}}

            else if(twoJ == 0){
                Eigen::MatrixXcd unit_vec(4,1);
                unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
                for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    
                    Sub_with_pol.insert(make_pair( (it->first) , (it->second) * unit_vec ));
                    
                }}

            else{
                    cerr << "Tensors and fermions not coded" << endl; exit(1);
                }

            break;
            }
            
        //In flight    
        default:{

            if(twoJ%2){cerr << "Fermions not coded" << endl; exit(1);} //fermions
            
            else{
                Sub = Subd::subduce_lg_boson(irrep, little_group);

                if(twoJ == 2){
                    for(map< int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){

                        Sub_with_pol.insert(make_pair( (it->first) , (it->second) * KfUt::Gmunu() * PolVec::getPol4(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi) ));}}

                else if(twoJ == 0){
                    Eigen::MatrixXcd unit_vec(4,1);
                    unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
                    for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                        Sub_with_pol.insert(make_pair( (it->first) ,(it->second) * unit_vec )) ;}}


                else{cerr << "Tensors not coded" << endl; exit(1);}     

            }  //bosons
            break;
            
            }
            
            
    }
    return Sub_with_pol;
    
};
