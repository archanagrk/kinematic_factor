
#include "subduction.h"


/* Register all  */

bool linkageHack(void) __attribute__ ((constructor));

bool linkageHack(void){

  bool foo = true;

  foo &= Hadron::SubduceTableEnv::registerAll();
  foo &= Hadron::CubicGroupsEnv::registerAll();
  foo &= Hadron::IrrepsCubicEnv::registerAll();

  return foo;
};

/* Find the embeddings Jo's code */


int Subd::find_n_subduced_embeddings(const string& group, const string& irrep, int twoJ, int eta_tilde){
  int embeds = 1;  //1 based labelling of embeddings

  if(group == "Oh"){ // at rest
    bool finished = false;
    while(!finished){
      stringstream tmp; tmp << "J" << J_name(twoJ) << "->" << irrep << "," << embeds;
      string irrep_label = tmp.str();
      if( Hadron::TheSubduceTableFactory::Instance().exist(irrep_label) ){ embeds++; }
      else{ finished = true; }
    }
    return embeds -1;
  }

  // in flight
  if((twoJ%2 == 0)){
    //bosons
    for(int two_abs_lam = 0; two_abs_lam <= twoJ; two_abs_lam += 2){
      stringstream tmp; tmp << "H" << J_name(two_abs_lam);
      if(two_abs_lam == 0){ tmp << sign(eta_tilde); }

      tmp << "->H" << J_name(two_abs_lam) << group << irrep << ",1"; //always only ever one embedding of each |lambda| (for bosons)
      string irrep_label = tmp.str();

      //cout << "looking for irrep_label = " << irrep_label << endl;

      if(Hadron::TheSubduceTableFactory::Instance().exist(irrep_label)){ embeds++; }
    }//next abs_lam
  }
  else{
    //fermions
    for(int two_abs_lam = 1; two_abs_lam <= twoJ; two_abs_lam += 2){
      stringstream tmp; tmp << "H" << J_name(two_abs_lam);
      tmp << "->H" << J_name(two_abs_lam) << group << irrep << ",1"; //always only ever one embedding of each |lambda|
      string irrep_label = tmp.str();

      //cout << "looking for irrep_label = " << irrep_label << endl;

      if(Hadron::TheSubduceTableFactory::Instance().exist(irrep_label)){ embeds++; }

    }
  }

  return embeds - 1;  //1 based labelling
};





/* Jo's Subduction code - in helicity basis -  Representations of Irreps taken from Appendix E - Table VII - X - Helicity Ops for Mesons */


//in flight subductions - for bosons
map< int, complex<double> > Subd::subduce_lg_boson(const irrep_label& irrep, const string& little_group){

  map< int, complex<double> > out;
  int eta_tilde = irrep.P;
  if (((irrep.twoJ/2)%2) ) { eta_tilde *= -1; }

  //cout << "eta_tilde = " << eta_tilde << endl;

  // this is the eta_tilde flag = P(-1)^J
  // P = Parity

  // names get swapped when eta_tilde flips
  // if they aren't all the irrep names will be swapped (A_1 <-> A_2, B_1 <-> B_2, E_2[row1] <-> E_2[row_2])

  const string substr = little_group + irrep.irrep + ",1";//only ever one embedding for each |lambda|

    for(int abs_lam = 0; abs_lam <= irrep.twoJ/2; abs_lam++){
      //get the subduction table
      int count_embedding = 0;
      stringstream tmp; tmp << "H" << abs_lam;
      if(abs_lam == 0){ tmp << sign(eta_tilde); }
      tmp << "->H" << abs_lam << substr;
      string label = tmp.str();


      if( !(Hadron::TheSubduceTableFactory::Instance().exist(label) )){ continue; } // weeds out unwanted abs_lam(only +- piece remains)
      count_embedding++; if( irrep.n != count_embedding ){ continue; } // wrong embedding

      //**************************************************************************
      // cout << "count_embedding = " << count_embedding << endl;

      ADAT::Handle< Hadron::SubduceTable > V = Hadron::TheSubduceTableFactory::Instance().createObject(label);

      complex<double> sub = (*V)(irrep.row, 1); //these should be real according to Table II of the in-flight paper
      if(abs(sub) > 0.0 ){out.insert(make_pair(2*abs_lam,sub));}  //lam +ve piece

      if(abs_lam != 0){ // always contains a factor of eta_tilde
        sub = (*V)(irrep.row, 2)* double(eta_tilde);
        if(abs(sub) > 0.0 ){out.insert(make_pair(-2*abs_lam,sub));} // lam -ve piece
      }

    }//next abs_lam
    //complex<double> s = complex<double>( toDouble(real(sub)), toDouble(imag(sub)) );

  return out;
};



//in flight subductions - for fermions
// -- uses the phase choices suggested by Robert, taken from Christopher - NOT CHECKED YET
map< int, complex<double> > Subd::subduce_lg_fermion(const irrep_label& irrep, const string& little_group){

  map< int, complex<double> > out;
  int parity = irrep.P;  // parity = P

  double phi = double(irrep.twoJ - 1) * Subd::PI / double(2.0);
  complex<double> phase = double(parity) * complex<double>( cos(phi), sin(phi) ); //cout << "phase = " << phase << endl;


    for(int two_abs_lam = 1; two_abs_lam <= irrep.twoJ; two_abs_lam += 2){
      int count_embedding = 0;
      //get the subduction table
      stringstream tmp; tmp << "H" << J_name(two_abs_lam)
                            << "->H" << J_name(two_abs_lam) << little_group << irrep.irrep << ",1"; //only ever one embedding for each |lambda|
      string label = tmp.str(); //cout << "label = " << label << endl;


      // this appears to order the embeddings in a sensible way
      if( !(Hadron::TheSubduceTableFactory::Instance().exist(label) )){ continue; } // don't want this subduction
      count_embedding++; if( irrep.n != count_embedding ){ continue; } // wrong embedding

      //**************************************************************************

      // cout << "count_embedding = " << count_embedding << endl;

      ADAT::Handle< Hadron::SubduceTable > V = Hadron::TheSubduceTableFactory::Instance().createObject(label);

      complex<double> sub = (*V)(irrep.row, 1); //these should be real according to Table II of the in-flight paper
       //lam +ve piece
      if(abs(sub) > 0.0 ){out.insert(make_pair(two_abs_lam,sub));}

      //lam -ve piece
      sub = (*V)(irrep.row, 2) * phase;                       //v = real( (*V)(irrep.row, 2) );
      if(abs(sub) > 0.0 ){out.insert(make_pair(-two_abs_lam,sub));}


    }//next abs_lam

    //    complex<double> s = complex<double>( toDouble(real(sub)), toDouble(imag(sub)) );

  return out;
};


// at rest subductions
map< int, complex<double> > Subd::subduce_oct(const irrep_label& irrep){

  map< int, complex<double> > out;

  //get the subduction table
  stringstream tmp; tmp << "J" << J_name(irrep.twoJ) << "->" << irrep.irrep << "," << irrep.n;
  string label = tmp.str();   //cout << "label = " << label << endl;;

  if( Hadron::TheSubduceTableFactory::Instance().exist(label) ){
    ADAT::Handle< Hadron::SubduceTable > V = Hadron::TheSubduceTableFactory::Instance().createObject(label);

    for(int two_lam = -irrep.twoJ; two_lam <= irrep.twoJ; two_lam += 2){
      complex<double> s = (*V)(irrep.row, irrep.twoJ/2 - two_lam/2 + 1) ;   // some strange indexing choice in ADAT-- should be real
      //complex<double> s = complex<double>( toDouble(real(v)), toDouble(imag(v)) );
      if(abs(s) > 0.0){out.insert(make_pair(two_lam,s));}  // lam = mu at rest
    }
  }
  return out;
};




/* Perform the subductions by adding over the J_z and abs_lams  */

sub_hel Subd::Subduce_all(double& mom_sq, double& mass_sq,  int& twoJ, const irrep_label& irrep,
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
            sum +=  (it->second) * PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi); two_hel = it->first;}}
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
              sum +=  (it->second) * PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi);  two_hel = it->first;}}
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

map< int, Eigen::MatrixXcd > Subd::Subduce_with_phases(double& mom_sq, double& mass_sq,  int& twoJ, const irrep_label& irrep,
                                   const string& little_group,double R1_phi, double R1_theta, double R1_psi){
    
    
    map< int, complex<double> > Sub;
    map< int, Eigen::MatrixXcd > Sub_with_pol;
    typedef std::complex<double> cd;
    
    switch(int(mom_sq)){
            
        case 0:{
            Sub = Subd::subduce_oct(irrep);
            if(twoJ == 2){
                for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    //cout << "Polstart" << PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi)<< "Polend    "; cout << irrep.irrep << "\n";
                    Sub_with_pol.insert(make_pair( (it->first) , (it->second) * PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi)));}}
            else{
                Eigen::MatrixXcd unit_vec(4,1);
                unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
                for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    
                    Sub_with_pol.insert(make_pair( (it->first) , (it->second) * unit_vec ));
                    
                }}
            break;}
            
            
        default:{
            if(twoJ%2){
                Eigen::MatrixXcd unit_vec(4,1);
                unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
                Sub  = Subd::subduce_lg_fermion(irrep,little_group);
                for(map< int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                    Sub_with_pol.insert(make_pair( (it->first) ,  (it->second) * unit_vec ));}
            } //fermions
            
            else{
                Sub = Subd::subduce_lg_boson(irrep, little_group);
                if(twoJ == 2){
                    for(map< int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                        //cout << "Polstart" << PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi)<< "Polend"; cout << irrep.irrep << "\n";
                        Sub_with_pol.insert(make_pair( (it->first) , (it->second) * PolVec::getPolarization(mom_sq, (it->first), mass_sq, R1_phi, R1_theta, R1_psi) ));}}
                else{
                    Eigen::MatrixXcd unit_vec(4,1);
                    unit_vec << cd(1,0),cd(1,0),cd(1,0),cd(1,0);
                    for(map<  int, complex<double> >::iterator  it = Sub.begin(); it != Sub.end(); it++){
                        Sub_with_pol.insert(make_pair( (it->first) ,(it->second) * unit_vec )) ;}}
            }  //bosons
            break;}
            
            
    }
    return Sub_with_pol;
    
};
