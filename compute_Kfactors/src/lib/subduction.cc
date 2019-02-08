
#include "subduction.h"

//**********************************************************************************************************************

/* Register all  */

//**********************************************************************************************************************

//bool linkageHack(void) __attribute__ ((constructor));

bool linkageHack(void){

  bool foo = true;

  foo &= Hadron::SubduceTableEnv::registerAll();
  foo &= Hadron::CubicGroupsEnv::registerAll();
  foo &= Hadron::IrrepsCubicEnv::registerAll();

  return foo;
};


//**********************************************************************************************************************

/* Find the embeddings */

//**********************************************************************************************************************


int Subd::find_n_subduced_embeddings(const string& group, const string& irrep, int twoJ, int eta_tilde){//returns the embedding which is int

  linkageHack();
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


      if(Hadron::TheSubduceTableFactory::Instance().exist(irrep_label)){ embeds++; }
    }//next abs_lam
  }



  else{
    //fermions
    for(int two_abs_lam = 1; two_abs_lam <= twoJ; two_abs_lam += 2){
      stringstream tmp; tmp << "H" << J_name(two_abs_lam);
      tmp << "->H" << J_name(two_abs_lam) << group << irrep << ",1"; //always only ever one embedding of each |lambda|
      string irrep_label = tmp.str();


      if(Hadron::TheSubduceTableFactory::Instance().exist(irrep_label)){ embeds++; }

    }
  }

  return embeds - 1;  //1 based labelling
};



//**********************************************************************************************************************

/* Subduction code - in helicity basis -  Representations of Irreps taken from Appendix E - Table VII - X - Helicity Ops for Mesons */

//**********************************************************************************************************************


//returns a map of helicity and subduction coeff


//in flight subductions - for bosons
map< int, complex<double> > Subd::subduce_lg_boson(const irrep_label& irrep, const string& little_group){


  linkageHack();
  map< int, complex<double> > out;
  int eta_tilde = irrep.P;


  if (((irrep.twoJ/2)%2) ) { eta_tilde *= -1; }


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

//**********************************************************************************************************************

//in flight subductions - for fermions
// -- uses the phase choices suggested by Robert, taken from Christopher - NOT CHECKED YET

map< int, complex<double> > Subd::subduce_lg_fermion(const irrep_label& irrep, const string& little_group){

  //returns a map of helicity and subduction coeff

  linkageHack();
  map< int, complex<double> > out;
  int parity = irrep.P;  // parity = P

  double phi = double(irrep.twoJ - 1) * PI / double(2.0);
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

//**********************************************************************************************************************

// at rest subductions
map< int, complex<double> > Subd::subduce_oct(const irrep_label& irrep){

  linkageHack();
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

//**********************************************************************************************************************




