//
//  main.cpp
//  Subduction
//
//  Created by Archana Radhakrishnan on 7/11/18.
//  Copyright © 2018 Archana Radhakrishnan. All rights reserved.
//

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



/* Find Little Groups - Jo's code  */

string generateLG(ArrayXXd& mom)
{
  //ArrayXXd mom(1,3);

  std::string littleGroup = "";

  if (mom(0,2) == 0){
    if (mom(0,1) == 0){
      if (mom(0,0) == 0)   // 0 0 0
        littleGroup = "Oh";

      else                  // n 0 0
        littleGroup = "D4";
    }

    else{
      if (mom(0,0) == mom(0,1))   // n n 0
        littleGroup = "D2";

      else
        littleGroup = "C4";         // n m 0
    }
  }

  else{
    if ( (mom(0,0) == mom(0,1)) && (mom(0,0) == mom(0,2)) )   // n n n
      littleGroup = "D3";

    else if ( (mom(0,0) == mom(0,1)) || (mom(0,0) == mom(0,2)) || (mom(0,1) == mom(0,2)) )   // m m n
      littleGroup = "C4";

    else   // n m p
      littleGroup = "C2";
  }

  return littleGroup;
};



/* Get the irrrep */

std::vector<std::string>  getIrrep(int& twoJ, int& P, string& lg)
{ 

 std::vector<std::string> Irrep;				// Irreps for twoJ=0 and Parity = -1 
 if (lg == "Oh" && twoJ == 0 && P == -1){
	Irrep.push_back("A1");
	}

 else if(lg != "Oh" && twoJ == 0 && P == -1){
 	Irrep.push_back("A2");  
	}


 else if(lg == "Oh" && twoJ == 2){				// Irreps for J=1
         Irrep.push_back("T1");
	}

 else if(lg == "D4" && twoJ == 2){
        Irrep.push_back("E2");
        }

 else if(lg == "D2" && twoJ == 2){
        Irrep.push_back("B1"); Irrep.push_back("B2");
        }

 else if(lg == "D3" && twoJ == 2){
	Irrep.push_back("E2");
        }

 else if(lg == "C4" && twoJ == 2){
 Irrep.push_back("A1"); Irrep.push_back("A2");
        }

 else{
	cout << "Irrep not in database";
	exit(1);
	}

 return Irrep;
};


/* Get the max rows of each irrep */

int irrepRows(string& irrep)
{
 int rows;
 switch(irrep[0])
 {
   case 'A':
	    rows = 1;			// A,A1,A2 irrep
   	    break;
   case 'B':
            rows = 1;			// B,B1,B2 irrep
            break;
   case 'E':
            rows = 2;			// E,E1,E2,E3 irrep
            break;
   case 'T':
            rows = 3;			// T1,T2 irrep
            break;
   case 'G':
            rows = 2;			// G1,G2 for fermions
            break;
   case 'H':
            rows = 4;			// H irrep for fermions
            break;
 }

 return rows;

};




/* Jo's code - Euler matrix - converted from itpp to Eigen3 - z-y-z convention */

MatrixXd eulerRotMat(double alpha, double beta, double gamma)
{
  // to act upon cartesian 4-vectors
  // R(a,b,g) = Rz(a) Ry(b) Rz(g)

  MatrixXd Rzg = MatrixXd::Zero(4,4);
  Rzg(1,1) = cos(gamma); Rzg(1,2) = -1 * (sin(gamma));
  Rzg(2,1) = -1 * (Rzg(1,2)); Rzg(2,2) = Rzg(1,1);
  Rzg(3,3) = 1; Rzg(0,0) = 1;

  MatrixXd Ryb = MatrixXd::Zero(4,4);
  Ryb(1,1) = cos(beta); Ryb(1,3) = sin(beta);
  Ryb(2,2) = 1; Ryb(0,0) = 1;
  Ryb(3,1) = -1 * (Ryb(1,3)); Ryb(3,3) = Ryb(1,1);

  MatrixXd Rza = MatrixXd::Zero(4,4);
  Rza(1,1) = cos(alpha); Rza(1,2) = -1 * (sin(alpha));
  Rza(2,1) = -1 * (Rza(1,2)); Rza(2,2) = Rza(1,1);
  Rza(3,3) = 1; Rza(0,0) = 1;

  return Rza * Ryb * Rzg;
};




/* To compute the polarization given the momentum and helicity - Convention from Appendix B of Helicity Ops for Mesons  */

 MatrixXcd getPolarization(double& abs_mom,int& helicity, double& mass, double& phi, double& theta, double& psi)
 {
                                                                 // For spin 1 particles extra polarization degree of freedom
   MatrixXcd pol_z(1,4);
   MatrixXcd pol(1,4);

   if(helicity == 0){
     if(mass == 0){pol_z << (0,0),(0,0),(0,0),(0,0);}            // Photon has only two physical polarizations

     else{
       double energy = sqrt(pow(abs_mom,2)+pow(mass,2));
       pol_z << (abs_mom/mass,0),(0,0),(0,0),(energy/mass,0);   // k.pol = 0
      }
   }

   else if(helicity == 1){pol_z <<  (0,0),(-sqrt(0.5),0),(0,-sqrt(0.5)),(0,0);}        // k.pol = 0
   else if(helicity == -1){pol_z << (0,0),(sqrt(0.5),0),(0,-sqrt(0.5)),(0,0);}

   pol = eulerRotMat(phi,theta,psi)*pol_z;                              // multiplies by the euler matrix to convert p_ref to p_canonical
   return pol;
 };






/* Get the reference angles for each LG - copied from Jo - Appendix E Table VI - Helicity ops for Mesons - z-y-z  John-Wick convention */


std::vector<float> refAngles(string little_group){
  std::vector<float> ref;
  double R1_phi, R1_theta, R1_psi; //reference rotation angles

  if(little_group == "Oh"){ R1_phi = 0.0; R1_theta = 0.0; R1_psi = 0.0;}
  else if(little_group == "D4"){ R1_phi = 0.0; R1_theta = 0.0; R1_psi = 0.0;} // (00n)
  else if(little_group == "D2"){ R1_phi = PI/2.0; R1_theta = PI/4.0; R1_psi = -PI/2.0;}// (0nn)
  else if(little_group == "D3"){ R1_phi = PI/4.0; R1_theta = acos(1.0/sqrt(3.0)); R1_psi = 0.0;}// (nnn)
  else if(little_group == "C4"){ cerr << "C4 not coded" << endl; exit(1); }// ???????
  else if(little_group == "C2"){ cerr << "C2 not coded" << endl; exit(1); }// ???????
  ref.push_back(R1_phi); ref.push_back(R1_theta); ref.push_back(R1_psi);

  return ref;
  };




/* J's for string stream Jo's code */



string J_name( int twoJ ){
  stringstream ss;

  if( twoJ % 2 ){ //fermion
    ss << twoJ << "o2";
  }

  else{ //boson
    ss << (twoJ/2);
  }

  return ss.str();
}

/* For sign of eta tilde  stolen from Jo  */


string sign(int x){

  if(x == -1){ return "-"; }
  else{ return "+"; }

}



/* Find the embeddings Jo's code */


int find_n_subduced_embeddings(const string& group, const string& irrep, int twoJ, int eta_tilde){
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




/* Jo's Subduction code - Representations of Irreps taken from Appendix E - Table VII - X - Helicity Ops for Mesons */




//in flight subductions - for bosons
map< int, complex<double> > subduce_lg_boson(const irrep_label& irrep, const string& little_group,
                                             double R1_theta, double R1_phi, double R1_psi){

  map< int, complex<double> > out;
  int eta_tilde = irrep.P1P2;
  if (((irrep.twoJ/2 + irrep.ell)%2) ) { eta_tilde *= -1; }

  //cout << "eta_tilde = " << eta_tilde << endl;

  // this is the eta_tilde flag = P(-1)^J
  // P = P_1 * P_2 * (-1)^ell
  // so eta_tilde = P_1 * P_2 * (-1)^(ell+J)

  // names get swapped when eta_tilde flips
  // if they aren't all the irrep names will be swapped (A_1 <-> A_2, B_1 <-> B_2, E_2[row1] <-> E_2[row_2])

  //djw: stringstreams are slow, this bit never changes
  const string substr = little_group + irrep.irrep + ",1";//only ever one embedding for each |lambda|

  for(int mu = -irrep.twoJ/2; mu <= irrep.twoJ/2; mu++){
    int count_embedding = 0;
    complex<double> sub(0., 0.);

    for(int abs_lam = 0; abs_lam <= irrep.twoJ/2; abs_lam++){
      //get the subduction table
      stringstream tmp; tmp << "H" << abs_lam;
      if(abs_lam == 0){ tmp << sign(eta_tilde); }
      tmp << "->H" << abs_lam << substr;
      string label = tmp.str();

      //    cout << "subduce_lg_boson::label = " << label << endl;
      //    cout << "irrep.n = " << irrep.n << endl;

      // this appears to order the embeddings in a sensible way
      if( !(Hadron::TheSubduceTableFactory::Instance().exist(label) )){ continue; } // don't want this subduction
      count_embedding++; if( irrep.n != count_embedding ){ continue; } // wrong embedding

      //**************************************************************************
      // cout << "count_embedding = " << count_embedding << endl;

      ADAT::Handle< Hadron::SubduceTable > V = Hadron::TheSubduceTableFactory::Instance().createObject(label);
      complex<double> DR1 = Hadron::Wigner_D(irrep.twoJ, 2*mu, 2*abs_lam, R1_phi, R1_theta, R1_psi);

      complex<double> v = (*V)(irrep.row, 1); //these should be real according to Table II of the in-flight paper

      sub += v * DR1; //lam +ve piece

      if(abs_lam != 0){ // always contains a factor of eta_tilde
        DR1 = Hadron::Wigner_D(irrep.twoJ, 2*mu, -2*abs_lam, R1_phi, R1_theta, R1_psi);
        v = (*V)(irrep.row, 2);
        sub += v * DR1 * double(eta_tilde);
      }

    }//next abs_lam

    //complex<double> s = complex<double>( toDouble(real(sub)), toDouble(imag(sub)) );

    if(abs(sub) > 0.0 ){ out.insert( pair<int, complex<double> >(2*mu, sub) ); };
  }//next mu and mu is the irrep row and J_z

  return out;
};





//in flight subductions - for fermions
// -- uses the phase choices suggested by Robert, taken from Christopher - NOT CHECKED YET
map< int, complex<double> > subduce_lg_fermion(const irrep_label& irrep, const string& little_group,
                                               double R1_phi, double R1_theta, double R1_psi){

  map< int, complex<double> > out;
  int parity = irrep.P1P2; if (!((irrep.ell)%2)) { parity *= -1; } // parity = P_1 * P_2 * (-1)^ell

  double phi = double(irrep.twoJ - 1) * PI / double(2.0);
  complex<double> phase = double(parity) * complex<double>( cos(phi), sin(phi) ); //cout << "phase = " << phase << endl;

  for(int twoMu = -irrep.twoJ; twoMu <= irrep.twoJ; twoMu += 2){
    int count_embedding = 0;
    complex<double> sub(0., 0.);

    for(int two_abs_lam = 1; two_abs_lam <= irrep.twoJ; two_abs_lam += 2){
      //get the subduction table
      stringstream tmp; tmp << "H" << J_name(two_abs_lam)
                            << "->H" << J_name(two_abs_lam) << little_group << irrep.irrep << ",1"; //only ever one embedding for each |lambda|
      string label = tmp.str(); //cout << "label = " << label << endl;

      //    cout << "subduce_lg_boson::label = " << label << endl;
      //    cout << "irrep.n = " << irrep.n << endl;

      // this appears to order the embeddings in a sensible way
      if( !(Hadron::TheSubduceTableFactory::Instance().exist(label) )){ continue; } // don't want this subduction
      count_embedding++; if( irrep.n != count_embedding ){ continue; } // wrong embedding

      //**************************************************************************

      // cout << "count_embedding = " << count_embedding << endl;

      ADAT::Handle< Hadron::SubduceTable > V = Hadron::TheSubduceTableFactory::Instance().createObject(label);
      complex<double> DR1 = Hadron::Wigner_D(irrep.twoJ, twoMu, two_abs_lam, R1_phi, R1_theta, R1_psi);

      complex<double> v = (*V)(irrep.row, 1); //these should be real according to Table II of the in-flight paper
      sub += v * DR1; //lam +ve piece

      //lam -ve piece
      DR1 = Hadron::Wigner_D(irrep.twoJ, twoMu, -two_abs_lam, R1_phi, R1_theta, R1_psi);
      v = (*V)(irrep.row, 2);                       //v = real( (*V)(irrep.row, 2) );
      sub += v * DR1 * phase;


    }//next abs_lam

    //    complex<double> s = complex<double>( toDouble(real(sub)), toDouble(imag(sub)) );
    if(abs(sub) > 0.0 ){ out.insert( pair<int, complex<double> >(twoMu, sub) ); };
  }//next mu

  return out;
};


// at rest subductions
map< int, complex<double> > subduce_oct(const irrep_label& irrep){

  map< int, complex<double> > out;

  //get the subduction table
  stringstream tmp; tmp << "J" << J_name(irrep.twoJ) << "->" << irrep.irrep << "," << irrep.n;
  string label = tmp.str();   //cout << "label = " << label << endl;;

  if( Hadron::TheSubduceTableFactory::Instance().exist(label) ){
    ADAT::Handle< Hadron::SubduceTable > V = Hadron::TheSubduceTableFactory::Instance().createObject(label);

    for(int two_m = -irrep.twoJ; two_m <= irrep.twoJ; two_m += 2){
      complex<double> s = (*V)(irrep.row, irrep.twoJ/2 - two_m/2 + 1) ;   // some strange indexing choice in ADAT-- should be real
      //complex<double> s = complex<double>( toDouble(real(v)), toDouble(imag(v)) );
      if(abs(s) > 0.0){ out.insert( pair<int, complex<double> >(two_m, s) ); }
    }
  }
  return out;
};



/* Automate the subductions for all types of particles */

map< int, complex<double> > Subduce_all(float& mom_sq, int& twoJ, const irrep_label& irrep, const string& little_group,
                                             double R1_phi, double R1_theta, double R1_psi){

  map< int, complex<double> > out;

  if(mom_sq == 0){out = subduce_oct(irrep);}
  else if(!(twoJ%2)){out  = subduce_lg_fermion(irrep,little_group,R1_phi,R1_theta,R1_psi); } //fermions

       else{out =  subduce_lg_boson(irrep, little_group,R1_phi,R1_theta,R1_psi);} //bosons

  return out;

};

  


/* Antisymmetric Tensor */

double LeviCivita(int ei,int ej, int ek, int el){

    if(ei == ej || ej == ek || ek == el || el == ei || ei == ek || ej == el){
       return 0;
       }

    int cnt = 0;
    std::vector<int> e = {ei,ej,ek,el};
    std::vector<int> per = {0,1,2,3};

    for (int i = 0; i < 4; ++i){
        while (i != e[i]){
            ++cnt;
            std::swap (per[i], per[e[i]]);
            std::swap (e[i], e[e[i]]);
	}
    }

    return pow(-1, cnt);
};




/* Loop over Lorentz indices to get the coefficient*/


complex<double> KinematicFactor(ArrayXXd& qp, ArrayXXd& qm, map< int, complex<double> >& Sub1 , map< int, complex<double> > SubCurr , map< int, complex<double> > Sub3 ){

    complex<double> Coeff = 0;
    auto sub_1 = Sub1.begin();
    auto sub_curr = SubCurr.begin();
    auto sub_3 = Sub3.begin();

    for(int i = 0; i < 4; i++ ){for(int j = 0; j < 4; j++ ){for(int k = 0; k < 4; k++ ){for(int l = 0; l< 4; l++ ){

      Coeff=Coeff+(qp(0,j))*(qm(0,k))*(sub_1->second)*(sub_curr->second)*(sub_3->second)*LeviCivita(i,j,k,l);}}}}

      return Coeff;
};
     
  


/* main code to return the irreps to compute */

int main(int argc, char** argv){
  if( argc != 14 ){
    cerr << "get_irreps <twoJ1> <P1> <m1sq> <twoJ2> <P2> <m2sq> <twoJ3> <P3> <m3sq> <max_mom1> <max_mom2> <max_mom3> <wave> \n ";
    exit(1); }

  int J1;  {istringstream a(argv[1]); a >> J1;};
  int P1;   {istringstream a(argv[2]); a >> P1;};
  int m1;   {istringstream a(argv[3]); a >> m1;};
  int J2;  {istringstream a(argv[4]); a >> J2;};
  int P2;   {istringstream a(argv[5]); a >> P2;};
  int m2;   {istringstream a(argv[6]); a >> m2;};
  int J3;  {istringstream a(argv[7]); a >> J3;};
  int P3;   {istringstream a(argv[8]); a >> P3;};
  int m3;   {istringstream a(argv[9]); a >> m3;};
  int max_mom1;   {istringstream a(argv[10]); a >> max_mom1;};
  int max_mom2;   {istringstream a(argv[11]); a >> max_mom2;};
  int max_mom3;   {istringstream a(argv[12]); a >> max_mom3;};
  int wave;	{istringstream a(argv[13]); a >> wave;};
  
  complex<double> Coeff;
  
  if((J1==0) && (P1==-1) && (J2==2) && (P2==-1) && (J3==2) && (P3==-1)){
    ArrayXXd  mom1(1,3);
    ArrayXXd mom3(1,3);
    ArrayXXd mom_curr(1,3);
    
    for(int i=0; i < sqrt(max_mom1); i++){
       for(int j = 0; j <= i; j++){						// Looping over all mom_1
	 for(int k = 0; k <= j; k++){


           float mom1_sq = (pow(i,2)+pow(j,2)+pow(k,2));

           while(max_mom1 >= mom1_sq){	

	     for(int l = 0; l < sqrt(max_mom3); l++){
                for(int m=0; m <= l ; m++){					// Looping over all mom_3
           	   for(int n=0; n<=m; n++){

		   float mom3_sq = (pow(l,2)+pow(m,2)+pow(n,2));

		   while(max_mom3 >= mom3_sq){

		    mom1 << i,j,k;
		    mom3 << l,m,n;
		    mom_curr << (i-l),(j-m),(k-n);

		    float mom_curr_sq = (pow((i-l),2)+pow((j-m),2)+pow((k-n),2));

		    cout << mom1 << "\n";
        	    cout << mom3 << "\n";

		    ArrayXXd  qp(1,4);  ArrayXXd  qm(1,4);  // q+ = (p1-p2) q- = (p1+p2)

		    qp  << (sqrt(m1+mom1_sq)+sqrt(m3+mom3_sq)),(i+l),(j+m),(k+n);
		    qm  << (sqrt(m1+mom1_sq)-sqrt(m3+mom3_sq)),(i-l),(j-m),(k-n);

  		    string LG1 = generateLG(mom1);
		    string LG3 = generateLG(mom3);
		    string LG_curr = generateLG(mom_curr);

		    std::vector<float> r1 = refAngles(LG1);
		    std::vector<float> r_curr = refAngles(LG_curr);
		    std::vector<float> r3 = refAngles(LG3);
		    std::vector<std::string> irrep1 = getIrrep(J1,P1,LG1);

		    //cout << *(irrep1.begin());
	       	    std::vector<std::string> irrep_curr = getIrrep(J2,P2,LG_curr);
                    std::vector<std::string> irrep3 = getIrrep(J3,P3,LG3);
                    irrep_label rep1; irrep_label rep_curr; irrep_label rep3;

		    rep1.ell = rep_curr.ell = 0;
		    rep3.ell = wave;

		    rep1.twoS = J1; rep3.twoS = J3; rep_curr.twoS = J2;
		    rep1.twoJ = J1; rep3.twoJ = J3; rep_curr.twoJ = J2;
		    rep1.P1P2 = P1; rep3.P1P2 = P3; rep_curr.P1P2 = P2;


		    for(auto p = irrep1.begin(); p != irrep1.end(); p++){
			for(auto q = irrep3.begin(); q != irrep3.end(); q++){         // Looping over all irreps at the source, sink & isertion
			   for(auto r = irrep_curr.begin(); r != irrep_curr.end(); r++){

				rep1.irrep = *p; rep3.irrep = *q; rep_curr.irrep = *r;
				cout << *p << *q << *r << "\n";
				rep1.n = find_n_subduced_embeddings(LG1, rep1.irrep, rep1.twoJ, (rep1.P1P2*pow(-1,(J1/2))));
                                rep_curr.n = find_n_subduced_embeddings(LG_curr, rep_curr.irrep, rep_curr.twoJ, (rep_curr.P1P2*pow(-1,(J2/2))));
                                rep3.n = find_n_subduced_embeddings(LG3, rep3.irrep, rep3.twoJ, (rep3.P1P2*pow(-1,rep3.ell + (J3/2))));


				for(int row1 = 1; row1 <= irrepRows(*p); row1++){     // Looping over all irrep rows at source, sink & insertion
				  for(int row_curr = 1; row_curr <= irrepRows(*r); row_curr++){
				    for(int row3 = 1; row3 <= irrepRows(*q); row3++){      


					rep1.row = row1; rep3.row = row3; rep_curr.row = row_curr;
					cout << row1 << row3 << row_curr << "\n";
					map< int, complex<double> > Sub1 = Subduce_all(mom1_sq, J1 , rep1, LG1, r1[0], r1[1], r1[2]);
					map< int, complex<double> > Sub3 = Subduce_all(mom3_sq, J3 , rep3, LG3, r3[0], r3[1], r3[2]);
					map< int, complex<double> > SubCurr = Subduce_all(mom_curr_sq, J2 , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2]);


					Coeff = KinematicFactor(qp,qm,Sub1,SubCurr,Sub3);
					cout << "This is" << Coeff << "\n";
		      }}}
		   }}}
		 }	
	     }}}
  	   }
	 }}}
  }
};  
