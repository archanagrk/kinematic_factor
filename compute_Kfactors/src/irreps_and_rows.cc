/* Get the irrrep */

#include "irreps_and_rows.h"


std::vector<std::string>  IrrepName::getIrrep(int& twoJ, int& P, string& lg)
{

 std::vector<std::string> Irrep;                                // Irreps for twoJ=0 and Parity = -1
 if (lg == "Oh" && twoJ == 0 && P == -1){
        Irrep.push_back("A1");
        }

 else if(lg != "Oh" && twoJ == 0 && P == -1){
        Irrep.push_back("A2");
        }


 else if(lg == "Oh" && twoJ == 2){                              // Irreps for J=1
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

int IrrepName::irrepRows(string& irrep)
 {
 int rows;
 switch(irrep[0])
 {
   case 'A':
            rows = 1;                   // A,A1,A2 irrep
            break;
   case 'B':
            rows = 1;                   // B,B1,B2 irrep
            break;
   case 'E':
            rows = 2;                   // E,E1,E2,E3 irrep
            break;
   case 'T':
            rows = 3;                   // T1,T2 irrep
            break;
   case 'G':
            rows = 2;                   // G1,G2 for fermions
            break;
   case 'H':
            rows = 4;                   // H irrep for fermions
            break;
   default :
            cerr << "Irrep not coded" << endl; exit(1);
 }

 return rows;

 };

