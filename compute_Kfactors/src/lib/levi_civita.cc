/* Antisymmetric Tensor */

#include "levi_civita.h"

  //**********************************************************************************************************************

double LevCiv::LeviCivita(int arr[], int n){ //modifies the input and sorts it
  
  int cnt =0;
  int i, j;
    
    for (i = 0; i < n-1; i++){
      for (j = 0; j < n-i-1; j++){
        
        if (arr[j] == arr[j+1] ){return 0;}
        else if (arr[j] > arr[j+1]){swap(arr[j], arr[j+1]); ++cnt;}}}; //gives the result based on the number of swaps
 

   //cout << "lev" << pow(-1,cnt) << "lev";  
   return pow(-1,cnt);
 

};

  //**********************************************************************************************************************
