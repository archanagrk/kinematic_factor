/* Antisymmetric Tensor */

#include "levi_civita.h"

  //**********************************************************************************************************************

double LevCiv::LeviCivita(int arr[], int n){ //modifies the input and sorts it
  
  int cnt =0;
  int i, j;
  int arr_copy[n];

  for (int r = n-1; r >= 0; r--) {arr_copy[r] = arr[r];}
    
    for (i = 0; i < n-1; i++){
      for (j = 0; j < n-i-1; j++){
        
        if (arr_copy[j] == arr_copy[j+1] ){return 0;}
        else if (arr_copy[j] > arr_copy[j+1]){swap(arr_copy[j], arr_copy[j+1]); ++cnt;}}}; //gives the result based on the number of swaps
 
   return pow(-1,cnt);
 

};

  //**********************************************************************************************************************
