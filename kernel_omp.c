#include <omp.h>
#include<stdio.h>
#include "kernel.h"

void kernel_omp(int *input, int *ref, int64_t rows, int64_t cols, int penalty) {
  /* int i, j; 

#pragma omp parallel for private (j)
for(int i = 1; i < rows; i++){
    for(int j = 1; j <= i; j++){
      int64_t idx = i * cols + j;                                                                                                                                                
      int sub = idx -cols;                                                                                                                                                       
      int64_t idxNW = sub - 1;                                                                                                                                                   
      int64_t idxN = sub;                                                                                                                                                        
      int64_t idxW = idx - 1;                                                                                                                                                    
      int r = ref[idx];                                                                                                                                                          
      int inputNW = input[idxNW];                                                                                                                                                
      int inputW = input[idxW];                                                                                                                                                  
      int inputN = input[idxN];                                                                                                                                                                                                             input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty);   
    }
 }

#pragma omp barrier
#pragma omp parallel for private(i)
  for(int j = cols -1; j > 0; j--){
    for(int i = 2; i < j; i++){
      int64_t idx = i * cols + j;                                                                                                                                                
      int sub = idx -cols;                                                                                                                                                       
      int64_t idxNW = sub - 1;                                                                                                                                                   
      int64_t idxN = sub;                                                                                                                                                        
      int64_t idxW = idx - 1;                                                                                                                                                    
      int r = ref[idx];                                                                                                                                                          
      int inputNW = input[idxNW];                                                                                                                                                
      int inputW = input[idxW];                                                                                                                                                  
      int inputN = input[idxN];                                                                                                                                                  
                                                                                                                                                                                 
      input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty);   
    }
  }


}*/ 

 int TILESIZE = 4;
 int i, k, l;
 for(int j = 1; j < cols-TILESIZE; j += TILESIZE){
    #pragma omp parallel for private (k,l)
    for(i = 1; i <= j; i+= TILESIZE){
	for (k = i; k < i + TILESIZE; k++){
	  for(l = j; l < j + TILESIZE; l++){
	    int64_t idx = k * cols + l;           
	    int sub = idx - cols;
	    int64_t idxNW = sub - 1;                                                  
	    int64_t idxN = sub;                                                       
	    int64_t idxW = idx - 1;                                                          
	    int r = ref[idx];                                                                
	    int inputNW = input[idxNW];                                                      
	    int inputW = input[idxW];                                                        
	    int inputN = input[idxN];                                                         
	    input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty); 
	  }
	}
      }
   }
}


/*#pragma omp parallel for collapse(3)
 for (int i = 1; i < rows; i++) {
    for(int initialI = i; initialI > 0; initialI--){
      for(int j = 1; j <= initialI; j++){
      int64_t idx = i * cols + j;
      int sub = idx -cols;
      int64_t idxNW = sub - 1;
      int64_t idxN = sub;
      int64_t idxW = idx - 1;
      int r = ref[idx];
      int inputNW = input[idxNW];
      int inputW = input[idxW];
      int inputN = input[idxN];

      input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty);
      }
    }
 }
}*/


 /* #pragma parallel for collapse(2)
  for (int i = 1; i < rows; ++i) {
    for (int j = 1; j < cols; ++j) {
      int64_t idx = i * cols + j;
      int64_t temp = idx - cols;
      int64_t idxNW = temp - 1;
      int64_t idxN = temp;
      int64_t idxW = idx - 1;
      int r = ref[idx];
      int inputNW = input[idxNW];
      int inputW = input[idxW];
      int inputN = input[idxN];

      input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty);
    }
    }
 }*/
