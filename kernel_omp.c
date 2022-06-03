#include <omp.h>
#include<stdio.h>
#include "kernel.h"

void kernel_omp(int *input, int *ref, int64_t rows, int64_t cols, int penalty) {
  /*  int i, j, k, l; 

  //#pragma omp parallel for private (k)
  for(i = 1, k = 1; i < rows, k <=1; i++, k++){
    // for(k  = 1; k <= i; k++){
      int64_t idx = i * cols + k;                                                                                                                                                
      int sub = idx -cols;                                                                                                                                                       
      int64_t idxNW = sub - 1;                                                                                                                                                   
      int64_t idxN = sub;                                                                                                                                                        
      int64_t idxW = idx - 1;                                                                                                                                                    
      int r = ref[idx];                                                                                                                                                          
      int inputNW = input[idxNW];                                                                                                                                                
      int inputW = input[idxW];                                                                                                                                                  
      int inputN = input[idxN];                                                                                                                                                                                                             input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty);   
      // }
 }

//#pragma omp barrier
//#pragma omp parallel for private(l)
  for(j = cols, l = 2; j > 0, l < j; j--, l++){
    // for(l = 2; l < j; l++){
      int64_t idx = l * cols + j;                                                                                                                                                
      int sub = idx -cols;                                                                                                                                                       
      int64_t idxNW = sub - 1;                                                                                                                                                   
      int64_t idxN = sub;                                                                                                                                                        
      int64_t idxW = idx - 1;                                                                                                                                                    
      int r = ref[idx];                                                                                                                                                          
      int inputNW = input[idxNW];                                                                                                                                                
      int inputW = input[idxW];                                                                                                                                                  
      int inputN = input[idxN];                                                                                                                                                  
                                                                                                                                                                                 
      input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty);   
      // }
  }


  } 
  */
 
 int TILESIZE = 2;
 int i, k, l;
#pragma omp parallel
 for(int j = 1; j < cols-TILESIZE; j += TILESIZE){
#pragma omp for private (i,l)
    for(i = j; i >= 1; i-= TILESIZE){
      printf("i: %d j: %d\n ", i, j);
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

 int x,y,z;
#pragma omp parallel
 for(int j = cols - TILESIZE; j > 0; j -= TILESIZE){
#pragma omp for private (y,z)
   for(x = TILESIZE; x <= rows-TILESIZE; x+= TILESIZE){
     for (y = x; y < x + TILESIZE; y++){
       for(z = j; z < j + TILESIZE; z++){
	 int64_t idx = y * cols + z;
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

  /* for (int i = 1; i < rows; i++) {
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

/*
 #pragma parallel for
  for (int i = 1; i < rows; ++i) {
   #pragma omp for
    for (int j = i; j < cols; ++j) {
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
  /*
  #pragma parallel for
  for (int i = 2; i < rows; ++i) {
   #pragma omp for
    for (int j = cols; j > 0; --j) {
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
    }*/

 }
