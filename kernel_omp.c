#include <omp.h>
#include<math.h>
#include<stdio.h>
#include "kernel.h"

void kernel_omp(int *input, int *ref, int64_t rows, int64_t cols, int penalty) {
  const int TILESIZE = 128;
  int i;
  for(i = 1; i <= rows-TILESIZE; i += TILESIZE){
#pragma omp parallel for 
    for(int j = 1; j <= i; j += TILESIZE){
      int m = i;      

      int temp = j;
      if(m > 0){
      for (int k = m; k < m + TILESIZE; ++k){                                                                                                                                        
	for(int l = temp; l < temp + TILESIZE; ++l){                                                                                                                    
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
         m -= TILESIZE;
      }
    }
  }

#pragma omp barrier
   int x,y,z;
    for(int j = 1 + TILESIZE; j <= cols - TILESIZE; j+= TILESIZE){
     #pragma omp parallel for 
      for(x = rows-TILESIZE; x <= j; x-= TILESIZE){
	int n = j;
       	for (y = x; y < x + TILESIZE; ++y){                       
           for(z = n; z < n + TILESIZE; ++z){
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
	
	n += TILESIZE;
      }
    }
 }
