#include <omp.h>
#include<math.h>
#include<stdio.h>
#include "kernel.h"

void kernel_omp(int *input, int *ref, int64_t rows, int64_t cols, int penalty) {
  const int TILESIZE = 128;
  int j,k,l; 
  for(int i = 1; i < rows; i += TILESIZE){
    int m = i;
    for(j = 1; m > 0; j += TILESIZE){
      for (k = m; k < m + TILESIZE; ++k){                                                                                                                                        
	for(l = j; l < j + TILESIZE; ++l){                                                                                                                                        
	  int64_t idx = k * cols + l;                                                                                                                                             	  int sub = idx - cols;                                                                                                                                                   
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
      printf("(%d,%d)\n", m,j);
      m -= TILESIZE;
    }
  }

    int x,y,z;
    for(int j = 1 + TILESIZE; j <= cols; j+= TILESIZE){
      int n = j;
      for(x = rows; x <= n; x-= TILESIZE){
       	for (y = x; y < x + TILESIZE; ++y){                                                                                                                                                                                                    for(z = n; z < n + TILESIZE; ++z){                                                                                                                                                                                          
            int64_t idx = y * cols + z;                                                                                                                                                                                      
            int sub = idx - cols;                                                                                                                                                                                                                 int64_t idxNW = sub - 1;                                                                                                                                                                                                  
            int64_t idxN = sub;                                                                                                                                                                                                       
            int64_t idxW = idx - 1;                                                                                                                                                                                                   
            int r = ref[idx];                                                                                                                                                                                                         
            int inputNW = input[idxNW];                                                                                                                                                                                               
            int inputW = input[idxW];                                                                                                                                                                                                
            int inputN = input[idxN];                                                                                                                                                                                                 
            input[idx] = maximum(inputNW + r, inputW - penalty, inputN - penalty);                                                                                                                                                    
          }                                                                                                                                                                                                                           
	  } 
	printf("(%d,%d)\n", x,n);
	n += TILESIZE;
      }
    }



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
  /* 
 int TILESIZE = 16;
 int j, k, l;
#pragma omp parallel
 for(int i = 1; i < rows-TILESIZE; i += TILESIZE){
   #pragma omp for private (k,l,j)
    for(j = i; j >= 1; j-= TILESIZE){
      printf("Thread#: %d      i: %d     j: %d\n ",omp_get_thread_num(), i, j);
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
    }*/
   /*
#pragma omp barrier
 int x,y,z;
#pragma omp parallel
 for(int f = cols - TILESIZE; f > 0; f -= TILESIZE){
  #pragma omp for private (y,x,z)
   for(x = TILESIZE; x <= rows-TILESIZE; x+= TILESIZE){
     printf("Thread#: %d       x: %d        f: %d\n",omp_get_thread_num(),x,f);
     for (y = x; y < x + TILESIZE; y++){
       for(z = f; z < f + TILESIZE; z++){
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
 }*/

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
