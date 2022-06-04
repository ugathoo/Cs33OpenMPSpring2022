#include <omp.h>
#include<math.h>
#include<stdio.h>
#include "kernel.h"

void kernel_omp(int *input, int *ref, int64_t rows, int64_t cols, int penalty) {
  const int TILESIZE = 4;
  const int rminT = rows - TILESIZE;
  //const int rplusT = rows + TILESIZE;
  int i;
  for(i = 1; i < rminT; i += TILESIZE){
    #pragma omp parallel for 
    for(int j = 1; j < i + TILESIZE; j += TILESIZE){
      int row = i + (1 - j);
	for (int k = row; k < row + TILESIZE; ++k){                                                      
	  for(int l = j; l < j + TILESIZE; ++l){                                                                //printf("i: %d  row: %d  j: %d  k: %d l: %d\n", i,row,j,k,l);                                      
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
	//m -= TILESIZE;
    }
  }

  #pragma omp barrier
  int j;
  for(j = 1; j < cols; j+= TILESIZE){
     #pragma omp parallel for 
    for(int x = rminT; x > j - TILESIZE; x-= TILESIZE){
      int col = j + rminT - x;
      for (int y = x; y < x + TILESIZE; ++y){                       
	for(int z = col; z < col + TILESIZE; ++z){
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
      // n += TILESIZE;
    }
  }
}
/*
void kernel_omp(int *input, int *ref, int64_t rows, int64_t cols, int penalty) {
  const int TILESIZE = 16;
  int i;
  for(i = 1; i < rows; i += TILESIZE){
    //#pragma omp parallel for 
    int m = i;
    for(int j = 1; j <= i; j += TILESIZE){
      int temp = j;
      if(m > 0){
	for (int k = m; k < m + TILESIZE; ++k){                                                                      for(int l = temp; l < temp + TILESIZE; ++l){                                                  
	    printf("(%d, %d) k: %d  l: %d\n", m, j, k, l);
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
	printf ("(%d,%d)\n", m, j);
	m -= TILESIZE;
      }
    }
  }

  //#pragma omp barrier
  int j;
  for(j = 1 + TILESIZE; j < cols; j+= TILESIZE){
    int n = j;
    // #pragma omp parallel for 
    for(int x = rows-TILESIZE; x <= j; x-= TILESIZE){
      for (int y = x; y < x + TILESIZE; ++y){                       
	for(int z = n; z < n + TILESIZE; ++z){
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
*/
