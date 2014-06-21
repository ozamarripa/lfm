//------------------------------------------------------------------------------
//
// kernel:  vadd  
//
#pragma OPENCL EXTENSION cl_amd_printf : enable

__kernel void reset_fuz( __global int* fuz,__global int* del )               
{  
	del[0]=-1;                                        
   int i = get_global_id(0);               
   fuz[i] = 0;                 
}                                          


__kernel void detect_overlap(
	__global int* fuz,
	__global int* pos_size_big,
	__global int* communities,
	__global int* membership 
	)               
{                                          
   int i = get_global_id(0);
    int pos=i*2;  
                  
   	 if(pos_size_big[i*2]!=-1 ){
        
	     for (int j=0; j<pos_size_big[(i*2)+1]; j++){
            //fuz[i]=pos_size_big[(i*2)+1];
			if (membership[ communities[pos_size_big[i*2]+j ] ] >1){ 	
					fuz[i]=fuz[i]+1;//sets the number of nodes of each community that belongs to more than one community
			}    
		}             
	}	
} 

__kernel void detect_overlap_record(
	__global int* fuz,
	__global int* pos_size_big,
	__global int* communities,
	__global int* membership,
	__global int* del 
	)               
{                                          
   int i = get_global_id(0);
    int pos=pos_size_big[i*2];  
                  
   	 if(pos!=-1 ){
         int size=pos_size_big[(i*2)+1];
         //int nodes[size];
         //for (int j=0; j<size; j++){
         //	nodes[j]=communities[pos+j ];
         //}
	     for (int j=0; j<size; j++){
            
			if (membership[ communities[pos+j ] ] >1){ 	
					fuz[i]=fuz[i]+1;//sets the number of nodes of each community that belongs to more than one community
			}
			if ( fuz[i]==size ){
				del[0]=i;
			}    
		}             
	}	
} 

__kernel void erase_community(
	__global int* pos_size_big,
	__global int* communities,
	__global int* membership,
	const int id_com,
	const int pos_com 
	)
{
	pos_size_big[id_com*2]=-1;
	int i = get_global_id(0); //size of community
	membership[ communities[pos_com+i] ]--;
}

