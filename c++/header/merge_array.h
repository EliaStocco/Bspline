#pragma once
#include <array>


#include "eliastocco_namespace_polynomial_tensor_product.h"

namespace eliastocco {
    
    namespace polynomial {

        namespace tensor_product {   
        
        
           template < class type >
                class merge_array<type,1,1> {
                
                    public :
                        merge_array<type,1,1>(const t_var<type,1> A,
                                                        const t_var<type,1> B){
                          
                            value={A,B};
                                                        
                         };
                        public :
                         t_var<type,1+1> value;
                         
                    };
        
      
                    
            template < class type , t_dim A_dim >
                class merge_array<type,A_dim,1> {
                
                    public :
                        merge_array<type,A_dim,1>(const t_var<type,A_dim> A,
                                                        const t_var<type,1> B){
                                                        
                            t_var< type , A_dim + 1 > C ;
                         for (auto i = 0 ; i < A_dim; ++i){ C[i]=A[i]; }
                         /*for (auto i = 0 ; i < B_dim; ++i){ */C[A_dim]=B;// }
                         
                         value= C;
                                                        
                         };
                         
                        public :
                            t_var<type,A_dim+1> value;
                         
                    };
                    
            template < class type , t_dim B_dim >
                class merge_array<type,1,B_dim> {
                
                    public :
                        merge_array<type,1,B_dim>(const t_var<type,1> A,
                                                        const t_var<type,B_dim> B){
                                                        
                            t_var< type , B_dim + 1 > C ;
                            C[0]=A;
                             for (auto i = 0 ; i < B_dim; ++i){ C[i+1]=B[i]; }
                          
                             value= C;
                                                        
                         };
                         
                        public :
                            t_var<type,1+B_dim> value;
                         
                    };
                    
                    
                          template < class type , t_dim A_dim , t_dim B_dim >
                class merge_array {
                
                    public :
                        merge_array<type,A_dim,B_dim>(const t_var<type,A_dim> A,
                                                        const t_var<type,B_dim> B){
                                                        
                             t_var< type , A_dim + B_dim > C ;
                             for (auto i = 0 ; i < A_dim; ++i){ C[i]=A[i]; }
                             for (auto i = 0 ; i < B_dim; ++i){ C[i+A_dim]=B[i]; }
                         
                              value = C;
                                                        
                         };
                         
                         public :
                            t_var< type , A_dim + B_dim > value;                         
                         
                    };
                    
         
                         
                          
           
        }
        
    }
    
}
