#pragma once
#include <array>
#include <tuple>


#include "eliastocco_namespace_polynomial_tensor_product.h"

namespace eliastocco {
    
    namespace polynomial {

        namespace tensor_product {  
        
        
          template < class type >
             class split_array <type,1,1> {
             
                 public :
                     split_array<type,1,1> (const t_var <type,1+1> C){
                     
                         t_var <type,1> A = C[0];
                         t_var <type,1> B = C[1];

                    
                         
                        value = std :: make_tuple(A,B);
                     
                     }
                     
                public :
                    std :: tuple <
                    t_var <type,1> ,
                    t_var <type,1> > value;
             
             
             };
        
        
         template < class type , t_dim A_dim >
             class split_array<type,A_dim,1> {
             
                 public :
                     split_array<type,A_dim,1> (const t_var <type,A_dim+1> C){
                     
                         t_var <type,A_dim> A ;
                         t_var <type,1> B ;
                         
                         B=C[0];
                         
                         for (auto i = 1 ; i < A_dim; ++i){ A[i-1]=C[i]; }

                         
                        value = std :: make_tuple(A,B);
                     
                     }
                     
                public :
                    std :: tuple <
                    t_var <type,A_dim> ,
                    t_var <type,1> > value;
             
             
             };
        
        template < class type , t_dim B_dim >
             class split_array<type,1,B_dim> {
             
                 public :
                     split_array<type,1,B_dim> (const t_var <type,1+B_dim> C){
                     
                         t_var <type,1> A ;
                         t_var <type,B_dim> B ;
                         
                         A=C[0];
                        
                         for (auto i = 1 ; i < B_dim; ++i){ B[i-1]=C[i]; }
                         
                        value = std :: make_tuple(A,B);
                     
                     }
                     
                public :
                    std :: tuple <
                    t_var <type,1> ,
                    t_var <type,B_dim> > value;
             
             
             };
        
         template < class type , t_dim A_dim , t_dim B_dim >
             class split_array {
             
                 public :
                     split_array<type,A_dim,B_dim> (const t_var <type,A_dim+B_dim> C){
                     
                         t_var <type,A_dim> A ;
                         t_var <type,B_dim> B ;

                         for (auto i = 0 ; i < A_dim; ++i){ A[i]=C[i]; }
                         for (auto i = 0 ; i < B_dim; ++i){ B[i]=C[i+A_dim]; }
                         
                        value = std :: make_tuple(A,B);
                     
                     }
                     
                public :
                    std :: tuple <
                    t_var <type,A_dim> ,
                    t_var <type,B_dim> > value;
             
             
             };
        
        }
        
    }
        
}
