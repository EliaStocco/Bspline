#pragma once

#include "tensor_product.h"
#include "eliastocco_namespace_polynomial.h"
#include "basis.h"
#include "evaluate.h"


namespace eliastocco {

    namespace polynomial {
    
    
        template < class type , t_dim dimension >
            basis < type , dimension > monomial_basis (const t_deg<dimension> N ){
                
                //costruttore
                return basis < type , dimension > ( make_monomials<type,dimension>(N).value );

            };
        
        
        template < class type >
            class make_monomials<type,1> {
                
                public :
                    make_monomials<type,1> ( const t_deg<1> N ){                        
                        
                         value = std :: make_shared < t_basis_function < type , 1 > >();
                     
                         for (auto i=0;i<=N;++i){
                    
                            (*value).insert(std :: make_pair(i,[=](const t_var<type,1> x){ return pow(x,i); }));
                         }
                        
                    }
                
                public :
                    t_ptr < type , 1 > value;
                
            };
        
    
        template < class type , t_dim dimension >
            class make_monomials {
                
                public :
                    make_monomials<type,dimension> ( const t_deg<dimension> N ){
                        
                        t_var<t_degree,1> N_0;
                        t_var<t_degree,dimension-1> N_1;
                                
                        tensor_product::split_array<t_degree,1,dimension-1> X (N);
                        std :: tie (N_0,N_1) =  X.value;
            
                        auto basis_1  = *(make_monomials<type,1> (N_0).value);
                        
                        auto basis_N_inf = *(make_monomials <type,dimension-1>(N_1).value);
           
                        auto T = tensor_product :: tensor_product<type,1,dimension-1>( basis_1 , basis_N_inf ).value;
                        
                        value = std :: make_shared<decltype(T)>(T);
                        
                    }
                
                public :
                    t_ptr < type , dimension > value;
                        
                        
            };  
      
    
    }
    
}
