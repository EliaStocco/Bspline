#pragma once

#include "dynamic_vector.h"
#include "tensor_product.h"
#include "eliastocco_namespace_polynomial.h"
#include "tensor_product.h"
#include "basis.h"
#include "polynomial.h"


namespace eliastocco {

    namespace polynomial {
    
        template < class type , t_dim dimension = 1 >
            basis < type , dimension > bernstein_basis (const t_deg<dimension> N );
            
        template < class type , t_dim dimension = 1 >
            class make_bernstein ;
            
        template < class type , t_dim dimension = 1 >
            polynomial<type,dimension> associated_bernstein_polynomial 
                ( const t_func<type,dimension> function , basis<type,dimension> & base );
    
            
    }
            
}
        
 namespace eliastocco {

    namespace polynomial {
        
    
        template < class type , t_dim dimension >
            basis < type , dimension > bernstein_basis (const t_deg<dimension> N ){
                
                //costruttore
                return basis < type , dimension > ( make_bernstein<type,dimension>(N).value );

            };
        
        
        template < class type >
            class make_bernstein<type,1> {
                
                public :
                    make_bernstein<type,1> ( const t_deg<1> N ){      
                    
                    
                        //variabile di output
                        value = std :: make_shared < t_basis_function<type,1> > ();
        
                        //calcolo i coefficienti binomiali
                        //https://cp-algorithms.com/combinatorics/binomial-coefficients.html
                        container :: dynamic :: dynamic_vector < 2 , t_deg<1> > C ( { N+1 , N+1 } );

                        //std :: cout << "coefficienti binomiali" << std :: endl ; 
                        C[0][0] = 1;
                        for ( t_deg<1> n = 1; n <= N; ++n) {
                            C[n][0] = C[n][n] = 1;
                            for ( t_deg<1> k = 1; k < n; ++k){
                                C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
                            }
                        }       

                        //std :: cout << "calcolo" << std :: endl ; 
                         //
                        for (auto k=0;k<=N;++k){
                        
                            auto f = [=](const t_var<type,1> x){

                                return C [N][k] * pow ( 1 - x , N - k ) * pow ( x , k );

                            }; 
                        
                            (*value).insert(std :: make_pair(k,f));
                        }
                        
                    }
                
                public :
                    t_ptr < type , 1 > value;
                
            };
        
    
        template < class type , t_dim dimension >
            class make_bernstein {
                
                public :
                    make_bernstein<type,dimension> ( const t_deg<dimension> N ){
                        
                        t_var<t_degree,1> N_0;
                        t_var<t_degree,dimension-1> N_1;
                                
                        tensor_product::split_array<t_degree,1,dimension-1> X (N);
                        std :: tie (N_0,N_1) =  X.value;
            
                        auto basis_1  = *(make_bernstein<type,1> (N_0).value);
                        
                        auto basis_N_inf = *(make_bernstein <type,dimension-1>(N_1).value);
           
                        auto T = tensor_product :: tensor_product<type,1,dimension-1>( basis_1 , basis_N_inf ).value;
                        
                        value = std :: make_shared<decltype(T)>(T);
                        
                    }
                
                public :
                    t_ptr < type , dimension > value;
                        
                        
            };  
      
    
    }
    
}

namespace eliastocco {
    
    namespace polynomial {
        
        template < class type , t_dim dimension >
            polynomial<type,dimension> associated_bernstein_polynomial 
                ( const t_func<type,dimension> function , basis<type,dimension> & base ){
                                
                static_assert ( dimension == 1 , "associated_bernstein_polynomial not implemented for dimension > 2");
                
                auto N = base.size();
                
                
                t_list<type,dimension> coeff (N);
        
                unsigned int k=0;
                for ( auto it = coeff.begin(); it != coeff.end(); ++it){ 
                    //attenzione, N è il numero di polinomi, N-1 è il grado polinomiale
                    (*it) = std :: make_pair ( k , function(static_cast<type>(k)/static_cast<type>(N-1)));
                    k++;
                }
        
                return polynomial<type,dimension>(coeff,&base);
                
                
        }
    }
}
