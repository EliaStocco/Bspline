#pragma once

#include <iostream>
#include <functional>
#include <tuple>

#include "split_array.h"
#include "merge_array.h"

#include "eliastocco_namespace_polynomial.h"
#include "eliastocco_namespace_polynomial_tensor_product.h"

namespace eliastocco {
    
    namespace polynomial {

        namespace tensor_product {   
           
             template < class type , t_dim A_dim, t_dim B_dim >
                 class tensor_product {
                 
                     public :
                         tensor_product<type,A_dim,B_dim> ( const t_basis_function<type,A_dim> A ,
                                        const t_basis_function<type,B_dim> B ) : value(){
                                        
                           
                            typedef typename decltype(A)::mapped_type lambda_A;
                            typedef typename decltype(B)::mapped_type lambda_B;                         

                            typedef typename lambda_A :: result_type result_A;
                            typedef typename lambda_B :: result_type result_B;
                            static_assert ( std::is_same<result_A, result_B>::value ,
                                "EliaStocco: different result types" );
                            typedef result_A result_type;                            

                            typedef result_type Type;

                            static_assert ( std::is_same<Type, type>::value ,
                                "EliaStocco: result type different form template type" );

                            typedef typename std :: array < type , A_dim + B_dim > parameter;  

                            typename decltype(value)::key_type key;
                            typename decltype(value)::mapped_type lambda;

                            typedef typename decltype(value)::mapped_type::argument_type argument;

                            static_assert ( std::is_same<argument, parameter>::value ,
                                "EliaStocco: argument type different form parameter type" );
                                

                            for ( auto a = A.begin(); a!=A.end();++a){
                                for ( auto b = B.begin(); b!=B.end();++b){
                                    
                                    merge_array<t_degree,A_dim,B_dim> app ((*a).first , (*b).first);
                                    key =app.value;

                                    
                                    auto a_2 = (*a).second;
                                    auto b_2 = (*b).second;

                                    lambda = [=](const parameter x ){

                                        t_var<type,A_dim> x_0;
                                        t_var<type,B_dim> x_1;

                                        split_array<type,A_dim,B_dim> X (x);
                                        std :: tie (x_0,x_1) =  X.value;


                                        auto y_0 = a_2(x_0);
                                        auto y_1 = b_2(x_1);
                                        return y_0*y_1;

                                    };                                   

                                    value.insert(std :: make_pair (key,lambda));

                                }
                            
                            }                                        
                                        
                        };
                                        
                public :
                    t_basis_function <type,A_dim+B_dim> value;
                    
                    
                 };
                    
                
            };
           
        }
        
    }
