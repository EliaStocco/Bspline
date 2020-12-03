#pragma once

#include <iostream>

#include "eliastocco_namespace_polynomial.h"
#include <valarray>
#include <iterator>

namespace eliastocco {

    namespace polynomial {   
            
            //template function
            template < class type , t_dim dimension , t_dim codimension >
                /*std :: valarray<type>*/t_var<type,codimension>   default_polynomial_evaluation
                    ( const t_coeff_array<type,dimension,codimension> & _coeffs ,
                        const t_ptr < type , dimension > _basis_function , //non usare la reference qui !!!
                            const t_var < type,dimension > x ){

                            //std :: cout << "default_polynomial_evaluation" << std :: endl;
       
                            typedef typename std :: conditional<  std :: less<t_degree>{}(1,codimension) ,
                                std :: valarray <type> , type > :: type coeff_type;

                            //attenzione ad inizializzare !!
                            coeff_type output = coeff_type ();
                            coeff_type coefficient = coeff_type (); 
                            
                            if constexpr ( codimension > 1 ){
                            
                                output.resize(codimension);
                                coefficient.resize(codimension);                                    
                                       
                            } 


                            try {

                               // basis_function deve essere già dotata di tutte le std::functional necessarie
                                for ( auto it = _coeffs.begin(); it != _coeffs.end(); ++it ) {                    


                                    if constexpr ( codimension > 1 ){
                                    
                                        auto & c = it->second;
                                        std :: copy( std :: begin(c), std :: end(c), std :: begin(coefficient) );
                                    
                                    } else {
                                    
                                        coefficient = it->second ;
                                    
                                    }
                                    
                                    output+= /*_coeffs.at*/coefficient * (*_basis_function)[it->first](x);
                                }

                                //return output;

                            } catch (  std::out_of_range e ) {
                                std :: cout << "error : default_polynomial_evaluation" << std :: endl;
                                std :: cout << e.what() << std :: endl ;      
                                //throw std::out_of_range();

                            }
                            
                            if constexpr ( codimension > 1 ){
                                t_var<type,codimension> Output;                            
                                std :: copy( std :: begin(output), std :: end(output), std :: begin(Output) );
                                return Output;
                            } else {                                
                                return output;                            
                            }                            

            };
       
    }
   
}
