
#pragma once

#include "eliastocco_namespace_piecewise-polynomial.h"

namespace eliastocco {

    namespace polynomial {

        template < class type >
            class square_range<type,1> {
            
                public :
                    square_range ( const t_var<type,1> left , const t_var<type,1> right ){
                    
                        assert ( left >= right && "square_range error : left < right" );                    
                        
                        lambda = [=] 
                            ( const t_var<type,1> x ){ 
                                     
                                return ( x < left or x >= right ) ? false : true ;
                  
                            };        
                        
                        };                
                    
                public :
                    domain_t<type,1> lambda;
            
            };
            
        template < class type , t_dim dimension >
            class square_range {
            
                public :
                    square_range ( const t_var<type,dimension> left , const t_var<type,dimension> right ){
                    
                         for ( auto i = 0 ; i < dimension ; ++i ){                                         
                                assert ( left[i] >= right[i] && 
                                    "square_range error : left < right" );                    
                          }
                    
                        lambda = [=] 
                            ( const t_var<type,dimension> x ){                                                
                                for ( auto i = 0 ; i < dimension ; ++i ){  
                                     if ( x[i] < left[i] or x[i] >= right[i] ){ return false ;} 
                                }
                                return true ;                                
                            };        
                        
                        };                
                    
                public :
                    domain_t<type,dimension> lambda;
            
            };
            
    }
        
}           
