#pragma once

#include "eliastocco_namespace_Bspline.h"


namespace eliastocco {

    namespace polynomial {

        namespace Bspline {
            
            template < class type , ut_dim dimension > 
                std :: array<type,dimension> operator* 
                    (const std :: array<type,dimension> left , const type right ){
                    
                        std :: array<type,dimension> output;

                        for(auto i=0; i<dimension ; ++i){
                            output[i]=left[i]*right;
                        }
                        return output;
                    
                };
            
            template < class type , ut_dim dimension > 
                std :: array<type,dimension> operator* 
                    (const type left , const std :: array<type,dimension> right ){
                    
                        std :: array<type,dimension> output;

                        for(auto i=0; i<dimension ; ++i){
                            output[i]=left*right[i];
                        }
                        return output;
                    
                };
            
            template < class type , ut_dim dimension > 
                std :: array<type,dimension> operator+ 
                    (const std :: array<type,dimension> left , const std :: array<type,dimension> right ){
                    
                        std :: array<type,dimension> output;

                        for(auto i=0; i<dimension ; ++i){
                            output[i]=left[i]+right[i];
                        }
                        return output;
                    
                };
            
        }
        
    }
    
}
