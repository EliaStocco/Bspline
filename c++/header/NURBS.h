#pragma once

#include "Bspline.h"

namespace eliastocco {

    namespace polynomial {
        
        namespace NURBS {
            
            template < class type , t_dim dimension = 1 , t_dim codimension = 1 >
                class NURBS 
                        : public eliastocco :: polynomial :: Bspline :: Bspline<type,dimension,codimension+1> {
                    
                        typedef container::dynamic::dynamic_vector < dimension , type > weight_t;
                        
                        #define motherclass eliastocco::polynomial::Bspline::Bspline<type,dimension,codimension+1>
                        
                        using motherclass :: Bspline;
                            
                        private :
                            weight_t weight;
                            
                        #undef motherclass
                    
                };
            
        }
        
    }

}
