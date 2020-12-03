#pragma once
//#include "eliastocco_namespace_piecewise-polynomial.h"
#include "dynamic_vector.h"

namespace eliastocco {

    namespace polynomial {    
    
        namespace Bspline {
        
            //typedef /*long*/ unsigned int ut_dim;
            //typedef /*long unsigned*/ int t_dim;
            typedef eliastocco :: container :: dynamic :: t_dim  t_dim;
            typedef eliastocco :: container :: dynamic :: ut_dim ut_dim;
                    
            template < class type > class knot_vector;
            
            template < class type , ut_dim dimension , ut_dim codimension > class Bspline;
            
            
        }
        
    }
    
}
