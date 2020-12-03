#pragma once
#include "eliastocco_namespace_polynomial.h"

namespace eliastocco {

    namespace polynomial {
    
        template < class type , t_dim dimension >
            using domain_t = std :: function < bool ( const t_var<type,dimension> ) >;
            
        template < class type , t_dim dimension , t_dim codimension >
                class piecewise_polynomial;
                
        template < class type , t_dim dimension > class square_range ;
    
    }
    
}
