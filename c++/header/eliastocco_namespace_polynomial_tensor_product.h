#pragma once

#include <memory>
#include <functional>
#include <vector>
#include <cmath>
#include <array>
#include <map>

#include "eliastocco_namespace_polynomial.h"


namespace eliastocco {

    namespace polynomial {
    
        namespace tensor_product {
            
            template < class type , t_dim A_dim , t_dim B_dim > class merge_array ;
                
            template < class type , t_dim A_dim , t_dim B_dim > class split_array ;
                
            template < class type , t_dim A_dim, t_dim B_dim > class tensor_product ;
    

        }
        
    }
        
}
