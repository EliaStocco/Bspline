#pragma once

#include "eliastocco_namespace_polynomial.h"
#include "basis.h"
#include "monomial.h"

namespace eliastocco {

    namespace polynomial {
    
        //template variable
        //devo usare un puntatore perché è una classe astratta!!
        //non è più una classe astratta
        template < class type , t_dim dimension = 1 >
            basis < type , dimension > basis_function_default = monomial_basis<type,dimension>(t_deg<dimension>());

        template < class type , t_dim dimension = 1  >
            basis < type , dimension >* basis_function_default_ptr = &basis_function_default<type,dimension>;
    
    
    }
    
}
