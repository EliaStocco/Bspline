#pragma once

#include "eliastocco_namespace_polynomial.h"
//necessario includerlo
#include "evaluate.h"

namespace eliastocco {

    namespace polynomial {    

        template < class type , t_dim dimension > //, t_dim codimension = 1 >
            class basis : public t_ptr<type,dimension> {

                //default constructor
                public :
                    basis() : t_ptr<type,dimension> (){};
                    
                //copy constuctor            
                public :
                    basis ( const t_ptr<type,dimension> user_basis) 
                        : basis(user_basis , 
                        eliastocco::polynomial::default_polynomial_evaluation<type,dimension,1>){};
                        
                //copy constuctor            
                public :
                    basis ( const t_ptr<type,dimension> user_basis , 
                                t_eval<type,dimension,1> user_evaluate ) 
                        : t_ptr<type,dimension> (user_basis),
                            evaluate_lambda(user_evaluate) {};

                
                            
//*******************************************************************************************************************

                public:
                    virtual type evaluate 
                        ( const t_coeff_array_basis <type,dimension> _coeffs ,
                            const t_var<type,dimension> x ) {

                        return evaluate_lambda(_coeffs,(*this),x);

                    };
                    
                public :
                    virtual const t_func < type , dimension > operator [](const t_deg<dimension> index) const {                    
                        const auto & a = dereference();
                        return a.at(index);
                    };
                    
//*******************************************************************************************************************

                private :
                    inline const auto & dereference () const 
                        { return *( t_ptr<type,dimension> :: get () ) ; }

                

                public :
                    auto size () const {                 
                        const auto & a = dereference();
                        return a.size();                
                    };

                protected :
                    t_eval<type,dimension,1> evaluate_lambda;                  
              


            };
    
    }
        
}
