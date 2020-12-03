
#pragma once
#include <vector>
#include <stdexcept>
#include "eliastocco_namespace_Bspline.h"

namespace eliastocco {

    namespace polynomial {    
    
        namespace Bspline {
        
            template < class type >
                class knot_vector
                    : public std :: vector <type> {
                
                        #define thisclass knot_vector<type>
                        #define motherclass std :: vector <type>
                        
                        public :
                            thisclass(){};
                            
                        public :
                            thisclass( const unsigned int p,
                                       const unsigned int n,
                                       std :: vector <type> v )
                                       : pol_order(p),
                                           basis_card(n), 
                                            motherclass(v) {
                                            
                                if ( p+n+1 != v.size() ){
                                    std :: cout << "knot_vector" << std :: endl;
                                    std :: cout << "p : " << p << std :: endl ;
                                    std :: cout << "n : " << n << std :: endl ;
                                    std :: cout << "v : " << v.size() << std :: endl ;
                                    std :: cout << "knot_vector error : wrong vector size" << std :: endl ;
                                    throw std::exception ();
                                }
                                       
                                       
                            };
                            
                        public :
                            inline unsigned int p () const { return pol_order ; }
                            
                        public :
                            inline unsigned int n () const { return basis_card ; }
                            
                        public :
                            inline motherclass knots () const { return (*this) ; }
                            
                        private :
                            unsigned int pol_order; //grado polinomiale
                            unsigned int basis_card;//cardinalità della base
                        
                        #undef thisclass
                        #undef motherclass
                };
                
            }
            
        }

}
