
#pragma once

#include <valarray>
#include "simple-piecewise-polynomial.h"

namespace eliastocco {

    namespace polynomial {        
    
        template < class type , t_dim dimension = 1  , t_dim codimension = 1 >
                class piecewise_polynomial :
                    public std :: valarray <
                        simple_piecewise_polynomial<type,dimension,codimension> > {       
                        
                        #define thisclass piecewise_polynomial<type,dimension,codimension>
                        #define motherclass std::valarray<simple_piecewise_polynomial<type,dimension,codimension>>
                        
                        //uso di default i costruttori di std :: valarray
                        using motherclass :: valarray;
                        
                        public :                             
                            typedef simple_piecewise_polynomial<type,dimension,codimension> pol_t;
                            typedef typename pol_t :: Type_out Type_out;
                            typedef typename pol_t :: Type_in Type_in;
                                
                        //constructor
                        public: 
                            thisclass ( const motherclass val ) 
                                : motherclass( val ) {};
                                    
                       //constructor
                        public: 
                            thisclass ( const pol_t pol ) 
                                : motherclass( {pol} ) {};
                                    
//********************************************************************************************************************

                        //valutazione del polinomio in x
                        public :
                             Type_out operator () ( const Type_in x ) const {
                             
                                 //std :: cout << "piecewise-polynomial operator()" << std :: endl;
                             
                                typedef typename std :: conditional<  std :: less<t_degree>{}(1,codimension) ,
                                    std :: valarray <type> , type > :: type coeff_type;
                                
                                coeff_type output = coeff_type ();                                
                                const auto This = static_cast<motherclass>(*this);                                
                                //std :: cout << "size:" << motherclass :: size() << std :: endl;
                                 
                                for(auto it = std :: begin(This); it != std::end(This); ++it ){                                    
                                    output += (*it)(x);
                                }
                                
                                if constexpr ( codimension > 1 ){
                                    Type_out Output;                            
                                    std :: copy( std :: begin(output), std :: end(output), std :: begin(Output) );
                                    return Output;
                                } else {                                
                                    return output;                            
                                }                               
              
                            };

                        //vectorized polinomial evaluation
                        public :
                        template < template < class ... > class vect = std :: vector >
                            vect <Type_out> operator () ( const vect < Type_in > x ) const { 
                                vect <Type_out> y (x.size());
                                std :: transform( std :: begin(x), std :: end(x), std :: begin(y),
                                    [&]( auto a ) {return (*this)(a);});   
                                return y;

                            }; 

//********************************************************************************************************************
                                        
                    //definisco la somma tra polinomi definiti a tratti
                    public :
                        friend thisclass operator + ( const thisclass left , const thisclass right ){
                            thisclass output ( left.size() + right.size() );
                            auto it = std :: copy( std :: begin(left), std :: end (left) , std :: begin(output));
                            //attenzione al std :: prev !!
                            std :: copy ( std :: begin(right), std :: end (right) , std :: prev(it) );
                            return output;                            
                        };
                        
                    //definisco il prodotto di polinomi definiti a tratti
                    public :
                        friend thisclass operator * ( const pol_t left , const thisclass right ){   
                        
                            //constructor specifying dimension
                            thisclass output ( right.size() );                        
                            //contatore
                            unsigned int k = 0 ;
                            // ciclo sui polinomi di right
                            for ( auto r = std :: begin(right) ; r != std::end(right) ; ++r ){   
                                //moltiplicazione tra polinomi definiti a tratti
                                output[k++]= left*(*r);
                            }                  
                            return output;  
                        };
                        
                    //definisco il prodotto di polinomi definiti a tratti
                    public :
                        friend thisclass operator * ( const thisclass left , const thisclass right ){   
                        
                            //constructor specifying dimension
                            thisclass output ( left.size() * right.size() );                        
                            //contatore
                            unsigned int k = 0 ;                            
                             // ciclo sui polinomi di left
                            for ( auto l = std :: begin(left) ; l != std::end(left) ; ++l ){
                                // ciclo sui polinomi di right
                                for ( auto r = std :: begin(right) ; r != std::end(right) ; ++r ){   
                                    //moltiplicazione tra polinomi definiti a tratti
                                    output[k++]= (*l)*(*r);
                                }                        
                            }                            
                            return output;  
                        };
                        
                        public :
                            friend thisclass operator * ( const type factor , const thisclass right ){
                                                           
                                motherclass output (right.size());
                                
                                for(auto i = 0; i < right.size(); ++i ){ 
                                    output[i]=factor*static_cast<motherclass>(right)[i];
                                }
                                
                                return output;
                            
                            };   

                        #undef thisclass
                        #undef motherclass
                    
                    };               
        
    }
}
