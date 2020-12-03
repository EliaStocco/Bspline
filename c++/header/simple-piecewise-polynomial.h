
#pragma once

#include "polynomial.h"
#include "eliastocco_namespace_piecewise-polynomial.h"
#include "square_range.h"

namespace eliastocco {

    namespace polynomial {

        template < class type , t_dim dimension = 1  , t_dim codimension = 1 >
                class simple_piecewise_polynomial :
                    public polynomial<type,dimension,codimension> {

                    #define thisclass simple_piecewise_polynomial<type,dimension,codimension>
                    #define motherclass polynomial<type,dimension,codimension>
                    
                    using motherclass :: polynomial;
                    
                   public : 
                        typedef typename motherclass :: Type_out Type_out;
                        typedef typename motherclass :: Type_in Type_in;
                
//********************************************************************************************************************

                    //constructor
                    public: 
                            thisclass ( ) 
                                : domain ( [](auto x){return true;} ),
                                    motherclass() {};
                                    
                    //constructor
                    public: 
                            thisclass ( const motherclass pol ) 
                                : domain ( [](auto x){return true;} ),
                                    motherclass(pol) {};
                    
                    //constructor
                    public :
                        thisclass ( basis<type,dimension> * user_basis_function )
                            : thisclass ( motherclass ( user_basis_function ) ){};
                                            
                                            
                    //constructor
                    public: 
                            thisclass ( domain_t<type,dimension> user_domain , const motherclass pol ) 
                                : domain ( user_domain ),
                                    motherclass(pol) {};
                    
//********************************************************************************************************************

                    //valutazione del polinomio in x
                    public :
                         Type_out operator () ( const Type_in x ) const {
                         
                             if ( domain(x) ){   
                                 auto m = static_cast<motherclass>(*this);
                                 return m.operator()(x);
                             } else {
                                 return Type_out();
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
                    
                    //piccola ottimizzazione nel caso si moltiplichi per un polinomio non definito a tratti
                    //in realtà la presenza di questo operatore è necessaria per non creare confusione
                    //tra l'operator* definito qui sotto e quello definito in polynomial.h
                    public : 
                        friend thisclass operator * 
                            (   const thisclass left ,
                                const motherclass right ) {                                
                                thisclass output = static_cast<motherclass>(left) * right;                            
                                output.set_domain( left.get_domain() );                                
                                return output;                                
                            };
                            
                     //piccola comodià: l'operatore è commutativo     
                     public : 
                        friend thisclass operator * 
                            (   const motherclass left ,
                                const thisclass right ) {   
                                return right*left;                                
                            };
                            
                     public :
                        friend thisclass operator * 
                            (   const thisclass left ,
                                const thisclass right ) {
                                
                                thisclass output = static_cast<motherclass>(left) * 
                                                   static_cast<motherclass>(right); 

                                auto intersection_range = [=]( const t_var<type,dimension> x ){         
                                    const auto left_domain  = left.get_domain();
                                    const auto right_domain = right.get_domain();                                    
                                    return left_domain(x)*right_domain(x);   
                                };
                                
                                output.set_domain( intersection_range );
                                
                                return output;
                                
                            };                            
                  
//********************************************************************************************************************

                    //using motherclass :: operator*(const type factor , const thisclass right);
                    // moltiplicazione per uno scalare
                    public :
                        friend thisclass operator * ( const type factor , const thisclass right ){                        
                            
                            motherclass mother = factor * static_cast<motherclass>(right);                            
                            thisclass output (right.get_domain(),mother);                                
                            return output;                            
                            
                        };
                        
//********************************************************************************************************************
                    
                    public :
                        friend thisclass operator + 
                            (   const thisclass left ,
                                const thisclass right ) = delete ;
                                
                    // somma tra due polinomi: non è definibile!
                        
//********************************************************************************************************************
                    public :
                        void set_domain ( domain_t<type,dimension> user_domain ) {   
                            domain=user_domain ; 
                        };
                        
                    public :
                        void set_domain ( const square_range<type,dimension> user_domain ) {   
                            domain=user_domain.lambda ; 
                        };
                        
                    public :
                        inline auto get_domain ( ) const { return domain ; };

//********************************************************************************************************************

                    private :
                        domain_t<type,dimension> domain;
                        
                    #undef thisclass
                    #undef motherclass
                    
                    };
                
        
    }
}
