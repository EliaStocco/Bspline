#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <string>
#include <set>

#include "eliastocco_namespace_polynomial.h"
#include "default_basis.h"
#include "basis.h"
#include "evaluate.h"
//#include "Rosetta-Code-fft.h"

namespace eliastocco {

    namespace polynomial {
    
            template < class type , t_dim dimension , t_dim codimension >
                class polynomial {

                    #define thisclass polynomial<type,dimension,codimension>    
                    
                    public :
                        typedef t_var<type,codimension> Type_out;
                        typedef t_var<type,  dimension> Type_in;

                    //default constructor
                    public :
                        thisclass () noexcept 
                            : thisclass( &basis_function_default < type,dimension > ){};  

                    public :
                        thisclass ( basis<type,dimension> * user_basis_function ) noexcept
                            : thisclass ( t_list<Type_out,dimension> () , 
                                            user_basis_function ){};
                                            
                 /*   public :
                        thisclass ( t_coeff_array<type,dimension,codimension> coe,
                                    basis<type,dimension> * user_basis_function ,
                                    t_eval<type,dimension,codimension> user_evaluate ) noexcept
                            : thisclass ( t_list<Type_out,dimension> () , 
                                            user_basis_function , user_evaluate ){                                            
                                            
                                        _coeffs = coe ;
                                        
                                    };    
                                         */   
                                        
                                            
                     public :
                        thisclass ( basis<type,dimension> * user_basis_function ,
                                    t_eval<type,dimension,codimension> user_evaluate ) noexcept
                            : thisclass ( t_list<Type_out,dimension> () , 
                                            user_basis_function , user_evaluate ){};

                    //list : lista dei coefficienti
                    //user constructor                 
                    public :
                        thisclass (const t_list<Type_out,dimension> list) noexcept 
                            : thisclass( list , 
                                    &eliastocco :: polynomial :: basis_function_default < type,dimension > ){};
                                    
 
                         
                    public :
                        thisclass ( const t_list<Type_out,dimension> list,
                                        basis<type,dimension> * user_basis_function ) noexcept 
                                : thisclass ( list , 
                                    user_basis_function ,
                                        eliastocco :: polynomial :: 
                                            default_polynomial_evaluation<type,dimension,codimension> ) {};
                        
                //user constructor
                    public :
                        thisclass ( const t_list<Type_out,dimension> list,
                                        basis<type,dimension> * user_basis_function ,
                                        const t_eval<type,dimension,codimension> user_evaluate ) noexcept 
                            : _basis_function (user_basis_function),
                                _coeffs(),
                                evaluate_lambda(user_evaluate) {  

                            for ( auto it = list.begin(); it != list.end(); ++it){ _coeffs.insert((*it)); }

                        };

//********************************************************************************************************************

                    public :
                        Type_out get ( const t_deg < dimension > index , const bool check = true ) const { 
                        
                            //Type_out a = Type_out();
                            try { 
                                return _coeffs.at(index);                                 
                            } catch ( std::out_of_range e ){  
                            
                                if ( check ){
                                    std :: cout << "error : polynomial.get" << std :: endl;
                                    std :: cout << e.what()<< std :: endl;  
                                }
                                return Type_out();
                            } 
                            //return a;
                            
                        };
                 
                    public : // piccola modifica
                        bool set (const t_deg < dimension > index, const Type_out value , const bool check = true ) {

                            try { 
                                (*_basis_function)[index]; 
                            } catch ( std::out_of_range e ){  
                                if ( check ){
                                    std :: cout << "error : polynomial.set" << std :: endl;
                                    std :: cout << e.what()<< std :: endl;  
                                    return false;
                                }
                            } 

                            _coeffs[index]=value; 

                            return true;

                        };
                        
//********************************************************************************************************************
                     public :
                        inline void clear () { _coeffs.clear () ; }
                        
//********************************************************************************************************************

                    public :
                        inline auto coeffs () const { return _coeffs; }
                        inline auto basis () const { return *_basis_function; }
                       
                        template < template < class ... > class vect = std :: vector >
                            vect < Type_out > all_coeffs () const {
                        
                                vect < Type_out >  output (_basis_function->size());

                                for ( auto i = 0 ; i < output.size() ; ++ i ){                            
                                    output[i] = get(i,false);                                
                                }

                                return output;
                        
                        };
//********************************************************************************************************************

                     public : 
                        friend thisclass multiply 
                            (   const thisclass left ,
                                const thisclass right )  {
                                                                      
                            //ridefinisco il prodotto tra polinomi
                            // I tried the FFT way but I failed

                            // dovrei controllare che i polinomi siano definiti rispetto alla stessa base
                            if ( left._basis_function != right._basis_function ) {
                                assert ("polynomial :: operator * error; coeff array with different basis" );
                            }
                            auto basis = left._basis_function ;
                            // questo metodo di moltiplicare i polinomi funziona solo per polinomi definiti
                            // rispetto ad una base di monomi
                            // per adesso non aggiungo nessun controllo sul "tipo di base"
                            // perché dovrei implementare ancora tutto

                            // la condizione su dimension può essere indebolita
                            // ma per ora implementiamo la moltiplicazione fra polinomi
                            // solo nel caso 1D 
                            static_assert ( codimension == 1 and dimension == 1 );

                            // creo un polinomio e gli assegno i coefficienti calcolati
                            // solo se questi non sono nulli, 
                            // per ottimizzare spazio e rendere l'output più leggibile
                            thisclass output ( basis ) ;

                            // ciclo sui coefficienti dei due polinomi
                            auto L_c = left.coeffs();
                            auto R_c = right.coeffs();

                            // ciclo sui coefficienti di left
                            for ( auto l = L_c.begin() ; l != L_c.end() ; ++l ){

                                auto l_deg   = (*l).first;
                                auto l_coeff = (*l).second;

                                // ciclo sui coefficienti di right
                                for ( auto r = R_c.begin() ; r != R_c.end() ; ++r ){

                                    auto r_deg   = (*r).first;
                                    auto r_coeff = (*r).second;

                                    auto degree = l_deg + r_deg;

                                    auto prev = output.get(degree,false);

                                    output.set(degree, prev + l_coeff*r_coeff );

                                }                        
                            }

                            return output;
                        
                        };                       
                            
                    //polynomial multiplication
                    public :
                        //typename std::enable_if < codimension == 1 and dimension == 1 , thisclass > :: type 
                        friend thisclass operator * ( const thisclass left , const thisclass right ){                        
                            return multiply(left,right) ;
                        }; 
                        
//********************************************************************************************************************
                    
                    // moltiplicazione per uno scalare
                    public :
                        friend thisclass operator * ( const type factor , const thisclass right ){                        
                            
                            auto output = right;                            
                            auto coeffs = output.coeffs();
                            
                            for ( auto it = coeffs.begin() ; it != coeffs.end() ; ++it ){                                
                                auto deg = (*it).first;
                                auto coe = (*it).second;
                                output.set(deg, coe*factor);                            
                            }
                            
                            return output;
                            
                        };
                            
//********************************************************************************************************************
                                            
                    // somma tra due polinomi
                    public :
                        friend thisclass operator + ( const thisclass left , const thisclass right ){ 
                            //std :: cout << "polynomial :: operator +" << std :: endl;
                            // dovrei controllare che i polinomi siano definiti rispetto alla stessa base
                            if ( left._basis_function != right._basis_function ) {
                                assert ("polynomial :: operator * error; coeff array with different basis" );
                            }
                            auto basis = left._basis_function ;
                        
                            std :: set < t_deg < dimension > > degree;
                            
                            auto L_c = left.coeffs();
                            auto R_c = right.coeffs();
                            
                            //left
                            for  (auto it = L_c.begin(); it != L_c.end(); ++it ){
                                degree.insert((*it).first);
                            }
                            //right
                            for  (auto it = R_c.begin(); it != R_c.end(); ++it ){
                                degree.insert((*it).first);
                            }
                            
                            //degree contiene tutti e solo i gradi del nuovo polinomio                            
                            thisclass output (basis);
                            
                            for (auto it = degree.begin(); it != degree.end();++it){
                                auto deg = (*it);
                                //auto left_c  = left.get(deg,false);
                                //auto right_c = right.get(deg,false);
                                output.set(deg,0.0/*left_c+right_c*/);
                            }
                            
                            return output;
                            
                        };

//********************************************************************************************************************
                    //valutazione del polinomio in x
                    public :
                        Type_out operator () ( const Type_in x ) const {
                        
                            //std :: cout << "polynomial operator()" << std :: endl;
                        
                            return evaluate_lambda(_coeffs,*_basis_function,x);
                            //auto a = _basis_function -> evaluate(_coeffs,x);
                        
                            //return _basis_function -> evaluate(_coeffs,x);
                            //return (*_basis_function).evaluate (_coeffs,x); 
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

                    public :                
                        eliastocco :: polynomial :: basis <type,dimension> * _basis_function;
                        t_coeff_array <type,dimension,codimension> _coeffs;
                        
                    private :
                        t_eval<type,dimension,codimension> evaluate_lambda;

                    #undef thisclass                

                };
                
       
    }
}