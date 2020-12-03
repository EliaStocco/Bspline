#pragma once
#include "eliastocco_namespace_Bspline.h"
#include "Bspline.h"

namespace eliastocco {

    namespace polynomial {    
    
        namespace Bspline {   
        
            template < class type , t_dim dimension , t_dim codimension >
                t_eval<type,dimension,codimension> make_deBoor
                    (const unsigned int N,
                         const unsigned int P,
                         const std :: vector<type> knots){
                         
                         //std :: cout << "make_de_Boor" << std :: endl ;
                         
                            //auto P = _basis_function.get_P();
                            //auto N = _basis_function.get_N();
                            //auto t = knots;
                            
                            std :: valarray<type> t (knots.size());
                            std :: copy(knots.begin(),knots.end(),std::begin(t));
                      
                      return [=]( const t_coeff_array<type,dimension,codimension> & coeffs ,
                        const t_ptr<type,dimension> _basis_function , //non usare la reference qui !!!
                            const t_var < type,dimension > _x )-> t_var<type,codimension> {
                            
                    
                     
                            
                    //******************************************
                    
                    std :: cout << "conversion of _coeffs" << std :: endl ;
                    
                        typedef typename std :: conditional<  std :: less<t_degree>{}(1,codimension) ,
                                matrix < dimension , std :: valarray<type> > , 
                                matrix < dimension , type > > :: type _coeff_type;
                    
                        _coeff_type _coeffs = _coeff_type();
                        
                            if constexpr ( codimension > 1 ){
                            
                                 auto from_t_var_2_valarray = [] ( const t_var<type,codimension> in ){

                                     std::valarray<type> out (codimension);
                                     std :: copy(in.begin(),in.end(),std::begin(out));
                                     return out;

                                 };
        
                                 //matrix < dimension , std :: valarray<type> > _coeffs ;
                                 for(auto it = coeffs.begin(); it!= coeffs.end();++it){
                                     _coeffs.insert(std :: make_pair(
                                     (*it).first,from_t_var_2_valarray((*it).second)));
                                 }
                                                
                            } else {
                            
                                _coeffs = coeffs;
                                
                            }
                    
                            
                    //******************************************        
                    
                    /*std :: cout << "conversion of _coeffs" << std :: endl ;
                    
                            typedef typename std :: conditional<  std :: less<t_degree>{}(1,codimension) ,
                                std :: valarray <type> , type > :: type coeff_type;
                            
                            //attenzione ad inizializzare !!
                            coeff_type output = coeff_type ();
                            coeff_type coefficient = coeff_type ();                            
                            
                            if constexpr ( codimension > 1 ){                            
                                output.resize(codimension);
                                coefficient.resize(codimension);   
                            } 
                            */
                    //******************************************
                    
                    std :: cout << "conversion of x and 1.0" << std :: endl ;
                    
                            typedef typename std :: conditional<  std :: less<t_degree>{}(1,dimension) ,
                                std :: valarray <type> , type > :: type input_type;

                            //modifico tipo della variabile posizione
                            input_type x = input_type();
                            input_type uno = input_type();
                            if constexpr ( dimension > 1 ){                                    
                                std :: copy( std :: begin(_x), std :: end(_x), std :: begin(x) );  
                                uno.fill(1.0);
                            } else {                                    
                                x = _x ;  
                                uno = 1.0;
                            }
                            
                     //******************************************
                          //trovo K : indice dell'intervallo a cui appartiene x
                          
                          
                            int K = -1;
                          
                          for(auto i = 0 ; i < t.size()-1;++i){
                              if (t[i]<=x and x<t[i+1]){
                                  K=i;
                                  break;
                              }
                          }
                          
                          if ( K == -1 ){ 
                              std :: cout << "domain error" << std :: endl;
                          }
                        
                        std :: cout << "K = " << K << std :: endl;
                        
                    //******************************************
                            
                    std :: cout << "conversion of d" << std :: endl ;
                                
                            typedef typename std :: conditional<  std :: less<t_degree>{}(1,codimension) ,
                            std :: valarray <std::valarray<type>> , std :: valarray <type> > 
                            :: type knot_type;
                                
                            knot_type d (P+1);
                            
                            if constexpr (codimension > 1 ){                            
                                for(auto i=0;i<P+1;++i){
                                    d[i]=std::valarray<type>(codimension);
                                }                            
                            }
                            
                            // qui inizia l'algoritmo vero e proprio
                                
                            try {
                            
                            std :: cout << "try" << std :: endl ;
                            
                            std :: cout << "P=" << P << std :: endl ;
                            std :: cout << "size=" << _coeffs.size() << std :: endl ;
                               
                                for( auto j=0; j <= P; ++j ){     
                                    std :: cout << "j : " << j ;
                                    d[j];
                                    std :: cout << " | ok" ;
                                    _coeffs.at(j+K-P);  
                                    std :: cout << " | ok" << std :: endl ;
                                    d[j] = _coeffs.at(j+K-P);                                
                                }
                                
                            std :: cout << "d" << std :: endl ;
                                
                                for( auto r=1; r <=P; ++r ){
                                    
                                    std :: cout << "r : " << r << std :: endl ;
                                    for( auto j=P; j>r; --r){
                                    
                                        input_type alpha = ( x - t[j+K-P] ) / ( t[j+1+K-r] - t[j+K-P] );
                                        d[j] = ( uno - alpha ) * d[j-1] + alpha * d[j];
                                    
                                    }                                
                                }
                            std :: cout << "finish" << std :: endl ;
                                
                                //return d.at(P);

                            } catch (  std::out_of_range e ) {
                                std :: cout << "error : deBoor" << std :: endl;
                                std :: cout << e.what() << std :: endl ;      
                                //throw std::out_of_range();

                            }
                            
                            if constexpr (codimension > 1 ){                                   
                                t_var<type,codimension> Output;                                                
                                for(auto i=0;i<codimension;++i){
                                    Output[i]=d[P][i];
                                }       
                                return Output;
                            } else {                            
                                return d[P];    
                            }
                            
                            };
                         
                    };
        
        
        
        }
        
    }
    
}
