#pragma once

#include <memory>
#include <functional>
#include <vector>
#include <cmath>
#include <array>
#include <map>


namespace eliastocco {

    namespace polynomial {
    
        typedef unsigned int t_dim;   //potrei modificarlo con long o short
        typedef unsigned int t_degree;//potrei modificarlo con long o short
        
        
        
        //grado polinomiale
        template < t_dim dimension >
          using t_deg = typename std::conditional < std :: less<t_degree>{}(1,dimension) ,
                                      std :: array < t_degree, dimension > , t_degree > :: type ;        
             
        //variabile indipendente 
        template < class type , t_dim dimension >
            using t_var = typename std::conditional < std :: less<t_degree>{}(1,dimension) ,
                                    std :: array < type, dimension > , type > :: type ;        
         
        //contenitore per i coefficienti e i polinomi
        //potrei modificarlo con un std :: map
        template < t_dim dimension  , class contained >
            using matrix = std :: map < t_deg < dimension > , contained >;
            
              
        //polinomio
        template < class type , t_dim dimension >
            using t_func = std :: function < type ( t_var < type , dimension > ) >; 

       /*template < class type , t_dim dimension >
            using t_list = std :: vector < std :: pair < t_deg < dimension > , type > >;*/
             
        


        //array di polinomi : base polinomiale
        //vettore,matrice, cubo, etc di polinomi
        template < class type , t_dim dimension >
            using t_basis_function = matrix < dimension , t_func < type , dimension > >;

        //puntatore alla base polinomiale
        template < class type , t_dim dimension >
            using t_ptr = std :: shared_ptr < t_basis_function< type , dimension > >;   
            
            
        template < class type , t_dim dimension > class basis ;
        
        template < class type , t_dim dimension = 1 , t_dim codimension = 1 > class polynomial ;
        
        /*template < class type , t_dim dimension , t_dim codimension >
            polynomial<type,dimension,codimension> operator * 
                ( const type factor , const polynomial<type,dimension,codimension> right ); 
                
        template < class type , t_dim dimension , t_dim codimension >
            polynomial<type,dimension,codimension> operator * 
                ( const polynomial<type,dimension,codimension> right , const type factor );
                
        template < class type , t_dim dimension , t_dim codimension >
            polynomial<type,dimension,codimension> operator + 
                ( const polynomial<type,dimension,codimension> right , const type factor );
        
        template < class type , t_dim dimension , t_dim codimension >
            polynomial<type,dimension,codimension>
                operator * 
                    ( const polynomial<type,dimension,codimension> left ,
                         const polynomial<type,dimension,codimension> right ) ;*/
                        
        //monomials        
                                        
        template < class type , t_dim dimension = 1 >
            basis < type , dimension > monomial_basis (const t_deg<dimension> );
            
        template < class type ,  t_dim dimension = 1 > class make_monomials;
        
        
        template < class contained , t_dim dimension >
                using t_list = std :: vector < std :: pair < t_deg < dimension > , contained > >;
                
        template < class type , t_dim dimension , t_dim codimension > 
                using t_coeff_array = matrix < dimension , t_var<type,codimension> >;  
                
        //array dei coefficienti di un polinomio
        //vettore,matrice, cubo, etc di polinomi
        template < class type , t_dim dimension > //, t_dim codimension > 
            using t_coeff_array_basis = t_coeff_array<type,dimension,1>;
            //using t_coeff_array = matrix < dimension , t_var<type,codimension> >;
        
            /*template < class type , t_dim dimension , t_dim codimension>
                using t_list_Bspline = 
                    std :: array < eliastocco :: polynomial :: t_list < type , dimension > , codimension >;*/

            template < class type , t_dim codimension> 
                using control_point = t_var < type,codimension>;
                
                
                  template < class type , t_dim dimension >
                    std :: vector <type> from_array_to_vector ( std :: array <type,dimension> input );
                    
            template < class type , t_dim dimension >
                    std :: vector <type> from_array_to_vector ( std :: array <type,dimension> input ){
                         
                        std :: vector < type > output;
                        std::copy(input.begin(), input.end(), std::back_inserter(output));                            
                        return output;
                    
                    };
                    
            //valutazione del punto
            template < class type , t_dim dimension , t_dim codimension >
                using t_eval = std :: function < t_var<type,codimension> 
                                            ( const t_coeff_array<type,dimension,codimension> & ,
                                                        const t_ptr < type , dimension > , 
                                                            const t_var < type , dimension > ) >;
                                                            
                    //valutazione del punto
            /*template < class type , t_dim dimension >
                using t_eval_basis = t_eval<type,dimension,1>;*/
                
            /*    std :: function < type 
                                            ( const t_coeff_array_basis<type,dimension> & ,
                                                        const t_ptr < type , dimension > & , 
                                                            const t_var < type , dimension > ) >; */
                                                            
                                                            
            template < class type , t_dim dimension , t_dim codimension >
                t_var<type,codimension>  default_polynomial_evaluation 
                       ( const t_coeff_array<type,dimension,codimension> & _coeffs ,
                            const t_ptr < type , dimension > _basis_function , 
                                const t_var < type,dimension > x );
        
                   
        }
        
}
