#pragma once


#include <string>
#include <array>
#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>
#include <tuple>
#include <fstream>
#include <type_traits>



#ifndef M_PI
#define M_PI acos (-1)
#endif

#include "vector\\vector.h"
#include "array\\array.h"
//#include "array\\singwolf_array.h"
#include "SphericalCoordinates\\SphericalCoordinates.h"


#define Print(a) std::cout<<#a<<" = "<<a<<std::endl


// ci sono da ridefinire i costruttori ,
// l'operatore di assegnamento = sia per std che per SingWolf
//
 
namespace SingWolf {

	template < class T , class ... Args >
		class tuple : public std :: tuple < T , Args ... > {

			public :
				tuple (  ) {} ;
				tuple ( const int a ) {} ;

	} ;

}

namespace SingWolf {

	template < class T , long unsigned int  N >
		class array : public std :: array < T , N > {

			public :
				array (  ) {} ;
				array ( const int a ) {} ;
				array (const SingWolf :: array<T,N> & copia){
					for ( unsigned int i = 0 ; i < N ; i ++){
						(*this)[i] = copia[i] ;
					}
				}

				array (const std :: array<T,N> & copia){
					for ( unsigned int i = 0 ; i < N ; i ++){
						(*this)[i] = copia[i] ;
					}
				}

			public :
				SingWolf::array<T,N>& operator=(SingWolf :: array<T,N> copia){
					swap(copia);
					return *this;
				}
				SingWolf::array<T,N>& operator=(std :: array<T,N> copia){
					swap(copia);
					return *this;
				}
			public :
				void swap (SingWolf :: array<T,N> & copia){

					for (unsigned int i = 0 ; i < N ; i ++ ){
						std :: swap <T> ((*this)[i],copia[i]) ;
					}
					return;
				}
				void swap (std :: array<T,N> & copia){

					for (unsigned int i = 0 ; i < N ; i ++ ){
						std :: swap <T> ((*this)[i],copia[i]) ;
					}
					return;
				}

	} ;

}
	

namespace SingWolf {

template <size_t n, typename... T>
typename std::enable_if<(n >= sizeof...(T))>::type
    print_tuple(std::ostream&, const std::tuple<T...>&)
{}

template <size_t n, typename... T>
typename std::enable_if<(n < sizeof...(T))>::type
    print_tuple(std::ostream& os, const std::tuple<T...>& tup)
{
    if (n != 0)
        os << "| ";
    os << std::get<n>(tup);
    print_tuple<n+1>(os, tup);
}

template <typename... T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& tup)
{
    //os << "[";
    print_tuple<0>(os, tup);
    return os ; // << "]";
}

}

namespace SingWolf  {

	template < bool first_line , class Type > 	
		bool save_on_file ( const std :: string File , const std :: string title , const Type & vett ) ;

	template < bool first_line , class Type > 
	
		bool save_on_file ( const std :: string File , const std :: string title , const Type & vett ){

		std :: ofstream out ( File ) ;

		if ( title . size () > 0 ) {

			out << title << std :: endl ;

		}

		const auto Size = vett . size () ; 

		//std :: cout << "ok , size : " << Size << std :: endl ; 

		for ( unsigned int i = 0 ; i < Size ; i ++ ) {

			//std :: cout << i << std :: endl ; 

			if /*constexpr*/ ( first_line ){

				out << i + 1 << "," ;

			}

			out << vett [ i ] << std :: endl ;

		}

		//std :: cout << "chiudo" << std :: endl ; 

		out . close () ;

		//std :: cout << "esco" << std :: endl ; 

		return true ;

	}


}

namespace SingWolf {

	template < class Type_In , class Type_Out > 
		Type_Out transpose_in_array ( const Type_In & matrix ) ;


	template < class Type_In , class Type_Out > 

		Type_Out transpose_in_array ( const Type_In & matrix ) {

			Type_Out ritorno ; 

			/*define_has_member(resize);

			if constexpr ( has_member ( Type_Out , resize ) ) {

				ritorno . resize ( matrix [ 0 ] . size () ) ;

			}*/

			const unsigned int Row = matrix . size () ;

			const unsigned int Col = matrix [ 0 ] . size () ; 

			for ( unsigned int i = 0; i < Col ; i ++ ) {

				ritorno [ i ] . resize ( Row ) ;

				for ( unsigned int j = 0 ; j < Row ; j ++ ) {

					ritorno [ i ] [ j ] = matrix [ j ] [ i ] ;

				}
				
			}

			return ritorno ;

		}
}