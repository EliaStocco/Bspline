/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#pragma once


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <array> 
#include <type_traits> 

#include "RandomVariable.h"
#include "Macro_Elia.h"


#ifndef M_PI
#define M_PI   acos(-1)
#endif


#define NUM 4096
const std::string SEED_OUT 	= "seed.out" 	;
const std::string SEED_IN 	= "seed.in" 	;
const std::string PRIMES 	= "Primes" 		;
const std::string KEY_SEED 	= "RANDOMSEED" 	;


static void SetFilesName ( void ) ;

namespace SingWolf {

static std::string dirnameOf(const std::string& fname) ;

std::string dirnameOf(const std::string& fname)
{
     size_t pos = fname.find_last_of("\\/");
     return (std::string::npos == pos)
         ? ""
         : fname.substr(0, pos);
}

}

static class StartBeforeMain {
	
	public :
		StartBeforeMain () {	SetFilesName () ; }
	
} StartBeforeMain_Variable ; 

void SetFilesName ( void ) {
	
	//std::cout<<"Hello 1"<<std::endl ;
	
	const_cast<std::string &>(SEED_OUT) 	= SingWolf :: dirnameOf(std::string(__FILE__)) + std::string("\\") + SEED_OUT 	;
	
	const_cast<std::string &>(SEED_IN)	 	= SingWolf :: dirnameOf(std::string(__FILE__)) + std::string("\\") + SEED_IN 	;
	
	const_cast<std::string &>(PRIMES) 		= SingWolf :: dirnameOf(std::string(__FILE__)) + std::string("\\") + PRIMES 	;

	//std::cout<<SEED_OUT<<std::endl ;

	//std::cout<<SEED_IN<<std::endl ;

	//std::cout<<PRIMES<<std::endl ;
	
	return ;	
	
}



namespace SingWolf {

template < class Type /* unsigned int N */ >

	class LCG : 

		virtual public RandomVariable< Type > {
	
	
		//satic_assert ( N != 0 , "Wrong Dimension" ) ; // errore se N == 0 
		
				

		protected:

			int l1;
			int l2;
			int l3;
			int l4;
			int n3;
			int n4;
  
			const int m1 = 502;
			const int m2 = 1521;
			const int m3 = 4071;
			const int m4 = 2107;
			const int n1 = 0;
			const int n2 = 0;



		public:

			LCG ( SingWolf :: RandomVariable < Type > & Gen ) : ClassGenerator (Gen) , ClassGenerator_Used ( true ) {

				LoadPrimes () ;  

				LoadSeed () ; 

			}
  
			LCG() : ClassGenerator (*this) , ClassGenerator_Used ( false ) { 

				LoadPrimes () ;  

				LoadSeed () ; 

			};// Elia
  
			LCG( std::string Primes_File , std::string Seed_File ) : ClassGenerator (*this) , ClassGenerator_Used ( false ) { 

				LoadPrimes (Primes_File) ;  

				LoadSeed (Seed_File) ; 

			};// Elia
  

			~LCG() { };
			

			void SetLCG(int [] , int, int);
			
			void SetLCG(int [] );
  
  
  
			inline int LoadPrimes( void ) 				{ 	return LoadPrimes 	( PRIMES.c_str() ) 	;	}
			inline int LoadPrimes( std::string File ) 	;			
			
  
			inline int LoadSeed( void ) 				{ 	return LoadSeed 	( SEED_IN.c_str() ) ;	}
			inline int LoadSeed( std::string File ) 	;			
			
  
			inline int SaveSeed( void ) 				{ 	return SaveSeed 	( SEED_OUT.c_str())	;	}
			inline int SaveSeed( std::string File ) 	;
			
  

		public :
		
			Type operator () (void) { 

				if ( ClassGenerator_Used ) { return ClassGenerator () ; }
		
				#define	twom12 0.000244140625 // ottimiziamo
  
				int i1,i2,i3,i4;
  

				i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
				i2 = l2*m4 + l3*m3 + l4*m2 + n2;
				i3 = l3*m4 + l4*m3 + n3;
				i4 = l4*m4 + n4;
				l4 = i4%NUM;
				i3 = i3 + i4/NUM;
				l3 = i3%NUM;
				i2 = i2 + i3/NUM;
				l2 = i2%NUM;
				l1 = (i1 + i2/NUM)%NUM;
				return twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4)))); // ottimiziamo
  
				#undef twom12
  
			} 
	

		public :
			inline Type generate( const Type min , const Type max ){  return min + ( max - min ) * (*this) ( ) ;	}

		public :
			void set_generator ( SingWolf :: RandomVariable < Type > & Gen ) { 

				ClassGenerator = Gen ;

				ClassGenerator_Used = true ;

			} ; 

		protected :
			SingWolf :: RandomVariable < Type > & ClassGenerator ; 

		protected :
			bool ClassGenerator_Used ; 


  
  
	} ;


#include "LCG_Functions.h"

}
