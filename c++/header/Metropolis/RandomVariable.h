#pragma once 

namespace SingWolf {

	template < class OutPut > 
	
		class RandomVariable {
		
			public :
				virtual OutPut operator ()  ( void ) = 0 ;

			public :
				virtual void set_generator ( SingWolf :: RandomVariable < OutPut > & ) = 0 ; 
	
		} ;
	
}