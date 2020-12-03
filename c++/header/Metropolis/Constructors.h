
#define TypeDistribution std :: function < Type ( const Coordinates ) >

#define TypeDistributionDivision std :: function < Type ( void ) >
	
#define TypeStepDistribution std :: function < Coordinates ( const Coordinates & ) >

#define TypeTransitionProbability std :: function < Type ( const Coordinates & , const Coordinates & ) >


template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	Metropolis < Type , Coordinates > :: Metropolis ( ) {

				
	}
	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	Metropolis < Type , Coordinates > :: Metropolis ( SingWolf :: LCG < Type > & Gen ) : 

		N_Step 				 ( 0 ) ,

		Accepted 			 ( 0 ) ,

		Distribution 		 ( 	[] ( const Coordinates a ) 							{ Type c = Type () ; return c ; } 	) ,

		DistributionDivision ( 	[] ( void ) 										{ Type a = Type () ; return a ; } 	) ,

		StepDistribution 	 (	[] ( const Coordinates a ) 							{ 					 return a ; } 	) ,

		TransitionProbability( [] ( const Coordinates a , const Coordinates b ) 	{ Type c = Type (); return c ; } 	) ,

		Point 				 () ,

		Generatore 			 ( Gen ) ,

		ClassSteps 			 () ,

		ClassPrintProgress 	 ( false ) ,

		ClassSymmetricTransition ( false ) ,

		ClassIsDistributionDivisionKnown ( false ) ,

		ClassLastSteps       	( ) ,

		Last_N_Step				( 0 ) ,

		Last_Accepted			( 0 ) ,

		ClassSaveSteps			( false ) ,

		ClassSingleStepFunction ( [] ( const Coordinates & a ) { return ; } ) ,

		ClassCallSingleStepFunction ( false ) ,  

		ClassSaveDistributionValue ( false ) ,

		ClassDistributionValue 	( ) ,

		single_accepted 		( false )

		{
				
		
	}

#undef TypeDistribution

#undef TypeDistributionDivision
	
#undef TypeStepDistribution

#undef TypeTransitionProbability
