
#define TypeDistribution std :: function < Type ( const Coordinates ) >

#define TypeDistributionDivision std :: function < Type ( void ) >
	
#define TypeStepDistribution std :: function < Coordinates ( const Coordinates & ) >

#define TypeTransitionProbability std :: function < Type ( const Coordinates & , const Coordinates & ) >


	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	void Metropolis < Type , Coordinates > :: set_point ( Coordinates NewPoint ) {
		
		Point = NewPoint ; 
		
		return ;
		
	}
	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	void Metropolis < Type , Coordinates > :: set_step_distribution ( TypeStepDistribution NewDistribution ) {
		
		StepDistribution = NewDistribution ; 
		
		return ;
		
	}
	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	void Metropolis < Type , Coordinates > :: set_target_distribution ( TypeDistribution NewDistribution ) {
		
		Distribution = NewDistribution ; 
		
		return ;
		
	}
	

template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	void Metropolis < Type , Coordinates > :: set_target_distribution_division ( TypeDistributionDivision NewDivision ) {
		
		DistributionDivision = NewDivision ; 
		
		return ;
		
	}

template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	void Metropolis < Type , Coordinates > :: set_transition ( TypeTransitionProbability NewTransition ) {
		
		TransitionProbability = NewTransition ; 
		
		return ;
		
	}


template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	bool Metropolis < Type , Coordinates > :: use_distribution_division ( const bool val ) {
		
		ClassIsDistributionDivisionKnown = val ; 
		
		return true ;
		
	}

	

template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	bool Metropolis < Type , Coordinates > :: is_transition_symmetric ( const bool val ) {
		
		ClassSymmetricTransition = val ; 
		
		return true ;
		
	}



template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	Coordinates Metropolis < Type , Coordinates > :: get_point ( void ) const {
		
		return Point ; 
		
	}
	

template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	unsigned int Metropolis < Type , Coordinates > :: get_n_steps ( void ) const {
		
		return Last_N_Step + N_Step ; 
		
	}

template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	void Metropolis < Type , Coordinates > :: clear ( void ) {
		
		Accepted = 0 ;
		
		N_Step = 0 ; 

		Last_Accepted = 0 ;
		
		Last_N_Step = 0 ; 

		ClassSteps . clear () ;

		ClassLastSteps . clear () ;
		
		return ; 
		
	}
	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	unsigned int Metropolis < Type , Coordinates > :: get_accepted ( void ) const {
			
		return Last_Accepted + Accepted ; 
		
	}
	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	double Metropolis < Type , Coordinates > :: get_acceptance ( void ) const {
			
		return static_cast<double> ( Last_Accepted + Accepted ) / static_cast<double> ( Last_N_Step + N_Step ) ; 
		
	}

template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	unsigned int Metropolis < Type , Coordinates > :: get_last_accepted ( void ) const {
			
		return Last_Accepted ; 
		
	}
	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	double Metropolis < Type , Coordinates > :: get_last_acceptance ( void ) const {
			
		return static_cast<double> ( Last_Accepted ) / static_cast<double> ( Last_N_Step ) ; 
		
	}

	
template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	Coordinates Metropolis < Type , Coordinates > :: move ( void ) {
		
		
		Coordinates NewPoint = StepDistribution ( Point ) ; // genero nuovo punto

						

		Type Q = Type () ;	

		Type A = Type () ; 

		Type B = Type () ; 

		Type C = Type () ; 

		Type D = Type () ; 

		

		if ( ClassIsDistributionDivisionKnown ) {

			Q = DistributionDivision () ;

		} else if ( ClassSymmetricTransition ) {


			B = Distribution ( NewPoint ) ;

			if ( ClassSaveDistributionValue && get_n_steps () > 0 ) {

				D = ClassDistributionValue ;

			} else {

				D = Distribution (   Point  )  ;

			}
			

			Q = B / D ;		
						
		} else {

			A = TransitionProbability ( Point , NewPoint ) ;

			B = Distribution ( NewPoint ) ;

			C = TransitionProbability ( NewPoint , Point ) ;


			if ( ClassSaveDistributionValue && get_n_steps () > 0 ) {

				D = ClassDistributionValue ;

			} else {

				D = Distribution (   Point  )  ;

			}			

			Q = 	A * B / 
		
						( C * D  ) ;

		}		
								
		
						
		if ( Q >= 1.0 ) {
			
				Point = NewPoint ; 
				
				Last_Accepted ++ ;	

				single_accepted = true ; 
		
		} else {
			
			if ( Generatore ( ) < Q ) {
				
				Point = NewPoint ; 	
				
				Last_Accepted ++ ;	

				single_accepted = true ;			
				
			} else {

				single_accepted = false ;			

			}
			
		}		
		
		Last_N_Step ++ ; 

		if ( ClassCallSingleStepFunction ) {

			ClassSingleStepFunction ( Point ) ;

		}

		if ( ClassSaveDistributionValue && single_accepted) {

			ClassDistributionValue = B ;

		} // altrimenti tengo quella vecchia 

	
		return Point ;
		
	}
	

template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

	std :: vector < Coordinates > Metropolis < Type , Coordinates > :: run ( const unsigned int Step ) {

		//

		ClassSteps . insert ( ClassSteps . begin () , ClassLastSteps . begin () , ClassLastSteps . end () ) ; 

		if ( ClassSaveSteps ) {	ClassLastSteps . resize ( Step ) ;	}

		//

		N_Step += Last_N_Step ;

		Accepted += Last_Accepted ;

		Last_N_Step = 0 ; 

		Last_Accepted = 0 ; 

		//

		if ( ClassSaveSteps && ClassPrintProgress ) {

			for ( unsigned int i = 0 ; i < Step ; i ++ ) {

				ClassLastSteps [ i ] = move () ;

				print_progress ( i / static_cast < Type > ( Step ) ) ; 

			}

			print_progress ( 1 ) ; 

			std :: cout << std :: endl ; 

		} else if ( ClassSaveSteps && !ClassPrintProgress ) {

			for ( unsigned int i = 0 ; i < Step ; i ++ ) {

				ClassLastSteps [ i ] = move () ;				

			}

		} else if ( !ClassSaveSteps && ClassPrintProgress ) {

			for ( unsigned int i = 0 ; i < Step ; i ++ ) {

				move () ;

				print_progress ( i / static_cast < Type > ( Step ) ) ; 

			}

			print_progress ( 1 ) ; 

			std :: cout << std :: endl ; 


		} else if ( !ClassSaveSteps && !ClassPrintProgress ) {

			for ( unsigned int i = 0 ; i < Step ; i ++ ) {	move () ;	}

		}

		return ClassLastSteps ;
		
	}


	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		std :: vector < Coordinates > Metropolis < Type , Coordinates > :: get_last_steps ( void ) const {

			return ClassLastSteps ;

		}

	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		std :: vector < Coordinates > Metropolis < Type , Coordinates > :: get_steps ( void ) const {

			return ClassSteps ;

		}


	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		bool Metropolis < Type , Coordinates > ::  save_steps ( const bool Save ) {

				ClassSaveSteps = Save ;

				return true ;
			}

	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		bool Metropolis < Type , Coordinates > ::  show_progress ( const bool Print ) {

				ClassPrintProgress = Print ;

				return true ;
			}


	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		bool Metropolis < Type , Coordinates > ::  get_show_progress ( void ) const {

				return ClassPrintProgress ;
			}



	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		void Metropolis < Type , Coordinates > :: set_single_step_function ( TypeSingleStepFunction Function ) {

				call_single_step_function ( true ) ;
			
				ClassSingleStepFunction = Function ; 
		
				return ;
		
			}


	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		void Metropolis < Type , Coordinates > ::  call_single_step_function ( const bool Call ) {

				ClassCallSingleStepFunction = Call ;

				return ;

			}

	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		bool Metropolis < Type , Coordinates > ::  have_accepted ( void ) {

				return single_accepted ;
			}

	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		Type Metropolis < Type , Coordinates > ::  get_distribution_value ( void ) {

				return ClassDistributionValue ;

			}

	template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

		bool Metropolis < Type , Coordinates > ::  save_distribution_value ( const bool Save ) {

				ClassSaveDistributionValue = Save ; 

				return true ;

			}

	



#undef TypeDistribution

#undef TypeDistributionDivision
	
#undef TypeStepDistribution

#undef TypeTransitionProbability