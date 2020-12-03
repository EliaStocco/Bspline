//%%file header/metropolis.h
#pragma once

#include <array>
#include <vector>
#include <typeinfo>
#include <functional>
#include <iostream>
#include <random>

//#include "Metropolis/LCG.h"
//#include "Metropolis/PrintProgress.h" // includendola qui, me la compila, ma essendo static è visibile solo qui dentro

namespace eliastocco {

      template < class Type , class Coordinates /*, unsigned int N_Coordinates*/ >

      class metropolis {

            //typedef std::array < Type , N_Coordinates > Coordinates ;

            typedef std :: function < Type ( const Coordinates ) > TypeDistribution ;
            typedef std :: function < Type ( void ) >TypeDistributionDivision ;
            typedef std :: function < Coordinates ( const Coordinates & ) > TypeStepDistribution ;
            typedef std :: function < Type ( const Coordinates & , const Coordinates & ) > TypeTransitionProbability ;
            typedef std :: function < void ( const Coordinates & ) >TypeSingleStepFunction ;
            typedef std :: function < Type ( void ) > Generator_t;

            static void print_progress (double percentage) {

                #define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
                #define PBWIDTH 60

                int val = (int) (percentage * 100);
                int lpad = (int) (percentage * PBWIDTH);
                int rpad = PBWIDTH - lpad;
                printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
                fflush (stdout);

                #undef PBSTR
                #undef PBWIDTH
            }

            //Methods

            public :
                metropolis () 
                    : metropolis(
                    [](void){    
                        static std :: mt19937 a ;
                        static double N = static_cast<double>(a.max());    
                        return static_cast<double>(a())/N;
                    }){};

            public : 
              metropolis ( const Generator_t Gen ) :
                N_Step  ( 0 ) ,
                Accepted  ( 0 ) ,
                Distribution  ( [] ( const Coordinates a ) { Type c = Type () ; return c ; } ) ,
                DistributionDivision ( [] ( void ) { Type a = Type () ; return a ; } ) ,
                StepDistribution  ([] ( const Coordinates a ) {  return a ; } ) ,
                TransitionProbability( [] ( const Coordinates a , const Coordinates b ) { Type c = Type (); return c ; } ) ,
                Point  () ,
                Generatore  ( Gen ) ,
                ClassSteps  () ,
                ClassPrintProgress  ( false ) ,
                ClassSymmetricTransition ( false ) ,
                ClassIsDistributionDivisionKnown ( false ) ,
                ClassLastSteps       ( ) ,
                Last_N_Step( 0 ) ,
                Last_Accepted( 0 ) ,
                ClassSaveSteps( false ) ,
                ClassSingleStepFunction ( [] ( const Coordinates & a ) { return ; } ) ,
                ClassCallSingleStepFunction ( false ) ,
                ClassSaveDistributionValue ( false ) ,
                ClassDistributionValue ( ) ,
                single_accepted ( false )
                {};


            public : double get_acceptance ( void ) const {
                  return static_cast<double> ( Last_Accepted + Accepted ) / static_cast<double> ( Last_N_Step + N_Step ) ;
            };

            public : unsigned int get_accepted ( void ) const {return Last_Accepted + Accepted ; };

            public : double get_last_acceptance ( void ) const {
                  return static_cast<double> ( Last_Accepted ) / static_cast<double> ( Last_N_Step ) ;
            };

            public : unsigned int get_last_accepted ( void ) const {
                  return Last_Accepted ;
            } ;

            public : void set_point ( Coordinates NewPoint ) {
                  Point = NewPoint ;
                  return ;
            }

            public : void set_target_distribution ( TypeDistribution NewDistribution ){
                  Distribution = NewDistribution ;
                  return ;
            } ;

            public : void set_step_distribution ( TypeStepDistribution NewDistribution ) {
                  StepDistribution = NewDistribution ;
                  return ;
            };

            public : void set_transition ( TypeTransitionProbability NewTransition ){
                  TransitionProbability = NewTransition ;
                  return ;
            } ;

            public : void set_single_step_function ( TypeSingleStepFunction Function ) {
                  call_single_step_function ( true ) ;
                  ClassSingleStepFunction = Function ;
                  return ;
            }; // funzione alla fine di ogni step

            public : void call_single_step_function ( const bool Call ){
                  ClassCallSingleStepFunction = Call ;
                  return ;
            } ;

            public : void set_target_distribution_division( TypeDistributionDivision NewDivision ){
                  DistributionDivision = NewDivision ;
                  return ;
            } ;

            public : bool is_transition_symmetric ( const bool val ){
                  ClassSymmetricTransition = val ;
                  return true ;
            } ;

            public : bool use_distribution_division ( const bool val ){
                  ClassIsDistributionDivisionKnown = val ;
                  return true ;
            } ;

            public : bool save_distribution_value( const bool Save ) {
                  ClassSaveDistributionValue = Save ;
                  return true ;};

            public : bool have_accepted( void ) {return single_accepted ;};

            public : Type get_distribution_value ( void ) {return ClassDistributionValue ;};

            public : void clear ( void ) {
                  Accepted = 0 ;
                  N_Step = 0 ;
                  Last_Accepted = 0 ;
                  Last_N_Step = 0 ;
                  ClassSteps . clear () ;
                  ClassLastSteps . clear () ;
                  return ;
            };

            public : Coordinates get_point ( void ) const {return Point ;};

            public :
            Coordinates move ( void ) {
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

                        Q = A * B /
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

            };

            public :
            std :: vector < Coordinates > run ( const unsigned int N_Step ) {
                  ClassSteps . insert ( ClassSteps . begin () , 
                                       ClassLastSteps . begin () , 
                                       ClassLastSteps . end () ) ;

                  if ( ClassSaveSteps ) {ClassLastSteps . resize ( N_Step ) ;}
                  //
                  N_Step += Last_N_Step ;
                  Accepted += Last_Accepted ;
                  Last_N_Step = 0 ;
                  Last_Accepted = 0 ;
                  //
                  if ( ClassSaveSteps && ClassPrintProgress ) {

                        for ( unsigned int i = 0 ; i < N_Step ; i ++ ) {
                              ClassLastSteps [ i ] = move () ;
                              print_progress ( i / static_cast < Type > ( N_Step ) ) ;
                        }

                        print_progress ( 1 ) ;
                        std :: cout << std :: endl ;
                  } else if ( ClassSaveSteps && !ClassPrintProgress ) {
                        for ( unsigned int i = 0 ; i < N_Step ; i ++ ) {
                              ClassLastSteps [ i ] = move () ;
                        }

                  } else if ( !ClassSaveSteps && ClassPrintProgress ) {
                        for ( unsigned int i = 0 ; i < N_Step ; i ++ ) {
                              move () ;
                              print_progress ( i / static_cast < Type > ( N_Step ) ) ; }

                        print_progress ( 1 ) ;
                        std :: cout << std :: endl ;
                  } else if ( !ClassSaveSteps && !ClassPrintProgress ) {
                        for ( unsigned int i = 0 ; i < N_Step ; i ++ ) {move () ;}
                  }

                  return ClassLastSteps ;
            };

            public :
            std :: vector < Coordinates > get_steps ( void ) const {return ClassSteps ;};

            public :
            std :: vector < Coordinates > get_last_steps ( void ) const {return ClassLastSteps ;};

            public :
            bool show_progress ( const bool Print ) {
                  ClassPrintProgress = Print ;
                  return true ;
            };

            public :
            bool save_steps ( const bool Save ) {
                  ClassSaveSteps = Save ;
                  return true ;
            };

            public :
            bool get_show_progress ( void ) const {
                  return ClassPrintProgress ;
            } ;

            public : unsigned int get_n_steps ( void ) const { return Last_N_Step + N_Step ; };


            //Members

            protected :
                unsigned int N_Step ;
            protected :
                unsigned int Accepted ;
            protected :
                TypeDistribution Distribution ; // distribuzione da generare ( può anche non essere normalizzata)
            protected :
                TypeDistributionDivision DistributionDivision ;
            protected ://funzione che mi fa fare un salto secondo la distribuzione specificato sotto        
                TypeStepDistribution StepDistribution ; 
            protected ://distribuzione dello step (deve essere normalizzata )
                TypeTransitionProbability TransitionProbability ; 
            protected :
                Coordinates Point ;
            protected :
                Generator_t Generatore ;
            protected :
                std :: vector < Coordinates > ClassSteps ;
            protected :
                bool ClassPrintProgress;
            protected :
                bool ClassSymmetricTransition ;
            protected :
                bool ClassIsDistributionDivisionKnown ;
            protected :
                std :: vector < Coordinates > ClassLastSteps ;
            protected :
                unsigned int Last_N_Step ;
            protected :
                unsigned int Last_Accepted ;
            protected :
                bool ClassSaveSteps ;
            protected :
                TypeSingleStepFunction ClassSingleStepFunction ;
            protected :
                bool ClassCallSingleStepFunction ;
            protected :
                bool ClassSaveDistributionValue ;
            protected :
                Type ClassDistributionValue ;
            protected :
                bool single_accepted ;

      } ;

}