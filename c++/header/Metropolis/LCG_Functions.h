


template < class Type /*, unsigned int N*/ > 
int LCG< Type /*, N*/ >::LoadPrimes ( std::string File ) {

	int Ritorno = 1 ;
	std::ifstream Primes ;
	Primes.open ( File ) ;
	
	if (Primes.is_open()){
      Primes >> n3 >> n4 ;
   } else {
	   std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
	   
	   Ritorno = -1 ;  
	   
   }   
   
   Primes.close();
   
   return Ritorno ; 
   
   return 1 ; 
   
}

template < class Type /*, unsigned int N*/ > 
int LCG< Type /*, N*/ > :: LoadSeed ( std::string File ) {
	
	
	int Ritorno = 1 ;
	int seed[4];
   std::ifstream input ;
   input.open ( File ) ;   
   std::string property;
   
   if (input.is_open()){
      while ( !input.eof() ){ // in realtà ho solo una riga
         input >> property;
         if( property == KEY_SEED ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetLCG(seed);
         }
      }
      input.close();
   } else {
	   
	   std::cerr << "PROBLEM: Unable to open "<<File<< std::endl;
	   
	   Ritorno = -1 ;
	   
   }
		
	return Ritorno ; 
}

template < class Type /*, unsigned int N*/ >
int LCG< Type /*, N*/ > :: SaveSeed( std::string File ){
	
	int Ritorno = 1 ;
   std::ofstream WriteSeed;
   WriteSeed.open(File);
   
   if (WriteSeed.is_open()){
	   
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << std::endl;
	  
   } else {
   
	   std::cerr << "PROBLEM: Unable to open "<<File<< std::endl;
	   
	   Ritorno = -1 ;
	   
   }
   
  WriteSeed.close();
  
  return Ritorno ;
}



template < class Type /*, unsigned int N*/ > 
void LCG< Type /*, N*/ > :: SetLCG(int s[4]){ 

  l1 = s[0]%NUM;
  l2 = s[1]%NUM;
  l3 = s[2]%NUM;
  l4 = s[3]%NUM;
  l4 = 2*(l4/2)+1;
  
  return ;
  
}

template < class Type /*, unsigned int N*/ > 	
void LCG< Type /*, N*/ > :: SetLCG(int s[4], int p1, int p2){
	
	SetLCG ( s ) ;	

	n3 = p1;
	n4 = p2;

  return;
}
