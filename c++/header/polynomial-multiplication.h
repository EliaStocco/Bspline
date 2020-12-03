
//#pragma once

//namespace eliastocco {

//    namespace polynomial {
    
//        template < class type , t_dim dimension , t_dim codimension >
//            polynomial<type,dimension,codimension>
//                /*polynomial<type,dimension,codimension> ::*/ operator * 
//                    ( const polynomial<type,dimension,codimension> left ,
//                         const polynomial<type,dimension,codimension> right ) 
                        {
                    
                        const unsigned int width = 15 ;
                        const unsigned int prec  = numDecimals ;
                                                
                        
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
                    
                        //ottengo gli array con i coeffcienti (reali)
                        auto left_coeff_real  = left.all_coeffs();
                        auto right_coeff_real = right.all_coeffs();  
                    
                        // eseguo più questo controllo perché l'algoritmo funziona anche se
                        // i coefficienti hanno lunghezza diversa
                        // l'importante è che allo stesso indice corrisponda lo stesso monomio
                    
                        //controllo che gli array abbiano la stessa dimensione
                        /*auto N = left_coeff_real.size();                    
                        if ( N != right_coeff_real.size() ){                            
                            assert ("polynomial :: operator * error; coeff array of different dimension" );                
                        }*/
                    
                        //modifico la dimensione degli array dei coefficienti 
                        // in modo tale che abbiano lunghezza pari ad una potenza di 2
                        unsigned int n = 1;
                        //while ( n < left_coeff_real.size() or n < right_coeff_real.size()){
                        while ( n < left_coeff_real.size() + right_coeff_real.size()){
                            n <<= 1;
                        }
                        
                        // non ne modifico la lunghezza,
                        // tanto i valori qui contenuti verranno copiati in due nuovi array
                    
                        //left_coeff_real.resize(n); // i nuovi valori saranno impostati a 0
                        //right_coeff_real.resize(n);// i nuovi valori saranno impostati a 0
                    
                        // efinisco gli array di coefficienti complessi
                        // che saranno gli input della funzione fft
                        typename Rosetta_Code :: CArray<type> left_coeff  (n);
                        typename Rosetta_Code :: CArray<type> right_coeff (n);

                        // copio i coefficienti e li trasformo in valori complessi
                        //left
                        for ( auto i = 0 ; i < left_coeff_real.size(); ++i ){
                            left_coeff[i] = Rosetta_Code :: Complex <type> ( left_coeff_real[i] , 0 );
                        }
                        //right
                        for ( auto i = 0 ; i < right_coeff_real.size(); ++i ){
                            right_coeff[i] = Rosetta_Code :: Complex <type> ( right_coeff_real[i] , 0 );
                        }
                        
                        if ( debug ) {   
                            std :: cout << "original arrays:" << std::setw(width)<< std :: endl;
                            for ( auto i = 0 ; i < n; ++i ){
                                std :: cout 
                                    <<left_coeff[i]  << std::setw(width)<< std::setprecision(prec)
                                    << " , " << std::setprecision(prec)
                                    <<right_coeff[i] << std::setw(width)<< std::setprecision(prec)
                                    << std :: endl ;                                     
                            }                        
                        }
                        
                        // eseguo le trasformate di Fourier dirette inplace  
                        Rosetta_Code :: fft ( left_coeff  );
                        Rosetta_Code :: fft ( right_coeff );         
                        
                        if ( debug ) {   
                            std :: cout << "fft arrays:" << std::setw(width)<< std :: endl;
                            for ( auto i = 0 ; i < n; ++i ){
                                std :: cout
                                    <<left_coeff[i]  << std::setw(width)<< std::setprecision(prec)
                                    << " , " << std::setprecision(prec)
                                    <<right_coeff[i] << std::setw(width)<< std::setprecision(prec)<< std :: endl ;                                     
                            }                        
                        }
                    
                        // moltiplico tra di loro i due vettori dei coefficienti nello spazio reciproco
                        // sfrutto il fatto che gli array di coefficienti siano dei std :: valarray
                        Rosetta_Code :: CArray<type> coeff = left_coeff * right_coeff;                    
                          
                        if ( debug ) {   
                            std :: cout << "fft multiplication:" << std::setw(width)<< std :: endl;
                            for ( auto i = 0 ; i < n; ++i ){
                                std :: cout
                                    <<coeff[i]  << std::setw(width)<< std::setprecision(prec)
                                    << std :: endl ;                                     
                            }                        
                        }
                    
                    
                        // eseguo la trasformata di Fourier inversa inplace                        
                        Rosetta_Code :: ifft( coeff );
                        
                        if ( debug ) {   
                            std :: cout << "new coeffs:" << std::setw(width) << std :: endl;
                            for ( auto i = 0 ; i < n; ++i ){
                                std :: cout 
                                    <<coeff[i]  << std::setw(width)<< std::setprecision(prec)
                                    << std :: endl ;                                     
                            }                        
                        }

                        // creo un polinomio e gli assegno i coefficienti calcolati
                        // solo se questi non sono nulli, 
                        // per ottimizzare spazio e rendere l'output più leggibile
                        polynomial<type,dimension,codimension> output ( basis ) ;
                        
                        
                        auto my_round = [=](const type& value ) {
                            const type multiplier = static_cast<type>(std::pow( 10.0, numDecimals ));
                            return static_cast<type>(std::ceil(value * multiplier) / multiplier );
                        };
                        
                        type c1 = type();
                        type c  = type();
                        for (auto i = 0; i < n; i++){
                            c1 = coeff [i] . real();
                            c = my_round(c1);
                            if ( c != type() ){                                
                                output.set(i,c,check);
                            }
                        }
                    
                        return output;
            
                } ; 
    
//    }
    
//}
