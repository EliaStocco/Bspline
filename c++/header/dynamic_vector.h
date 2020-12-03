//%%file header/dynamic_vector.h
#pragma once 

#include <vector>
#include <array>
//#include <initializer_list>

namespace eliastocco { // you can change namespace    
       
    namespace container {
    
        namespace dynamic {
        
            typedef /*long unsigned*/ int t_dim;
            typedef /*long*/ long unsigned int ut_dim;
        
            constexpr ut_dim decrease ( ut_dim a ) { return (a > 1 ) ? ( a - 1 ) : 0 ; }
        
            // do not change the order of these classes

            // do not use this class !!         
            template < class type , template <class...> class STL_container >
                class my_dynamic_vector : public STL_container < type > {

                    public :
                        my_dynamic_vector ( std :: array < t_dim , 1 > Lenght = { 0 } ) 
                            : STL_container < type > ( Lenght [ 0 ] ) {};

                };

            //eredito un "array" di dimensione Dimension -1 
            template < ut_dim Dimension , class type , 
                template <class...> class STL_container = std :: vector > 
                    class dynamic_vector 
                        : public my_dynamic_vector 
                            < dynamic_vector < decrease(Dimension) , type , STL_container > , STL_container > {
                        
                        #define thisclass dynamic_vector<Dimension,type,STL_container>
                        
                        #define motherclass my_dynamic_vector<dynamic_vector<decrease(Dimension),\
                        type,STL_container>,STL_container>
                        
                        #define childclass dynamic_vector < decrease(Dimension) , type , STL_container >
                        
                        static_assert ( Dimension > 0 , "Dimension should be > 0" ) ;
                        
                        public :
                            typedef std :: array < t_dim , Dimension > Index_t;
                        
                        
                        public :
                            dynamic_vector ( Index_t Lenght = { 0 } ) {

                                std :: array < t_dim , decrease(Dimension) > Other = { 0 } ;

                                for ( auto i = 0 ; i != Other . size () ; i ++ ) {
                                    Other [ i ] = Lenght [ i + 1 ] ;
                                }

                                my_dynamic_vector < 
                                    dynamic_vector < decrease(Dimension) , type , STL_container > ,
                                        STL_container > ::
                                resize ( Lenght [ 0 ], dynamic_vector < 
                                                decrease(Dimension) , type , STL_container > ( Other ) ) ;

                            };
                            
                //*************************************************************** 
                
                      /*  public : 
                            bool set ( const std :: array<t_dim,Dimension> index , const type value ){
                                (*this)[index]=value;
                                return true;
                            }*/
                
                         public :
                            childclass & operator[](const t_dim index ){
                                //std :: cout << "childclass & operator[](const t_dim index ) " << index << std :: endl;
                                //return motherclass::operator[](index);
                                return motherclass::at(index);
                            }
                        public :
                            childclass operator[](const t_dim index ) const {
                                //std :: cout << "childclass operator[](const t_dim index ) const" << index << std :: endl;
                                //return motherclass::operator[](index);
                                return motherclass::at(index);
                            }
                        
                //***************************************************************                
                
                        public :
                            template < ut_dim Dim >
                            auto & operator[](const std :: array<t_dim,Dim> index ) {
                                //std :: cout << "auto & operator[](const std :: array<t_dim,Dim> index ) Dim=" << Dim 
                                //<< std :: endl;
                                if constexpr (Dim>1){
                                    
                                    std :: array < t_dim , decrease(Dim) > new_index = { 0 } ;
                                    for ( auto i = 0 ; i < Dim-1 ; i ++ ) {
                                        new_index [ i ] = index [ i + 1 ] ;
                                    }
                                    
                                    //std :: cout << "ciao " << index[0] << std :: endl;
                                    //std :: cout << "ciao " << index[1] << std :: endl;
                                    //std :: cout << "ciao " << index . size () << std :: endl;
                                    //std :: cout << "ciao " << new_index[0] << std :: endl;
                                    
                                    //auto a = motherclass::operator[](index[0]);

                                    //return motherclass::operator[](index[0])[new_index];
                                    return (motherclass::at(index[0]))[new_index];
                                
                                }else {
                                    //std :: cout << "ciao2" << index[0] << std :: endl;
                                    return motherclass::at(index[0]);
                                }
                                
                            }
                            
                         public :
                            template < ut_dim Dim >
                            auto operator[](const std :: array<t_dim,Dim> index ) const {
                                //std :: cout << "auto operator[](const std :: array<t_dim,Dim> index ) const, Dim=" 
                                //<< Dim << std :: endl;
                                if constexpr (Dim>1){
                                    
                                    std :: array < t_dim , decrease(Dim) > new_index = { 0 } ;
                                    for ( auto i = 0 ; i < Dim-1 ; i ++ ) {
                                        new_index [ i ] = index [ i + 1 ] ;
                                    }

                                    //std :: cout << "ciao" << index[0] << std :: endl;
                                    
                                    //auto a = motherclass::operator[](index[0]);

                                    //return motherclass::operator[](index[0])[new_index];
                                    return (motherclass::at(index[0]))[new_index];
                                
                                }else {
                                    //std :: cout << "ciao2" << index[0] << std :: endl;
                                    return motherclass::at(index[0]);
                                }
                                
                            }
                            
                        
                //***************************************************************
                
                    public :
                        thisclass transpose (const t_dim left = 0, const t_dim right = 1 ) const {
                        
                            if(left==right){return *this;}
                        
                            //std :: cout << "transpose" << std :: endl ;
                            
                            auto e = this->end();
                            auto _e = this->last();
                            for(auto i=0;i<Dimension;++i){
                                _e[i]++;
                            }
                            std :: array < t_dim , Dimension > Lenght = _e;                            
                            Lenght[left] = _e[right];
                            Lenght[right]= _e[left];                            
                        
                            thisclass out (Lenght);           
                            
                            //std :: cout << "out" << std :: endl ;
                            t_dim appo =0;
                            
                            auto b = this->begin();                            
                            auto j = b; //solo per avere il tipo giusto                            
                            //int k = 0 ;
                            
                            //std :: cout << "cycle" << std :: endl ;
                            
                            for(auto i=b;i!=e;i=this->next(i)){
                            
                                //std :: cout << k ++  << std :: endl ;
                                
                                //for(auto kk=0;kk<Dimension;++kk){
                                //    std::cout << i[kk] << ",";
                                //}
                                //std :: cout << std :: endl;
                            
                                j=i;
                                j[left]=i[right];
                                j[right]=i[left];
                                //std::swap(j[left],i[right]);                                
                                out[j]=(*this)[i];
                            
                            }
                            
                            return out;
                        
                        }
                
                //***************************************************************
              
                    public :
                        Index_t begin() const {                            
                            Index_t index; 
                            index.fill(0);
                            return index;
                        }
                        
                    public :
                        Index_t next(const Index_t now ) const {
                        
                            auto size=this->size();
                            if(now[0]>=size){                                
                                return this->end();//now
                            }
                            
                            std :: array < t_dim , Dimension -1 > sub_index;
                            for(auto i=0;i<Dimension-1;++i){
                                sub_index[i]=now[i+1];
                            }
                            
                            auto A = (*this)[now[0]];
                            auto a = A.next(sub_index);
                            //auto a = motherclass::operator[](now[0]).next(sub_index);
                            
                            
                            if (a==sub_index){                                
                                if(now[0]==size-1){
                                    return this->end();//now
                                }else{
                                    Index_t index;
                                    index.fill(0);
                                    index[0]=now[0]+1;
                                    return index;
                                }                                
                            } else {
                                Index_t index;
                                std::copy(a.begin(),a.end(),std::next(index.begin()));
                                index[0]=now[0];
                                return index;
                            }
                        
                        }
                        
                    public :
                        Index_t last() const {
                        
                            Index_t now ;
                            Index_t last ;
                            Index_t index =this->begin(); 
                            do{
                                last=now;
                                now=index;
                                index = this->next(now);  
                            }while(index!=now);     
                            
                            return last;
                        }
                        
                        Index_t end() const {                        
                            Index_t index;
                            index.fill(-1);                 
                            return index;
                        }
                            
                          #undef motherclass 
                          #undef childclass
                          #undef thisclass
                        

                };

            // class specialization 
            template < class type , template <class...> class STL_container > 
                class dynamic_vector < 1 , type , STL_container > 
                    : public  my_dynamic_vector < type , STL_container > {

                    #define motherclass my_dynamic_vector<type,STL_container>
                    
                    public :
                        typedef std :: array < t_dim ,1> Index_t;
                      
                    public :             
                        dynamic_vector ( std :: array < t_dim , 1 > Lenght = { 0 } ) 
                            : my_dynamic_vector < type , STL_container > ( Lenght ) {}
                            
                    public :
                        Index_t begin() const {
                            Index_t index ;
                            index[0]=0;
                            return index;
                        }
                    
                    public :
                        Index_t last() const {
                            Index_t index ;
                            index[0]=this->size();
                            return index;
                        }
                        
                    public :
                        Index_t end() const {
                            Index_t index ;
                            index[0]=-1;
                            return index;
                        }
                            
                    public :
                        Index_t next(const Index_t now ) const {

                            auto size = this->size();
                            Index_t index ;
                            if (now[0]>=size){ //out of range
                                index[0]=-1;
                            }else if(now[0]==size-1){ // the last one
                                index[0]=size-1;
                            }else{
                                index[0]=now[0]+1; //return the next
                            }

                            return index;

                        }
                            
                    
                    public :
                        type & operator[](const t_dim index ){
                            //std :: cout << "type & operator[](const t_dim index ) " << index << std :: endl;   
                            //return motherclass::operator[](index);
                            return motherclass::at(index);
                        }
                        
                    public :
                        type operator[](const t_dim index ) const {
                            //std :: cout << "type operator[](const t_dim index ) const " << index << std :: endl;   
                            //return motherclass::operator[](index);
                            return motherclass::at(index);
                        }
                        
                    public :
                        type & operator[](const Index_t index ){
                            //std :: cout << "type & operator[](const std :: array < t_dim,1> index ) " 
                            //<< index[0] << std :: endl;   
                            //return motherclass::operator[](index[0]);
                            return motherclass::at(index[0]);
                        }
                        
                    public :
                        type operator[](const Index_t index ) const {
                            //std :: cout << "type operator[](const std :: array < t_dim,1> index ) const " 
                            //<< index[0] <<std :: endl;  
                            //return motherclass::operator[](index[0]);
                            return motherclass::at(index[0]);
                        }
                        
                    /*public : 
                            bool set ( const std :: array<t_dim,1> index , const type value ){
                                (*this)[index]=value;
                                return true;
                            }*/
                        
                    #undef motherclass  

                };  
        
        }
                
    }           
}
            
 