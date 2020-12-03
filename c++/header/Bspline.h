//%%file header/Bspline.h
#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <string>
#include <iterator>

#include "eliastocco_namespace_Bspline.h"
#include "knot_vector.h"
#include "dynamic_vector.h"
#include "Bspline_operator.h"


namespace eliastocco {

    namespace polynomial {

        namespace Bspline {        
 
            template < class type , ut_dim dimension = 1 , ut_dim codimension = 1 >
                class Bspline {

                    #define thisclass Bspline<type,dimension,codimension>    

                    public :
                        typedef std :: array<type,codimension> Type_out;
                        typedef std :: array<type,dimension>   Type_in; 
                        typedef std :: array<Type_out,dimension>   jacobian_t;
                    
                        template < class Type , ut_dim Dim >
                            using square_matrix_t = std :: array<std :: array < Type, Dim > , Dim>;

                        typedef std :: array<t_dim,dimension>  Index_t;

                        template < ut_dim Dim >
                            using template_knots_t = std :: array < knot_vector<type> , Dim>;
                        typedef template_knots_t<dimension> knots_t;
                        
                        typedef std :: array<type,codimension> single_cp_t;
                        template < ut_dim Dim >
                            using template_cp_t = container::dynamic::dynamic_vector < Dim , single_cp_t >;                    
                        
                        typedef template_cp_t<dimension> cp_t;

//********************************************************************************************************************
                    
                    private :
                        thisclass(){};
                        
                    public :
                        thisclass ( const knots_t knots , const cp_t cp )
                            :thisclass(knots){
                                _cp = cp ;
                            };
                        
                    //non sto controllando la correttezza degli indici dei control points
                    //non sto controllande che ogni control point abbia la lunghezza giusta
                    //per comodità ho deciso di optare per dei control points std :: valarray
                    public :
                        thisclass ( const knots_t knots )
                            : _kv(knots){

                                //std :: cout << "constructor" << std :: endl;
                                Index_t init = Index_t();

                                for(auto i=0; i<dimension; ++i){

                                    auto p = _kv[i].p();
                                    auto n = _kv[i].n();
                                    auto k = _kv[i].knots();
                                    auto v = k.size();

                                    init[i]=n;

                                    assert ( std::is_sorted(k.begin(),k.end()) %% "knot vector not sorted" );
                                    assert ( p > 0 %% "wrong polyniomial degree" );
                                    assert ( n > 1 %% "wrong basis cardinality" );
                                    assert ( p+n+1 == v %% "wrong knot vector size" );

                                }

                                //std :: cout << "_cp" << std :: endl;
                                _cp = cp_t (init);
                            
                            };
                            
//********************************************************************************************************************

                    public :
                        Type_out get ( const Index_t index , const bool check = true ) const { 

                            try { 
                                return _cp[index];
                            } catch ( std::out_of_range e ){
                                if ( check ){
                                    std :: cout << "error : Bspline.get" << std :: endl;
                                    std :: cout << e.what()<< std :: endl;  
                                }
                                return Type_out();
                            } 
                            
                        };                    

//********************************************************************************************************************

                    public : // piccola modifica
                        bool set (const Index_t index, const Type_out value , const bool check = true ) {
                            //std :: cout << "set" << std :: endl;
                            try {
                                _cp[index]=value; 
                            } catch ( std::out_of_range e ){  
                                if ( check ){
                                    std :: cout << "error : Bspline.set" << std :: endl;
                                    std :: cout << e.what()<< std :: endl;  
                                    return false;
                                }
                            } 

                            return true;

                        };

//********************************************************************************************************************

                     protected :
                        static int find_k( const type x , const knot_vector<type> t ) {
                            
                            //std :: cout << "indide find_k" << std :: endl;
                            //std :: cout << "size =" << t.size() << std :: endl;
                            int output = -1;
                            for(auto i = 0 ; i < t.size()-1; ++i){
                                //std :: cout << "i: =" << i << std :: endl;
                                if ( t[i] <= x and x < t[i+1] ){ output = i ;}   
                            }

                            return output;

                        };

//********************************************************************************************************************
                    
                    public :
                        Type_out deBoor ( const type x , const bool derivative = false ) const {
                            
                            // https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
                            
                            //std :: cout << "x:" << x << std :: endl;

                            static_assert ( dimension == 1 , "dimension for deBoor algorithm is greater than 1");

                            //
                            auto v = knots_vector();
                            if ( v.size() > 1 ){
                                std :: cout << "deBoor error : v.size()>1" << std :: endl ;
                                throw std::exception ();
                            }  
                            
                            //
                            knot_vector<type> t = v[0];
                            int p = t.p();
                            template_cp_t<1> c = control_points();                                                     
                            
                            //
                            int k = find_k(x,t); 
                            if ( k < 0 ){        
                                std :: cout << "deBoor error : k<0" << std :: endl ;
                                return Type_out();
                            }

                            //
                            std :: array<t_dim,1> init ;
                            init[0]=p+1;
                            template_cp_t<1> d   (init);
                            //template_cp_t<1> q   (init);
                            
                        //****************************
                            if (derivative){
                                
                                //template_cp_t<1> q   (init);
                            
                                //
                                //std :: cout << "derivative" << std :: endl ;
                                auto jj=0;
                                for(auto j=0; j!= p ; ++j){
                                    //std :: cout << "j : " << j << std :: endl ;
                                    //std :: cout << "left" << std :: endl ;
                                    auto a = j+k-p+1;
                                    auto left  = a>=0 and a<c.size() ? c[a] : Type_out();
                                    //std :: cout << "right" << std :: endl ;
                                    auto b = j+k-p;
                                    auto right = b>= 0 and b<c.size() ? c[b] : Type_out();
                                    //std :: cout << "factor" << std :: endl ;
                                    auto factor = p/(t[j+k+1] - t[j+k-p+1]);
                                    //std :: cout << "q" << std :: endl ;
                                    d[j]=factor*(left+(-1.0)*right);
                                    
                                }
                                
                                 //
                                //std :: cout << "deBoor" << std :: endl ;
                                type alpha  =type();
                                //type dalpha =type();
                                for(auto r = 1 ; r !=p+1; ++r){
                                    //std :: cout << "r : " << r << std :: endl ;
                                    for(auto j = p-1; j!=r-1; --j){
                                        //std :: cout << "r : " << r << std :: endl ;
                                        auto right = j+1+k-r;
                                        auto left = j+k-(p-1);
                                        alpha = (x - t[left]) / (t[right] - t[left]);
                                        d[j] = (1.0 - alpha) * d[j-1] + alpha * d[j];

                                    }
                                }                                
                                
                                Type_out output;
                                auto out = d[p-1];                                
                                std::copy(std::begin(out),std::end(out),std::begin(output));                                
                                return output;                                
                                
                            } else {
  
                    //****************************       

                                //
                                auto jj=0;
                                for(auto j=0; j!= p+1 ; ++j){
                                    if(j+k-p>=0 and j+k-p < c.size() ){
                                        d[j]=c[j+k-p];
                                    }else {
                                        d[j] = Type_out();
                                    }
                                }

                                //
                                std :: cout << "k : " << k << std :: endl ; 
                                std :: cout << "c.size() : " << c.size() << std :: endl ; 
                                std :: cout << "t.size() : " << t.size() << std :: endl ; 
                                type alpha = type();
                                for(auto r = 1 ; r !=p+1; ++r){
                                    std :: cout << "r : " << r << std :: endl ; 
                                    for(auto j = p; j!=r-1; --j){
                                        std :: cout << "j : " << j << std :: endl ; 
                                        std :: cout << "j+k-p : " << j+k-p << std :: endl ; 
                                        std :: cout << "j+1+k-r : " << j+1+k-r << std :: endl ; 
                                        std :: cout << "t(" << j+1+k-r <<") :" << t[j+1+k-r] << std :: endl ; 
                                        alpha  = (x - t[j+k-p]) / (t[j+1+k-r] - t[j+k-p]);
                                        std :: cout << "alpha(" << j <<") :" << alpha << std :: endl ; 
                                        d[j] = (1.0 - alpha) * d[j-1] + alpha * d[j];                                       
                                    }
                                }

                                //
                                Type_out output;
                                single_cp_t out = d[p];                   
                                std::copy(std::begin(out),std::end(out),std::begin(output));
                                return output;
                            
                            }


                        };
                    
//********************************************************************************************************************

                    protected :
                       Bspline<type,1,codimension> iterative_deBoor (const Type_in x) const {
                           
                            // https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/surface/bspline-de-boor.html

                                //
                                std :: array<type,dimension-1> x_surf;
                                std :: copy(std::next(x.begin()),x.end(),x_surf.begin());

                                //
                                template_knots_t<dimension-1> curve_knot = {} ;
                                std :: copy ( std :: next (_kv.begin()) , _kv.end() , curve_knot.begin() );

                                //
                                t_dim m  = _kv[0].n();
                                    
                                //
                                template_cp_t<1> surface_cp = {{m}};  
                                Bspline<type,dimension-1,codimension> curve(curve_knot);

                                //                                    
                                for(auto i=0; i<m ; ++i){                                       
                                    curve.set_control_points(_cp[i]);                                        
                                    surface_cp[i]=curve(x_surf);
                                }
                                    
                                //
                                template_knots_t<1> sub_kv = {_kv[0]};
                                    
                                //                           
                                return Bspline<type,1,codimension>(sub_kv,surface_cp);
                           
                       }
                    
//********************************************************************************************************************

                    //valutazione del polinomio in x
                    public :
                        Type_out operator () ( const Type_in x ) const {
                                if constexpr ( dimension == 1 ){                                    
                                    return deBoor(x[0]);                                       
                                } else {
                                    auto curve_final = iterative_deBoor(x); 
                                    return curve_final.deBoor(x[0]);
                                }
                        };

                    //vectorized polinomial evaluation
                    public :
                        template < template < class ... > class vect = std :: vector >
                            vect <Type_out> operator () ( const vect < Type_in > x ) const { 
                                vect <Type_out> y (x.size());
                                std :: transform( std :: begin(x), std :: end(x), std :: begin(y),
                                    [&]( const auto a ) { return (*this)(a);});   
                                return y;

                            };
                    
//********************************************************************************************************************
                    
                    public :
                        thisclass transpose (const t_dim left , const t_dim right)const{
                            
                            if (left==right){return *this;}
                            
                            knots_t new_kv = this->knots_vector();
                            std::swap(new_kv[left],new_kv[right]); 
                            
                            auto new_cp = this->control_points();
                            new_cp = new_cp.transpose(left,right);
                            
                            return thisclass (new_kv,new_cp);
                            
                        }
                 
                    public :
                        jacobian_t jacobian ( const Type_in x ) const {
                            
                            //curve costanti a tratti: lo jacobian è nullo
                            if(_kv[0].p()<=0){ return jacobian_t();}
                            
                            if constexpr (dimension==1){   
                                
                                //questo funziona
                                auto out = deBoor(x[0],true);                                
                                jacobian_t output;                                
                                for(auto i=0;i<codimension;++i){
                                    output[0][i]=out[i];
                                }                                 
                                return output;
                                    
                            } else {
                                
                                jacobian_t output;                                
                                thisclass newBspline;
                                
                                for(auto i=0;i<dimension;++i){                                    
                                   
                                    auto newBspline = this->transpose(0,i);
                                    Bspline<type,1,codimension> curve_final = newBspline.iterative_deBoor(x);
                                    
                                    std :: array<type,1> X = {x[0]};                                    
                                    std::array<std::array<type,codimension>,1> a = curve_final.jacobian(X);                                    
                                    std::array<double,codimension> b = a[0];                                    
                                    output[i]=b;
                                }
                                
                                return output;            
                                
                            }
                            
                        }
                    
                    public :
                        template < template < class ... > class vect = std :: vector >
                            vect <jacobian_t> jacobian ( const vect < Type_in > x ) const { 
                                vect <jacobian_t> y (x.size());
                                std :: transform( std :: begin(x), std :: end(x), std :: begin(y),
                                    [&]( const auto a ) { return jacobian(a);});   
                                return y;

                            };
                    
//********************************************************************************************************************                    
                    
                    public :
                        template < ut_dim Dim >
                            type determinant ( const square_matrix_t<type,Dim> mat ) const {

                                std :: cout << "determinant " << Dim << std :: endl;
                                if constexpr ( Dim == 1 ){
                                    return mat[0][0];
                                } else if ( Dim == 2 ){
                                    return (mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]);
                                } else {
                                    
                                    square_matrix_t<type,Dim-1> submat;
                                    type d = type ();   
                                    for (auto c = 0; c < Dim; ++c){                                        
                                        int subi = 0; //submatrix's i value                                        
                                        for (auto i = 1; i < Dim; ++i){
                                            int subj = 0;
                                            for (auto j = 0; j < Dim; ++j){
                                                if (j == c){
                                                    continue;
                                                }

                                                submat[subi][subj] = mat[i][j];
                                                subj++;                
                                            }          
                                            subi++;
                                        }
                                        d += std::pow(-1, c) * mat[0][c] * determinant(submat);
                                    }      
                                    
                                   return d ; 
                                }

                            };
                    
                    public :
                        type jacobian_determinant ( const Type_in x ) const { 
                            
                            static_assert(codimension>=dimension,"error:codimension < dimension");
                            
                            auto jac = this->jacobian(x);
                            
                            if constexpr (dimension==codimension){
                                
                                std :: cout << "jacobian_determinant " << "dimension==codimension" << std :: endl;
                                
                                auto d = determinant(jac);
                                return std :: abs (d);
                                
                            } else if (dimension == 1){
                                
                                std :: cout << "jacobian_determinant " << "dimension == 1" << std :: endl;
                                              
                                type det = type();                                
                                for(auto i=0;i<codimension;++i){
                                   det += std :: pow(jac[i][0],2.0);
                                }
                                return std :: sqrt(det);
                                
                            } else {           
                                
                                std :: cout << "jacobian_determinant " << "else" << std :: endl;
                                return type();
                            }
                            
                        };
                    
                    public :
                        template < template < class ... > class vect = std :: vector >
                            vect <type> jacobian_determinant ( const vect < Type_in > x ) const { 
                                vect <type> y (x.size());
                                std :: transform( std :: begin(x), std :: end(x), std :: begin(y),
                                    [&]( const auto a ) { return jacobian_determinant(a);});   
                                return y;

                            };
                    
//********************************************************************************************************************

                    public :
                        inline auto knots(const t_dim index) const { return _kv[index].knots(); }                    
                    public :
                        inline knots_t knots_vector() const { return _kv; }                    
                    public :
                        inline auto p() const { return _kv.p(); }                    
                     public :
                        inline auto n() const { return _kv.n(); }
                    public :
                        inline cp_t control_points() const { return _cp; }
                    
//********************************************************************************************************************
                    
                    public :
                        auto clear (){                            
                                Index_t init = Index_t();
                                for(auto i=0; i<dimension; ++i){                               
                                    init[i]=_kv[i].n();      
                                }
                                _cp = cp_t (init);
                            
                            return init;
                            
                        };
                        
//********************************************************************************************************************

                     protected :
                        inline void set_knots(const knots_t user_kv ) { _kv = user_kv ; }
                        
                     public :
                        inline void set_control_points(const cp_t user_cp ) { _cp = user_cp ; }


//********************************************************************************************************************

                    protected :
                        cp_t _cp; //control points

                    protected :
                        knots_t _kv; //knot vectors
                        

                    #undef thisclass
                };

            }

        }

}