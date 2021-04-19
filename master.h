#ifndef LIBRARIES_H_INCLUDED
#define LIBRARIES_H_INCLUDED

#define _USE_MATH_DEFINES

//{--------------------------------------------------Libraries--------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <math.h>
#include <algorithm>
#include <vector>
#include <map>
#include <memory>
#include <Eigen/SparseCore>
#include <boost/dynamic_bitset.hpp>
#include <mgl2/fltk.h>

// delete if unused!
#include <valarray>
//}

//{--------------------------------------------------Functors-------------------------------------------------
class Valmatrix {
public:
 std::valarray<double> data;
 int width, depth, height;

 Valmatrix(int x, int y, int z) : data(x*y*z), width(x), depth(y), height(z) {} //constructor
 Valmatrix(const Valmatrix& m) {data = m.data; width = m.width; depth = m.depth; height = m.height;} //copy constructor
 double& operator()(int x, int y, int z) {return data[x + width*(y + depth*z)];}
 int length() {return width*depth*height;}
 void circshift(char dim, int val, Valmatrix& m) const{
     int index = 0;
     int chunk = 0;
     int shift_mult = 0;
     switch (dim){
        case 'x':{
            chunk = width;
            shift_mult = 1;
            break;}
        case 'y':{
            chunk = width*depth;
            shift_mult = width;
            break;}
        case 'z':{
            chunk = width*depth*height;
            shift_mult = width*depth;
            break;}
     }
     while(index < width*depth*height){
        m.data[std::slice(index, chunk, 1)] = data[std::slice(index, chunk, 1)].cshift(-1*shift_mult*val);
        index += chunk;
     }
 }
};

//}

//{-------------------------------------------------Functions-------------------------------------------------
inline bool is_on_plate(int& x, int& y, int& z, int& max_x, int& max_y, int& max_z, std::array<int,3> (&i)){
 return((x+i[0]<=max_x-1) && (x+i[0]>=0) && (y+i[1]<=max_y-1) && (y+i[1]>=0) && (z+i[2]<=max_z-1) && (z+i[2]>=0));
}
//inline bool is_on_plate(int& i, int& max_i, int& j){
// return((i+j>0) && (i+j<=max_i));
//}

//}

//{--------------------------------------------------Classes--------------------------------------------------
class Plate{
public:
    int dt;
    int dx;
    int diff_cnt;
    std::string diff_type;
    int x;
    int y;
    int z;
    long plate_size;
    int t;
    int agar_height;
    boost::dynamic_bitset<> occupancy_space;
    std::map<std::string, std::vector<std::array<int,3>>> neighbours {
        {"Moore",{{-1,-1,-1},{0,-1,-1},{1,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{-1,1,-1},{0,1,-1},{1,1,-1},
                  {-1,-1,0},{0,-1,0},{1,-1,0},{-1,0,0},{1,0,0},{-1,1,0},{0,1,0},{1,1,0},
                  {-1,-1,1},{0,-1,1},{1,-1,1},{-1,0,1},{0,0,1},{1,0,1},{-1,1,1},{0,1,1},{1,1,1}}},
        {"Neumann",{{0,0,-1},{0,0,1},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0}}}
    };
//    std::map<std::string, std::vector<int>> neighbours {
//        {"Moore",{-1-x-x*y, 0-x-x*y, 1-x-x*y, -1+0-x*y, 0+0-x*y, 1+0-x*y, -1+x-x*y, 0+x-x*y, 1+x-x*y,
//                  -1-x+0, 0-x+0, 1-x+0, -1+0+0, 1+0+0, -1+x+0, 0+x+0, 1+x+0,
//                  -1-x+x*y, 0-x+x*y, 1-x+x*y, -1+0+x*y, 0+0+x*y, 1+0+x*y, -1+x+x*y, 0+x+x*y, 1+x+x*y}},
//        {"Neumann",{0+0-x*y, 0+0+x*y, -1+0+0, 1+0+0, 0-x+0, 0+x+0}}
//    };

    Plate(){}

    Plate(std::ifstream& fin, std::string& line){
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        getline(fin,line,'\n');
        while(!line.empty()){
            std::stringstream s(line);
            s >> std::skipws >> dt;
            s >> std::skipws >> dx;
            s >> std::skipws >> diff_cnt;
            s >> std::skipws >> diff_type;
            s >> std::skipws >> x;
            s >> std::skipws >> y;
            s >> std::skipws >> z;
            s >> std::skipws >> t;
            s >> std::skipws >> agar_height;
            agar_height-=1;
            plate_size = x*y*z;
            getline(fin,line,'\n');
        }
        occupancy_space.resize(x*y*z);
    }

    ~Plate(){}
};

class Nutrient{
public:
    std::string name;
    float init_density;
    float diff_const_agar; //diffusion constant in the agar
    float diff_const_air; //diffusion constant in air
    float diff_const_cell; //diffusion constant near cells, should be vector?
    Eigen::SparseVector<float> density_space;
    Eigen::SparseMatrix<float> laplace_mat;
    std::vector<Eigen::Triplet<float>> laplace_triplet_list;
    int laplace_nutrient_size; //stores the size of laplace_triplet_list after the initialization
                               //laplace_triplet_list's elements will be erased after this point in each update

    Nutrient(){}

    Nutrient(std::ifstream& fin, std::string& line, Plate*& plate){
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        getline(fin,line,'\n');
        while(!line.empty()){
            std::stringstream s(line);
            s >> std::skipws >> name;
            s >> std::skipws >> init_density;
            s >> std::skipws >> diff_const_agar;
            s >> std::skipws >> diff_const_air;
            s >> std::skipws >> diff_const_cell;
            getline(fin,line,'\n');
        }
        laplace_mat.resize(plate->plate_size,plate->plate_size);
        density_space.resize(plate->plate_size);
        init_matrices(plate);
    }

    ~Nutrient(){}

    void init_matrices(Plate* plate){
//        std::vector<Eigen::Triplet<float>> laplace_triplet_list;
        laplace_triplet_list.reserve(7*(plate->plate_size));
        density_space.reserve(plate->x*plate->y*((diff_const_agar>0)*plate->agar_height+(diff_const_air>0)*(plate->z-plate->agar_height))); //big brain move

        for(auto x=0; x<plate->x; x++){
            for(auto y=0; y<plate->y; y++){
                for(auto z=0; z<plate->z; z++){
                    float counter = 0;
                    float diff_const = (z<=plate->agar_height)*diff_const_agar + (z>plate->agar_height)*diff_const_air;
                    float diff_const_up = (z+1<=plate->agar_height)*diff_const_agar + (z+1>plate->agar_height)*diff_const_air;
                    float diff_const_down = (z-1<=plate->agar_height)*diff_const_agar + (z-1>plate->agar_height)*diff_const_air;
                    if(diff_const!=0){
                        for(auto i : plate->neighbours[plate->diff_type]){
                            if(is_on_plate(x,y,z,plate->x,plate->y,plate->z,i) && !((i[2]==1)&&(diff_const_up==0)) && !((i[2]==-1)&&(diff_const_down==0))){
                                laplace_triplet_list.push_back(Eigen::Triplet<float>(x+(plate->x)*(y+(plate->y)*z),x+i[0]+(plate->x)*(y+i[1]+(plate->y)*(z+i[2])),-1*diff_const));
                                counter+=diff_const;
                            }
                        }
                        laplace_triplet_list.push_back(Eigen::Triplet<float>(x+(plate->x)*(y+(plate->y)*z),x+(plate->x)*(y+(plate->y)*z),counter));
                        density_space.insert(x+(plate->x)*(y+(plate->y)*z)) = init_density;
                    }
                }
            }
        }

//        for(auto i=0; i<plate->plate_size; i++){
//            float counter = 0;
//            int z= i/(plate->x*plate->y);
//            float diff_const = (z<=plate->agar_height)*diff_const_agar + (z>plate->agar_height)*diff_const_air;
//            float diff_const_up = (z+1<=plate->agar_height)*diff_const_agar + (z+1>plate->agar_height)*diff_const_air;
//            float diff_const_down = (z-1<=plate->agar_height)*diff_const_agar + (z-1>plate->agar_height)*diff_const_air;
//            if(diff_const!=0){
//                for(auto j : plate->neighbours[plate->diff_type]){
//                    if(is_on_plate(i, plate->plate_size, j) && !((i[2]==1)&&(diff_const_up==0)) && !((i[2]==-1)&&(diff_const_down==0))){
//                        laplace_triplet_list.push_back(Eigen::Triplet<float>(x+(plate->x)*(y+(plate->y)*z),x+i[0]+(plate->x)*(y+i[1]+(plate->y)*(z+i[2])),-1*diff_const));
//                        counter+=diff_const;
//                    }
//                }
//                laplace_triplet_list.push_back(Eigen::Triplet<float>(x+(plate->x)*(y+(plate->y)*z),x+(plate->x)*(y+(plate->y)*z),counter));
//                density_space.insert(x+(plate->x)*(y+(plate->y)*z)) = init_density;
//            }

       laplace_mat.setFromTriplets(laplace_triplet_list.begin(),laplace_triplet_list.end());
       laplace_nutrient_size = laplace_triplet_list.size();
    }
};

class Cells{
public:
    Cells(){}

    ~Cells(){}

    //matrices for all cells, like positions and such

    class Species{
    public:
        std::string name;
        float init_E; //initial energy
        float div_threshold; //division threshold (in terms of energy)
        float div_distance; //division distance
        float death_threshold; //death threshold (in terms of energy)
        float g0_factor; //penalty in g0 state (constant multiplier)
        float g0_threshold; //g0 threshold (in terms of energy)
        float metab_E; //energy usage by metabolism
        float growth_type; //type of colony
        float branch_prob; //branching probability
        float div_dir_dev; //division direction deviation
        float init_cell_num; //initial cell number
        float init_cell_dev; //initial cell deviation
        float init_pos_x;
        float init_pos_y;
        float init_pos_z;
        float food_uptake; //nutrient uptake
        float food_uptake_eff; //nutrient uptake efficiency

        Species(){}

        Species(std::ifstream& fin, std::string& line, Plate*& plate){
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            getline(fin,line,'\n');
            while(!line.empty()){
                std::stringstream s(line);
                s >> std::skipws >> name;
                s >> std::skipws >> init_E;
                s >> std::skipws >> food_uptake;
                s >> std::skipws >> food_uptake_eff;
                s >> std::skipws >> div_threshold;
                s >> std::skipws >> div_distance;
                s >> std::skipws >> death_threshold;
                s >> std::skipws >> metab_E;
                s >> std::skipws >> growth_type;
                s >> std::skipws >> g0_factor;
                s >> std::skipws >> g0_threshold;
                s >> std::skipws >> branch_prob;
                s >> std::skipws >> div_dir_dev;
                div_dir_dev /= M_PI;
                s >> std::skipws >> init_cell_num;
                s >> std::skipws >> init_cell_dev;
                s >> std::skipws >> init_pos_x;
                s >> std::skipws >> init_pos_y;
                init_pos_z = plate->agar_height;
                getline(fin,line,'\n');
            }
        }

        ~Species(){}
    };
};

class mglEigenVec : public mglDataA{
public:
  long nx;
  long ny;
  long nz;
  Eigen::SparseVector<float>* a = nullptr;

  inline mglEigenVec(){}

  inline mglEigenVec(Plate* plate, Eigen::SparseVector<float>* d)
  { nx=plate->x; ny=plate->y; nz=plate->z; a=d; }

  ~mglEigenVec(){}

  inline long GetNx() const { return nx; }
  inline long GetNy() const { return ny; }
  inline long GetNz() const { return nz; }

  inline mreal Maximal() const  { return (a->coeffs()).maxCoeff(); }
  inline mreal Minimal() const  { return /*(a->coeffs()).minCoeff()*/0; }

protected:
  inline mreal v(long i,long j=0,long k=0) const
  { return a->coeff(i+nx*(j+ny*k)); }

  inline mreal vthr(long i) const
  { return a->coeff(i); }

  inline mreal dvx(long i,long j=0,long k=0) const
  { long i0=i+nx*(j+ny*k);
    float res=i>0? (i<nx-1? (a->coeff(i0+1)-a->coeff(i0-1))/2.:a->coeff(i0)-a->coeff(i0-1)) : a->coeff(i0+1)-a->coeff(i0);
    return res;  }

  inline mreal dvy(long i,long j=0,long k=0) const
  { long i0=i+nx*(j+ny*k);
    float res=j>0? (j<ny-1? (a->coeff(i0+nx)-a->coeff(i0-nx))/2.:a->coeff(i0)-a->coeff(i0-nx)) : a->coeff(i0+nx)-a->coeff(i0);
    return res;  }

  inline mreal dvz(long i,long j=0,long k=0) const
  { long i0=i+nx*(j+ny*k), n=nx*ny;
    float res=k>0? (k<nz-1? (a->coeff(i0+n)-a->coeff(i0-n))/2.:a->coeff(i0)-a->coeff(i0-n)) : a->coeff(i0+n)-a->coeff(i0);
    return res;  }

  inline mreal valueD(mreal x,mreal y=0,mreal z=0,mreal *dx=0,mreal *dy=0,mreal *dz=0) const
  { return 0;  }

  inline mreal value(mreal x,mreal y=0,mreal z=0) const
  { return 0;  }

};

class yeastDraw : public mglDraw{
public:
    int is_draw;
    int steps;
    int resolution;
    int* x = nullptr;
    int* y = nullptr;
    int* z = nullptr;
    std::vector<mglEigenVec> nutrients;

    yeastDraw(){}

    yeastDraw(std::ifstream& fin, std::string& line, Plate*& p, std::vector<Nutrient>& _nutrients){
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        getline(fin,line,'\n');
        while(!line.empty()){
            std::stringstream s(line);
            s >> std::skipws >> is_draw;
            s >> std::skipws >> steps;
            s >> std::skipws >> resolution;
            getline(fin,line,'\n');
        }

        nutrients.reserve(_nutrients.size()); //Hell yeah!
        for(auto& n:_nutrients){
            nutrients.emplace_back(p, &n.density_space);
        }

        x = &(p->x);
        y = &(p->y);
        z = &(p->z);
    }

    ~yeastDraw(){}

    int Draw(mglGraph* gr){
     Agar(gr);
     return 0;
    }

    void Agar(mglGraph* gr){
     gr->Title("MathGL Demo");
     gr->Rotate(60,40);
     gr->Aspect(*x,*y,*z);
     gr->Alpha(true);
     gr->Box();
     gr->Cloud(nutrients[0],"BBBBBB{xFFFFFF00}BbcyrR");//BBBBBB{xFFFFFF00}BbcyrR
    }
};

//}


#endif // LIBRARIES_H_INCLUDED
