#ifndef LIBRARIES_H_INCLUDED
#define LIBRARIES_H_INCLUDED

//{-----------------------------------Libraries-----------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <math.h>
#include <algorithm>
#include <vector>
#include <Eigen/SparseCore>

// delete if unused!
#include <map>
#include <valarray>
//}

//{-----------------------------------Functors-----------------------------------
class Valmatrix {
public:
 std::valarray<double> data;
 int width, depth, height;

 Valmatrix(int x, int y, int z) : data(x*y*z), width(x), depth(y), height(z) {} //constructor
 Valmatrix(const Valmatrix& m) {data = m.data; width = m.width; depth = m.depth; height = m.height;} //copy constructor
 double& operator()(int x, int y, int z) {return data[x + width*(y + depth*z)];}
 int length() {return width*depth*height;}
 void circshift(char dim, int val, Valmatrix& m) const{
     size_t index = 0;
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

//{-----------------------------------Classes------------------------------------
class Draw{
public:
    int is_draw;
    int steps;
    int resolution;

    Draw(){}

    Draw(std::ifstream& fin, std::string& line){
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        getline(fin,line,'\n');
        while(!line.empty()){
            std::stringstream s(line);
            s >> std::skipws >> is_draw;
            s >> std::skipws >> steps;
            s >> std::skipws >> resolution;
            getline(fin,line,'\n');
        }
    }

    ~Draw(){}
};

class Plate{
public:
    int dt;
    int dx;
    int diff_cnt;
    int max_x;
    int max_y;
    int max_z;
    int max_t;
    int agar_height;

    Plate(){}

    Plate(std::ifstream& fin, std::string& line){
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        getline(fin,line,'\n');
        while(!line.empty()){
            std::stringstream s(line);
            s >> std::skipws >> dt;
            s >> std::skipws >> dx;
            s >> std::skipws >> diff_cnt;
            s >> std::skipws >> max_x;
            s >> std::skipws >> max_y;
            s >> std::skipws >> max_z;
            s >> std::skipws >> max_t;
            s >> std::skipws >> agar_height;
            getline(fin,line,'\n');
        }
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

    Nutrient(){}

    Nutrient(std::ifstream& fin, std::string& line){
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
    }

    ~Nutrient(){}
    //diffusion matrices later
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

        Species(std::ifstream& fin, std::string& line, Plate* plate){
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
//}

#endif // LIBRARIES_H_INCLUDED
