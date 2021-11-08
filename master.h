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
#include <array>
#include <unordered_map>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <random>
#include <Eigen/SparseCore>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <mgl2/fltk.h>
#include <omp.h>

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

inline float inverse_square_root(float n){ //QUAKE
    float y = n;
    std::uint32_t i;
    std::memcpy(&i, &y, sizeof(float));
    i  = 0x5f3759df - ( i >> 1 );
    std::memcpy(&y, &i, sizeof(float));
    return y * ( 1.5F - ( (n * 0.5F )* y * y ) );
}

//}

//{--------------------------------------------------Structs--------------------------------------------------
struct GridCell{
    int x, y, z;

    GridCell(int x, int y, int z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    ~GridCell(){}

    bool operator==(const GridCell& other) const{
        return x == other.x && y == other.y && z == other.z;
    }
};

struct hash_fn{
    std::size_t operator()(const GridCell& other) const{
        std::size_t seed = 0;
        boost::hash_combine(seed, other.x);
        boost::hash_combine(seed, other.y);
        boost::hash_combine(seed, other.z);
        return seed;
    }
};
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
        occupancy_space.resize(plate_size);
    }

    ~Plate(){}

    void set_plate_dimensions(int gridcell_size){
        //update plate sizes to ensure the grids will fit
        x = ceil((float)x/(float)gridcell_size)*gridcell_size;
        y = ceil((float)y/(float)gridcell_size)*gridcell_size;
        z = ceil((float)z/(float)gridcell_size)*gridcell_size;
        plate_size = x*y*z;
    }
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

        laplace_mat.reserve(plate->plate_size*plate->plate_size);
        laplace_mat.resize(plate->plate_size,plate->plate_size);
        density_space.reserve(plate->plate_size);
        density_space.resize(plate->plate_size);
    }

    ~Nutrient(){}

    void set_matrices(Plate*& plate){
        laplace_mat.reserve(plate->plate_size*plate->plate_size);
        laplace_mat.resize(plate->plate_size,plate->plate_size);
        density_space.reserve(plate->plate_size);
        density_space.resize(plate->plate_size);
        init_matrices(plate);
    }

    void init_matrices(Plate*& plate){
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
    class Species{
    public:
        ///initial properties of cells
        std::string name;
        std::string color;
//        float init_E; //initial energy
        std::vector<float> init_E; //initial energy
//        float div_threshold; //division threshold (in terms of energy)
        std::vector<float> div_threshold; //division threshold (in terms of energy)
        float cell_radius; //size of the spherical agent
        float sqrd_radius; //squared cell radius for efficient distance comparisons
//        float death_threshold; //death threshold (in terms of energy)
        std::vector<float> death_threshold; //death threshold (in terms of energy)
        float g0_factor; //penalty in g0 state (constant multiplier)
//        float g0_threshold; //g0 threshold (in terms of energy)
        std::vector<float> g0_threshold; //g0 threshold (in terms of energy)
//        float metab_E; //energy usage by
        std::vector<float> metab_E; //energy usage by metabolism
        float growth_type; //type of colony
        float branch_prob; //branching probability
        float div_dir_dev; //division direction deviation
        int init_cell_num; //initial cell number
        float init_cell_dev; //initial cell deviation
        float init_pos_x;
        float init_pos_y;
        float init_pos_z;
//        float nutrient_uptake; //nutrient uptake
        std::vector<float> nutrient_uptake; //nutrient uptake
//        float uptake_eff; //nutrient uptake efficiency
        std::vector<float> uptake_eff; //nutrient uptake efficiency

        Species(){}

        Species(std::ifstream& fin, std::string& line, Plate*& plate, size_t nutrient_size){
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            getline(fin,line,'\n');
            while(!line.empty()){
                float value;
                std::stringstream s(line);
                s >> std::skipws >> name;
                s >> std::skipws >> color;
                for(size_t i=0;i<nutrient_size;++i){
                    s >> std::skipws >> value;
                    init_E.push_back(value);
                }
//                s >> std::skipws >> init_E;
                for(size_t i=0;i<nutrient_size;++i){
                    s >> std::skipws >> value;
                    nutrient_uptake.push_back(value);
                }
//                s >> std::skipws >> nutrient_uptake;
                for(size_t i=0;i<nutrient_size;++i){
                    s >> std::skipws >> value;
                    uptake_eff.push_back(value);
                }
//                s >> std::skipws >> uptake_eff;
                for(size_t i=0;i<nutrient_size;++i){
                    s >> std::skipws >> value;
                    div_threshold.push_back(value);
                }
//                s >> std::skipws >> div_threshold;
                s >> std::skipws >> cell_radius;
                sqrd_radius = pow(cell_radius,2);
                for(size_t i=0;i<nutrient_size;++i){
                    s >> std::skipws >> value;
                    death_threshold.push_back(value);
                }
//                s >> std::skipws >> death_threshold;
                for(size_t i=0;i<nutrient_size;++i){
                    s >> std::skipws >> value;
                    metab_E.push_back(value);
                }
//                s >> std::skipws >> metab_E;
                s >> std::skipws >> growth_type;
                s >> std::skipws >> g0_factor;
                for(size_t i=0;i<nutrient_size;++i){
                    s >> std::skipws >> value;
                    g0_threshold.push_back(value);
                }
//                s >> std::skipws >> g0_threshold;
                s >> std::skipws >> branch_prob;
                s >> std::skipws >> div_dir_dev;
                div_dir_dev /= M_PI;
                s >> std::skipws >> init_cell_num;
                s >> std::skipws >> init_cell_dev;
                s >> std::skipws >> init_pos_x;
                s >> std::skipws >> init_pos_y;
                init_pos_z = plate->agar_height+cell_radius;
                getline(fin,line,'\n');
            }
        }

        ~Species(){}
    };

    class Agent{
    public:
        Cells::Species *species; //pointer to species of the agent to acces it's attributes
//        float energy; //energy
        std::vector<float> energy; //energy
//        float nutrient_uptake;
        std::vector<float> nutrient_uptake;
        std::array<float,3> pos; //position
        std::vector<int> occupied_uniform_grids; //GridCells which the agent occupies
        std::vector<int> occupied_density_space; // density space matrix cells which the agent occupies
        bool state;
        bool divide;

        Agent(){}

        Agent(Plate *&plate, Cells &cells, Cells::Species *species, std::array<float,3> pos){
            this->species = species;
            this->energy = species->init_E;
            decide_state();
            can_divide(cells);
            this->pos = pos;

            occupied_uniform_grids = occupied_grid_cells({0,0,plate->agar_height},cells.conversion_factor,cells.gridcell_size,cells.gridmap_x,cells.gridmap_y);
            occupied_density_space = occupied_grid_cells({0,0,0},1,1,plate->x,plate->y);
            for(auto &n:species->nutrient_uptake){
                nutrient_uptake.push_back(static_cast<float>(n/occupied_density_space.size()));
            }
//            nutrient_uptake = static_cast<float>(species->nutrient_uptake/occupied_density_space.size());
        }

        //shallow copy constructor
        Agent(const Cells::Agent &old_agent){
            energy = old_agent.energy;
            state = old_agent.state;
            divide = old_agent.divide;
            pos = old_agent.pos;
            nutrient_uptake = old_agent.nutrient_uptake;
            occupied_uniform_grids = old_agent.occupied_uniform_grids;
            occupied_density_space = old_agent.occupied_density_space;
//            species = new Cells::Species;
//            *species = *old_agent.species;
            species = old_agent.species;
        }

        //move constructor
        Agent(Cells::Agent &&old_agent){
            energy = old_agent.energy;
            state = old_agent.state;
            divide = old_agent.divide;
            pos = old_agent.pos;
            nutrient_uptake = old_agent.nutrient_uptake;
            occupied_uniform_grids = old_agent.occupied_uniform_grids;
            occupied_density_space = old_agent.occupied_density_space;
            species = old_agent.species;
            old_agent.species = nullptr;
        }

        //copy assignment
//        Cells::Agent& operator=(const Cells::Agent &old_agent){
//            std::cout<<"shit's been called yo"<<std::endl;
//            energy = old_agent.energy;
//            divide = old_agent.divide;
//            state = old_agent.state;
//            pos = old_agent.pos;
//            nutrient_uptake = old_agent.nutrient_uptake;
//            occupied_uniform_grids = old_agent.occupied_uniform_grids;
//            occupied_density_space = old_agent.occupied_density_space;
//            species = old_agent.species;
//            return *this;
//        }

        //~Agent(){}
        ~Agent(){species = nullptr;}

        void decide_state(){ //decide whether an agent should be active or in G0 state
            state = (energy > species->g0_threshold);
        }

        void can_divide(Cells &cells){
            //check energy level first
            divide = (energy > species->div_threshold);
            if(!divide){return;}

            //check 'density' of grid cells occupied by agent
            //if all of them are 'full', then divide is 0, else 1
            divide = false;
            for(auto &i:occupied_uniform_grids){
                auto range = cells.agent_gridmap.equal_range(i);
                auto density = std::distance(range.first,range.second);
                divide = divide||(density < 24);
            }
        }

        std::vector<int> occupied_grid_cells(std::array<int,3> e, float c, int d, int length, int width){

            std::vector<int> grids;
            float r = species->cell_radius;

            //bounding box
            std::array<std::array<float,3>,2> bounding_box = {{{pos[0]-r,pos[1]-r,pos[2]-r},{pos[0]+r,pos[1]+r,pos[2]+r}}};
            std::array<std::array<int,3>,2> grid_corners;

            //get the grid cells in which box corners lie
            std::transform(bounding_box.begin(),bounding_box.end(),grid_corners.begin(),[&c, &e](std::array<float,3> &n)
                          {
                              std::array<int,3> m;
                              for(auto i=0;i<3;++i){
                                    m[i] = static_cast<int>(c*(n[i] - e[i]));
                              }
                              return m;
                          });

            //check distance between grid cells and agent center
            for(int z=grid_corners[0][2];z<=grid_corners[1][2];++z){
                for(int y=grid_corners[0][1];y<=grid_corners[1][1];++y){
                    int x_min = grid_corners[0][0];
                    int x_max = grid_corners[1][0];
//Expensive trimming of bounding box NEEDS FIXING!!! grid_corners ARE ALREADY TRANSFORMED, BUT TRIMMING ASSUMES ABSOLUTE DISTANCES!!!
//                    if(1.415*r>d){
//                        std::cout<<"find_x_min started"<<std::endl;
//                        find_x_min(x_min, y, z, x_max, d);
//                        //x_min was not found -> no grid cells belong to agent along this y,z -> skip
//                        if(x_min == x_max){
//                            continue;
//                        }
//
//                        std::cout<<"find_x_max started"<<std::endl;
//                        find_x_max(x_min, y, z, x_max, d);
//                    }
                    //save (each grid from x_min to x_max
                    for(auto x=x_min;x<=x_max;++x){
                        grids.push_back(x+length*(y+width*z));
                    }
                }
            }
            return grids;
        }

        inline void find_x_min(int &x_min, int &y, int &z, int &x_max, int &d){
            while(x_min != x_max){
                std::array<double,3> midpoint = {x_min+d*0.5,y+d*0.5,z+d*0.5};
                std::array<double,3> shifted_center = {pos[0]-midpoint[0],pos[1]-midpoint[1],pos[2]-midpoint[2]};
                std::array<double,3> shifted_corner = {x_min+d-midpoint[0],y+d-midpoint[1],z+d-midpoint[2]};
                float sqrd_dist = pow(std::max(0.0,std::abs(shifted_center[0])-shifted_corner[0]),2) +
                                  pow(std::max(0.0,std::abs(shifted_center[1])-shifted_corner[1]),2) +
                                  pow(std::max(0.0,std::abs(shifted_center[2])-shifted_corner[2]),2);

                if(sqrd_dist<species->sqrd_radius){
                    return;
                }

                x_min += d;
            }
        }

        inline void find_x_max(int &x_min, int &y, int &z, int &x_max, int &d){
            while(x_min != x_max){
                std::array<double,3> midpoint = {x_max-d*0.5,y+d*0.5,z+d*0.5};
                std::array<double,3> shifted_center = {pos[0]-midpoint[0],pos[1]-midpoint[1],pos[2]-midpoint[2]};
                std::array<double,3> shifted_corner = {x_max-midpoint[0],y+d-midpoint[1],z+d-midpoint[2]};
                float sqrd_dist = pow(std::max(0.0,std::abs(shifted_center[0])-shifted_corner[0]),2) +
                                  pow(std::max(0.0,std::abs(shifted_center[1])-shifted_corner[1]),2) +
                                  pow(std::max(0.0,std::abs(shifted_center[2])-shifted_corner[2]),2);

                if(sqrd_dist<species->sqrd_radius){
                    return;
                }

                x_max -= d;
            }
        }
    };

    std::list<Cells::Agent> agent_list;
    std::unordered_multimap<int, Cells::Agent*> agent_gridmap;
    std::map<Cells::Species*,std::map<bool,std::array<std::vector<double>,3>>> agent_positions; //will be linked to mglData objects to plot agents
    int gridmap_size;
    int gridmap_x, gridmap_y, gridmap_z;
    int gridcell_size; //size of a grid cell x, y and z-wise
    float conversion_factor;


    Cells(){}

    Cells(Plate*& plate, std::vector<Cells::Species>& species){
        size_of_gridcell(plate, species);//determine gridcell_size based on cell radii
        conversion_factor = 1/(float)gridcell_size;
        num_of_gridcells(plate); //determine gridmap_size and dimensions based on plate and size of one grid cell
        agent_gridmap.reserve(gridmap_size);//reserve space for agent_gridmap
        setup_position_container(species);
    }

    ~Cells(){}

    void setup_position_container(std::vector<Cells::Species>& species){
        for(auto &i:species){
            std::array<std::vector<double>,3> positions;
            agent_positions[&i][0] = positions;
            agent_positions[&i][1] = positions;
        }
    }

    void copy_positions(std::list<Cells::Agent>::iterator agent){
        for(size_t i=0;i<agent->pos.size();++i){
            agent_positions[agent->species][agent->state][i].push_back(agent->pos[i]);
        }
    }

    void clear_positions(){
        for(auto &species_map:agent_positions){
            for(auto &positions_map:species_map.second){
                for(auto &position_vector:positions_map.second){
                    position_vector.clear();
                    //preallocate memory
                    position_vector.reserve(agent_list.size()*2);
                }
            }
        }
    }

    void size_of_gridcell(Plate*& plate, std::vector<Cells::Species>& species){
        //calculate gridcell size based on largest agent in the model
        auto n = std::max_element(species.begin(), species.end(),
                                  [](const Cells::Species& a, const Cells::Species& b)
                                  {return a.cell_radius < b.cell_radius;});

        gridcell_size = ceil(2*2*n->cell_radius);
    }

    void num_of_gridcells(Plate*& plate){
        gridmap_x = ceil(plate->x*conversion_factor);
        gridmap_y = ceil(plate->y*conversion_factor);
        gridmap_z = ceil((plate->z-plate->agar_height-1)*conversion_factor);
        gridmap_size = gridmap_x*gridmap_y*gridmap_z;
    }

    void add_agent(Plate*& plate, Cells::Agent &yeast){
        agent_list.push_front(yeast); //list
        for(auto &i:yeast.occupied_uniform_grids){ //hashmap
            agent_gridmap.emplace(i,&(agent_list.front()));
        }
        for(auto &i:yeast.occupied_density_space){ //occupied cells in density matrices
            plate->occupancy_space.set(i);
        }
    }

    void delete_agent(Plate*& plate, std::list<Cells::Agent>::iterator &agent){
        for(auto j:agent->occupied_uniform_grids){ //go through every uniform grid it occupies
            auto k = agent_gridmap.equal_range(j); //find the grids in the agent_map
            for(auto l=k.first;l!=k.second;++l){ //go through all occupied grids and delete agent to be killed
                if(l->second == &(*agent)){
                    //delete from hashmap
                    agent_gridmap.erase(l);
                    break;
                }
            }
        }
        for(auto &i:agent->occupied_density_space){ //occupied cells in density matrices
            plate->occupancy_space.reset(i);
        }
        //delete from list
        agent = agent_list.erase(agent);
        --agent;
    }

    inline void init_positions(Plate*& plate, Cells::Species *s){
        std::vector<std::array<float,2>> positions;
        float x = s->init_pos_x /*+ s->init_cell_dev * (float) rand()/RAND_MAX*/;
        float y = s->init_pos_y /*+ s->init_cell_dev * (float) rand()/RAND_MAX*/;
        positions.push_back({x,y});
        Cells::Agent yeast(plate, *this, s, {x,y,s->init_pos_z});
        add_agent(plate, yeast);
        if(s->init_cell_num==1){return;}

        int counter = 1;
        float d = s->cell_radius*2;
        while(counter < s->init_cell_num){
            //get a trial point
            float x = s->init_pos_x + s->init_cell_dev * (float) rand()/RAND_MAX;
            float y = s->init_pos_y + s->init_cell_dev * (float) rand()/RAND_MAX;

            //see how close it is to existing points
            std::vector<float> dist;
            for(auto i:positions){
                float sqrd_dist = pow(x-i[0],2)+pow(y-i[1],2);
                dist.push_back(sqrd_dist);
            }
            float min_dist = *min_element(dist.begin(),dist.end());
            if(min_dist >= d){
                positions.push_back({x,y});
                Cells::Agent yeast(plate, *this, s, {x,y,s->init_pos_z});
                add_agent(plate, yeast);
                counter += 1;
            }
        }
    }

    void init_agents(Plate*& plate, std::vector<Cells::Species>& species){
        for(auto &i:species){
            init_positions(plate, &i);
        }
    }

    ///Agent lifecycle functions (in a different form for later parallelisation)

    //Feeding
    void feed(Plate*& plate, std::vector<Nutrient>& nutrients){
        #pragma omp parallel
        #pragma omp single
        {
            for(auto agent=agent_list.begin(); agent!=agent_list.end(); ++agent){
                #pragma omp task firstprivate(agent)
                {
                    //if agent is active: uptake and metabolism are unchanged, if it's in G0 state: both are weighted by g0_factor
                    float weight = (agent->state+!(agent->state)*agent->species->g0_factor);

                    auto idx=0;
                    for(auto &n:nutrients){
                        for(auto &i:agent->occupied_density_space){
                            //feeding
                            float uptake = std::min((agent->nutrient_uptake[idx])*(agent->species->uptake_eff[idx]),n.density_space.coeff(i));
                            uptake = weight*uptake;
                            agent->energy[idx] += uptake;
                            n.density_space.coeffRef(i) -=  uptake;
                        }
                        //metabolising
                        float burned_E = weight*agent->species->metab_E[idx];
                        agent->energy[idx] -= burned_E;
                        idx++;
                    }
                }
            }
            #pragma omp taskwait
        }
    }

    //State decision
    void decide_state(){
        #pragma omp parallel
        #pragma omp single
        {
            for(auto agent=agent_list.begin(); agent!=agent_list.end(); ++agent){
                #pragma omp task firstprivate(agent)
                //decide whether an agent should be active or in G0 state
                agent->state = (agent->energy > agent->species->g0_threshold);
            }
            #pragma omp taskwait
        }
    }

    //Division decision
    void can_divide(){
        #pragma omp parallel
        #pragma omp single
        {
            for(auto agent=agent_list.begin(); agent!=agent_list.end(); ++agent){
                #pragma omp task firstprivate(agent)
                //check energy level first
                agent->divide = (agent->energy > agent->species->div_threshold);
                if(!agent->divide){continue;}

                //check 'density' of grid cells occupied by agent
                //if all of them are 'full', then divide is 0, else 1
                agent->divide = false;
                for(auto &i:agent->occupied_uniform_grids){
                    auto range = agent_gridmap.equal_range(i);
                    auto density = std::distance(range.first,range.second);
                    agent->divide = agent->divide||(density < 24);
                }
            }
            #pragma omp taskwait
        }
    }

    //Cell division
    void cell_division(Plate*& plate, std::mt19937 &gen, std::normal_distribution<> &distr, int tries, float sigma, bool draw){
        for(auto agent=agent_list.begin(); agent!=agent_list.end(); ++agent){
            //skip agent if it's not qualified for division
            if(!agent->divide){continue;}

            float r = agent->species->cell_radius;
            float sqrd_r = agent->species->sqrd_radius;
            float a = agent->pos[2]-(plate->agar_height);
            if(a>=2*r){a=r-0.05;}
            float m = -1*a*inverse_square_root(4*r*r-a*a); //4*r*r = (2r)^2

            for(auto i=0;i<tries;++i){
                //create a possible position for daughter cell
                float x = (float) distr(gen);
                float y = (float) distr(gen);
                float z = (float) distr(gen) * sigma;

                float norm = inverse_square_root(x*x + y*y + z*z);
                x = x*norm;
                y = y*norm;
                z = z*norm;

//                //reject the ones that are outside of permitted area - classic rejection method
//                if((x>=plate->x-r) || (x<=0) ||
//                   (y>=plate->y-r) || (y<=0) ||
//                   (z>=plate->z-r) || (z<=plate->agar_height + r)){goto loop_end;}

                z = z+m*(-1*z+1);

                if(z<m){
                    z = m;
                    norm = inverse_square_root(x*x + y*y);
                    x = x*norm;
                    y = y*norm;
                }

                norm = inverse_square_root(x*x + y*y + z*z);
                x = x*norm*2*r + agent->pos[0];
                y = y*norm*2*r + agent->pos[1];
                z = z*norm*2*r + agent->pos[2];

                //if position is out of bounds, project it on the other side of the "division sphere"
                x = (x>=plate->x-r)*(x-2*(x-agent->pos[0]))+(x<=0)*(x+2*(agent->pos[0]-x))+(x<plate->x-r && x>0)*x;
                y = (y>=plate->y-r)*(y-2*(y-agent->pos[1]))+(y<=0)*(y+2*(agent->pos[1]-y))+(y<plate->y-r && y>0)*y;
                z = (z>=plate->z-r)*(z-2*(z-agent->pos[2]))+
                    (z<=plate->agar_height+r)*(z+2*(agent->pos[2]-z))+
                    (z<plate->z-r && z>plate->agar_height+r)*z;

                {
                //test each point, until a good one is found
                //see if occupied uniform grids here are "full"
                std::vector<int> occupied_uniform_grids = agent->occupied_grid_cells({0,0,plate->agar_height},conversion_factor,gridcell_size,gridmap_x,gridmap_y);
                for(auto grid:occupied_uniform_grids){
                   auto range = agent_gridmap.equal_range(grid);
                   auto density = std::distance(range.first,range.second);

                   if(density >= 24){goto loop_end;}

                   //if not, check the distance between the candidate coordinate and each agent already in this grid cell
                   for(auto& it=range.first; it!=range.second; ++it){
                        float sqrd_dist = pow(it->second->pos[0]-x,2) +
                                          pow(it->second->pos[1]-y,2) +
                                          pow(it->second->pos[2]-z,2);

                        if(sqrd_dist<(sqrd_r + it->second->species->sqrd_radius +
                           2*r*it->second->species->cell_radius)){goto loop_end;}
                    }
                }
                //if the distances are good, create daughter agent at (x,y,z)
                Cells::Agent new_agent(plate, *this, agent->species, {x,y,z});
                add_agent(plate, new_agent);
                if(draw){copy_positions(agent_list.begin());} //not so pretty to call it here as well :/ <- this may be redundant!
                goto next_agent;
                }

                loop_end:;
            }

            next_agent:;
        }
    }

    //Cell death
    void cell_death(Plate*& plate, std::vector<Nutrient>& nutrients, bool draw){
        for(auto agent=agent_list.begin(); agent!=agent_list.end(); ++agent){
            //skip agent if it's not qualified for death
            if(agent->energy >= agent->species->death_threshold){
                if(draw){copy_positions(agent);}
                continue;
            }

            //return the agent's respective energy value into the corresponding nutrient density matrix
            for(auto &i:agent->occupied_density_space){
                auto idx=0;
                for(auto &n:nutrients){
                    n.density_space.coeffRef(i) +=  agent->energy[idx];
                    idx++;
                }
            }
            //kill the agent
            delete_agent(plate, agent);
        }
    }
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
    bool is_draw;
    int steps;
    int resolution;
    int x;
    int y;
    int z;
    std::map<std::string,std::pair<const char*, const char*>> colormap;
    std::vector<mglEigenVec> nutrients;
    std::map<Cells::Species*,std::map<bool,std::array<mglData,3>>> agent_positions; //agent_positions in Cells will be linked to these

    yeastDraw(){}

    yeastDraw(std::ifstream& fin, std::string& line, Plate*& p, std::vector<Nutrient>& _nutrients, std::vector<Cells::Species>& species){
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        getline(fin,line,'\n');
        while(!line.empty()){
            std::stringstream s(line);
            s >> std::skipws >> is_draw;
            s >> std::skipws >> steps;
            s >> std::skipws >> resolution;
            getline(fin,line,'\n');
        }

        setup_colormap();
        setup_nutrient_container(p, _nutrients);
        setup_position_container(species);

        x = p->x;
        y = p->y;
        z = p->z;
    }

    ~yeastDraw(){}

    void setup_colormap(){
        colormap.emplace("red",std::make_pair("{xA93226}","{xE74C3C}"));
        colormap.emplace("orange",std::make_pair("{xE67E22}","{xF5B041}"));
        colormap.emplace("yellow",std::make_pair("{xF1C40F}","{xF7DC6F}"));
        colormap.emplace("green",std::make_pair("{x1E8449}","{x2ECC71}"));
        colormap.emplace("turquoise",std::make_pair("{x16A085}","{x76D7C4}"));
        colormap.emplace("blue",std::make_pair("{x2980B9}","{x5DADE2}"));
        colormap.emplace("purple",std::make_pair("{x7D3C98}","{xAF7AC5}"));
    }

    void setup_nutrient_container(Plate*& p, std::vector<Nutrient>& _nutrients){
        nutrients.reserve(_nutrients.size()); //Hell yeah!
        for(auto& n:_nutrients){
            nutrients.emplace_back(p, &n.density_space);
        }
    }

    void setup_position_container(std::vector<Cells::Species>& species){
        for(auto &i:species){
            std::array<mglData,3> positions;
            agent_positions[&i][0] = positions;
            agent_positions[&i][1] = positions;
        }
    }

    void link_agent_position(std::map<Cells::Species*,std::map<bool,std::array<std::vector<double>,3>>>& _agent_positions){
        for(auto &species_map:_agent_positions){
            for(auto &positions_map:species_map.second){
                for(size_t i=0;i<positions_map.second.size();++i){
                    agent_positions[species_map.first][positions_map.first][i].Link(&(positions_map.second[i])[0],positions_map.second[i].size());
                }
            }
        }
    }

    int Draw(mglGraph* gr){
     Agar(gr);
     YeastCells(gr);
     return 0;
    }

    void Agar(mglGraph* gr){
     //gr->Title("MathGL Demo");
     gr->Rotate(60,40);
     gr->Alpha(true);
     gr->Aspect(x,y,z);
     gr->SetRange('x',0,x); gr->SetRange('y',0,y); gr->SetRange('z',0,z);gr->SetRange('c',nutrients[0],true);
     gr->Axis("U");
     gr->Box();
     gr->Cloud(nutrients[0],"BbcyrR");//BBBBBB{xFFFFFF00}BbcyrR {xFFFFFF00}BbcyrR
     //gr->Colorbar("{xFFFFFF00}BbcyrR>I");
    }

    void YeastCells(mglGraph* gr){
        for(auto &species_map:agent_positions){
            //active agents
            gr->Dots(species_map.second[1][0],species_map.second[1][1],species_map.second[1][2],colormap[species_map.first->color].first,"size 10");
            //g0 agents
            gr->Dots(species_map.second[0][0],species_map.second[0][1],species_map.second[0][2],/*colormap[species_map.first->color].second*/"{x000000}","size 10");
        }
    }
};

class Statistics : public mglDraw{
public:
//    std::map<Cells::Species*,std::map<bool,std::vector<int>>> growth_curves;
//    std::map<Cells::Species*,std::map<bool,std::array<std::vector<double>,2>>> heatmaps_x_y;
    int species_num;
    std::map<Cells::Species*,std::map<bool,mglData>> growth_curves;
    std::map<Cells::Species*,std::map<bool,std::array<mglData,2>>> heatmaps_x_y;

    Statistics(){}

    //Initialize containers in constructor
    Statistics(int max_t, std::vector<Cells::Species>& species){
        species_num = species.size();
        for(auto &i:species){
            //growth curves
//            mglData population(max_t, 0);
            mglData population;
            population.Create(max_t);
            growth_curves[&i][0] = population;
            growth_curves[&i][1] = population;

            //heatmaps
            std::array<mglData,2> heatmap_positions;
            heatmaps_x_y[&i][0] = heatmap_positions;
            heatmaps_x_y[&i][1] = heatmap_positions;
        }
    }

    ~Statistics(){}

    void fill_growth_curve(int t, std::list<Cells::Agent>& agent_list){
        for(auto agent=agent_list.begin(); agent!=agent_list.end(); ++agent){
            growth_curves[agent->species][agent->state][t] += 1;
        }
    }

    //fill heatmaps_x_y(call it after the final positions have been copied)
    void fill_heatmaps(std::map<Cells::Species*,std::map<bool,std::array<std::vector<double>,3>>>& agent_positions){
        for(auto &species_map:agent_positions){
            for(auto &positions_map:species_map.second){
                //heatmaps
                heatmaps_x_y[species_map.first][positions_map.first][0] = agent_positions[species_map.first][positions_map.first][0];
                heatmaps_x_y[species_map.first][positions_map.first][1] = agent_positions[species_map.first][positions_map.first][1];
            }
        }
    }

    int Draw(mglGraph* gr){

        //growth curves
        std::string title;
        int counter = 1;
        int i = 0;
        //int x, y;
        for(auto &species_map:growth_curves){
            //active agents
            title = "Species-" + std::to_string(counter) + ", active agents";
            gr->SubPlot(species_num+1,2,i);  gr->Title(title.c_str());
            gr->SetRanges(0,species_map.second[1].GetNx(),0,species_map.second[1].Maximal());
            gr->Axis("U");
            gr->Box();
            gr->Plot(species_map.second[1]);

            //g0 agents
            i = i+species_num;
            title = "Species-" + std::to_string(counter) + ", G_0 agents";
            gr->SubPlot(species_num+1,2,i);  gr->Title(title.c_str());
            gr->SetRanges(0,species_map.second[0].GetNx(),0,species_map.second[0].Maximal());
            gr->Axis("U");
            gr->Box();
            gr->Plot(species_map.second[0]);

            counter ++;
            i++;
        }

//        //heatmaps
//        mglData z_values;
//        std::string title;
//        int counter = 1;
//        for(auto &species_map:heatmaps_x_y){
//            //active agents
//            int i = 0;
//            z_values.Create(species_map.second[1][0].GetNx());
//            title = "Species #" + std::to_string(counter) + ", active agents";
//            gr->SubPlot(species_num+1,2,i);  gr->Title(title.c_str());
//            gr->Dens(species_map.second[1][0],species_map.second[1][1], z_values);
//
//            //g0 agents
//            i = i+species_num;
//            z_values.Create(species_map.second[0][0].GetNx());
//            title = "Species #" + std::to_string(counter) + ", G_0 agents";
//            gr->SubPlot(species_num+1,2,i);  gr->Title(title.c_str());
//            gr->Dens(species_map.second[0][0],species_map.second[0][1], z_values);
//
//            counter += 1;
//        }

        return 0;
    }

};
//}

#endif // LIBRARIES_H_INCLUDED
