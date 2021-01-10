#include <iostream>
#include <math.h>
#include <algorithm>
#include "main.h"

using namespace std;

int main()
{
 //{ Read in and declare variables
 vector<vector<float>>var_table = read_in("3D_parameters.txt");
 int n = 0;

 /// Plate parameters ///
 int dt = static_cast<int>(var_table[n][0]);n++; //time step
 int dx = static_cast<int>(var_table[n][0]);n++; //space step
 int max_x = static_cast<int>(var_table[n][0]);n++; //size of space along x axis
 int max_y = static_cast<int>(var_table[n][0]);n++; //size of space along y axis
 int max_z = static_cast<int>(var_table[n][0]);n++; //size of space along z axis
 int max_t = static_cast<int>(var_table[n][0]);n++; //iteration number (the time of the simulation)
 float init_nutrient = var_table[n][0];n++; //initial nutrient density
 float diff_co_nutrient_agar = var_table[n][0];n++; //diffusion coefficient: nutrient diffuses in agar
 float diff_co_nutrient_air = var_table[n][0];n++; //diffusion coefficient: nutrient diffuses in air
 float diff_co_nutrient_cell = var_table[n][0];n++; //diffusion coefficient: nutrient diffuses in cell
 float diff_co_oxigen_agar = var_table[n][0];n++; //diffusion coefficient: oxigen diffuses in agar
 float diff_co_oxigen_air = var_table[n][0];n++; //diffusion coefficient: oxigen diffuses in air
 float diff_co_oxigen_cell = var_table[n][0];n++; //diffusion coefficient: oxigen diffuses in cell
 float agar_height = var_table[n][0];n++; //height of culture medium
 int diff_cnt = static_cast<int>(var_table[n][0]);n++; //number of diffusion steps/agent life cycle

 /// Visualization parameters ///
 int isdraw = static_cast<int>(var_table[n][0]);n++;
 int plotstep = static_cast<int>(var_table[n][0]);n++; //visualization step
 int resolution = static_cast<int>(var_table[n][0]);n++; //resolution of the slice plot

 /// Physics parameters ///
 float r_cutoff = var_table[n][0];n++; //cutoff distance(Lennard-Jones)
 float epsilon = var_table[n][0];n++; //for the Lennard-Jones force
 float sigma = var_table[n][0];n++; // -||-
 float max_force = var_table[n][0];n++; // -||-
 int t_lj = static_cast<int>(var_table[n][0]);n++; // iteration number for Lennard-Jones

 /// Cell parameters ///
 vector<float> init_E = var_table[n];n++; //initial energy
 vector<float> nutrient_uptake = var_table[n];n++; //nutrient uptake
 vector<float> nutrient_uptake_eff = var_table[n];n++; //nutrient uptake efficiency
 vector<float> div_th = var_table[n];n++; //division threshold (in terms of energy)
 vector<float> div_r = var_table[n];n++; //division distance
 vector<float> death_th = var_table[n];n++; //death threshold (in terms of energy)
 vector<float> metab_E = var_table[n];n++; //energy usage by metabolism
 vector<float> colony_type = var_table[n];n++; //type of colony
 vector<float> g0_factor = var_table[n];n++; //penalty in g0 state (constant multiplier)
 vector<float> g0_th = var_table[n];n++; //g0 threshold (in terms of energy)
 vector<float> branchprob = var_table[n];n++; //branching probability
 vector<float> div_dir_dev = var_table[n];n++; //division direction deviation
 transform(div_dir_dev.begin(), div_dir_dev.end(), div_dir_dev.begin(), [](float i){return i/M_PI;});
 vector<float> init_cellnum = var_table[n];n++; //initial cell number
 vector<float> initdev = var_table[n];n++; //initial cell deviation
 vector<float> initpos_x = var_table[n];n++;
 vector<float> initpos_y = var_table[n];n++;
 vector<float> initpos_z(initpos_x.size());
 fill(initpos_z.begin(),initpos_z.end(), agar_height);
 //}

 //{ Initialization

 //}

 return 0;
}
