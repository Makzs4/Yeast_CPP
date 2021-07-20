#include "master.h"

void feed(Plate*& plate, std::vector<Nutrient>& nutrients, std::list<Cells::Agent>::iterator &agent){
    for(auto &i:agent->occupied_density_space){
        for(auto &n:nutrients){
            float uptake = std::min((agent->nutrient_uptake)*(agent->species->uptake_eff),n.density_space.coeff(i));
            agent->energy += uptake;
            n.density_space.coeffRef(i) -=  uptake;
        }
    }
}
