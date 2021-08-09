#include "master.h"

void feed(Plate*& plate, std::vector<Nutrient>& nutrients, std::list<Cells::Agent>::iterator &agent){
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
