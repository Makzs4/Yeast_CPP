#include "master.h"

void cell_death(Plate*& plate, std::vector<Nutrient>& nutrients, Cells &cells, std::list<Cells::Agent>::iterator &agent){
    //return the agent's respective energy value into the corresponding nutrient density matrix
    for(auto &i:agent->occupied_density_space){
        auto idx=0;
        for(auto &n:nutrients){
            n.density_space.coeffRef(i) +=  agent->energy[idx];
            idx++;
        }
    }
    //kill the agent
    cells.delete_agent(plate, agent);
}
