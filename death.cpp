#include "master.h"

void cell_death(Plate*& plate, std::vector<Nutrient>& nutrients, Cells &cells, std::list<Cells::Agent>::iterator &agent){
    //return the agent's respective energy value into the corresponding nutrient density matrix
    auto num_girds = agent->occupied_density_space.size();
    for(auto &i:agent->occupied_density_space){
        auto idx=0;
        for(auto &n:nutrients){
            n.density_space.coeffRef(i) +=  static_cast<float>(agent->energy[idx])/num_girds;
            idx++;
        }
    }
    //kill the agent
    cells.delete_agent(plate, agent);
}
