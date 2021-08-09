#include "master.h"

//assume only a dividing agent gets here: check goes in main

void cell_division(Plate*& plate, Cells &cells, std::list<Cells::Agent>::iterator &agent, int points, int tries){
    //create a number of possible coordinates for daughter cell
    std::cout<<"Energy: "<<agent->energy[0]<<std::endl;

    auto total_points = points * tries;
    for(auto i=0;i<total_points;++i){
        float x = (float) rand()/RAND_MAX;
        float y = (float) rand()/RAND_MAX;
        float z = (float) rand()/RAND_MAX;
        float norm = inverse_square_root(x*x + y*y + z*z);

        x = x*norm*2*agent->species->cell_radius + agent->pos[0];
        y = y*norm*2*agent->species->cell_radius + agent->pos[1];
        z = z*norm*2*agent->species->cell_radius + agent->pos[2];

        std::cout<<"Candidate coordinates: "<<x<<", "<<y<<", "<<z<<std::endl;

        //reject the ones that are outside of permitted area
        if((x<=plate->x-agent->species->cell_radius) && (x>=0) &&
           (y<=plate->y-agent->species->cell_radius) && (y>=0) &&
           (z<=plate->z-agent->species->cell_radius) && (z>=plate->agar_height+agent->species->cell_radius)){

               std::cout<<"Point is in permitted area"<<std::endl;

               //test each point, until a good one is found
               //see if occupied uniform grids here are "full"
               std::vector<int> occupied_uniform_grids = agent->occupied_grid_cells({0,0,plate->agar_height},cells.conversion_factor,cells.gridcell_size,cells.gridmap_x,cells.gridmap_y);
               for(auto grid:occupied_uniform_grids){
                   auto range = cells.agent_gridmap.equal_range(grid);
                   auto density = std::distance(range.first,range.second);

                   std::cout<<"Local density: "<<density<<std::endl;
                   //if not, check the distance between the candidate coordinate and each agent already in this grid cell
                   if(density < 33){
                        bool suitable = false; //logic is wrong/maybe not
                        for(auto& it=range.first; it!=range.second; ++it){
                            float sqrd_dist = pow(it->second->pos[0]-x,2) +
                                              pow(it->second->pos[1]-y,2) +
                                              pow(it->second->pos[2]-z,2);

                            std::cout<<"Squared distance: "<<sqrd_dist<<std::endl;
                            std::cout<<"Permitted sqrd distance: "<<agent->species->sqrd_radius<<std::endl;
                            suitable = suitable || (sqrd_dist>=agent->species->sqrd_radius);
                            }

                        std::cout<<"Suitable: "<<suitable<<std::endl;
                        //if the distances are good, create daughter agent at (x,y,z)
                        if(suitable){
                            Cells::Agent new_agent(plate, cells, agent->species, {x,y,z});
                            cells.add_agent(plate, new_agent);
                            cells.copy_positions(cells.agent_list.begin()); //not so pretty to call it here as well :/ //HERE!
                            std::cout<<"Division happened!"<<std::endl;
                            return;
                        }
                   }
               }
           }
    }
    std::cout<<"No suitable place found, no division"<<std::endl;
}
