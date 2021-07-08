#include "master.h"
#include "main.h"

int main()
{
 srand(time(NULL));
 //{--------------------------------------------------Create initial objects--------------------------------------------------
 Plate* plate = nullptr;
 yeastDraw* draw = nullptr;
 std::vector<Nutrient> nutrients;
 std::vector<Cells::Species> species;
 //}

 //{--------------------------------------------------Read in and initialize--------------------------------------------------
 read_in(plate, draw, nutrients, species, "3D_parameters.txt");
 Cells cells(plate, species);
 cells.init_agents(plate, species);
// mglFLTK gr(draw, "MathGL run test");
 //}

 //{--------------------------------------------------------Testing-----------------------------------------------------------
 // Cells::Agent agent(plate, cells, &species[1], {15.6,29.2,9.3});
 std::cout<<cells.agent_gridmap.size()<<std::endl;
 for(auto i:cells.agent_gridmap){
    std::cout<<i.first<<", "<<i.second<<std::endl;
 }

 std::cout<<cells.agent_list.size()<<std::endl;
 for(auto& i:cells.agent_list){ //I AM A GOD
    std::cout<<&i<<std::endl;
 }

 auto a = cells.agent_list.begin(); //dummy example
 a->energy = 0.1;
 std::cout<<"Deleted agent: "<<&(*a)<<std::endl;

 for(auto it=cells.agent_list.begin();it!=cells.agent_list.end();++it){ //go through all agents
    if(it->energy<(it->species->death_threshold)){ //see if agent should die
        for(auto j:it->occupied_uniform_grids){ //go through every uniform grid it occupies
            auto k = cells.agent_gridmap.equal_range(j); //find the grids in the agent_map
            for(auto l=k.first;l!=k.second;++l){ //go through all occupied grids and delete agent to be killed
                if(l->second == &(*it)){
                    //delete from hashmap
                    cells.agent_gridmap.erase(l);
                    break;
                }
            }
        }
        //delete from list
        it = cells.agent_list.erase(it);
        --it;
    }
    if(it==cells.agent_list.begin()){
        Cells::Agent yeast(plate, cells, &species[1], {15.6,29.2,species[1].init_pos_z});
        cells.add_agent(yeast);
        std::cout<<"Added agent: "<<&(*cells.agent_list.begin())<<std::endl;
    }
 }

 std::cout<<cells.agent_gridmap.size()<<std::endl;
 for(auto i:cells.agent_gridmap){
    std::cout<<i.first<<", "<<i.second<<std::endl;
 }

 std::cout<<cells.agent_list.size()<<std::endl;
 for(auto& i:cells.agent_list){
    std::cout<<&i<<std::endl;
 }
 //}

 //{-------------------------------------------------------Diffusion----------------------------------------------------------
// update_laplace(plate, nutrients);
// diffusion(plate->diff_cnt, nutrients, gr);
// gr.Update();
 //}

 //{----------------------------------------------------Run FLTK window-------------------------------------------------------
// gr.Run();// This goes at the end of main!
 //}

 //{--------------------------------------------------------Cleanup-----------------------------------------------------------
 delete draw;
 //clean up agents
 delete plate;
 nutrients.clear();
 species.clear();
 //}

 return 0;
}
