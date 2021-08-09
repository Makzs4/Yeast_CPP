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
 mglFLTK gr(draw, "MathGL run test");
 //}

 //{--------------------------------------------------------Testing-----------------------------------------------------------
 std::cout<<cells.agent_gridmap.size()<<std::endl;
 for(auto i:cells.agent_gridmap){
    std::cout<<i.first<<", "<<i.second<<std::endl;
 }

 std::cout<<cells.agent_list.size()<<std::endl;
 for(auto& i:cells.agent_list){ //I AM A GOD
    std::cout<<&i<<std::endl;
 }

 auto a = cells.agent_list.begin(); //dummy example
 a->energy[0] = 0.1;
 std::cout<<"Deleted agent: "<<&(*a)<<std::endl;
 auto b = std::next(cells.agent_list.begin());
 b->energy[0] = 10;

 cells.clear_positions();
 for(auto it=cells.agent_list.begin();it!=cells.agent_list.end();++it){ //go through all agents
    feed(plate, nutrients, it);

    it->decide_state();
    it->can_divide(cells);

    if(it->state && it->divide){cell_division(plate, cells, it, 100, 10);}

    if(it->energy<(it->species->death_threshold)){ //see if agent should die
        cell_death(plate, nutrients, cells, it);
    }
    else{cells.copy_positions(it);}
 }
 draw->link_agent_position(cells.agent_positions);

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
// piramid_diffusion(plate);
 update_laplace(plate, nutrients);
 diffusion(plate->diff_cnt, nutrients, gr);
 gr.Update();
 //}

 //{----------------------------------------------------Run FLTK window-------------------------------------------------------
 gr.Run();// This goes at the end of main!
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
