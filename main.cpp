#include "master.h"
#include "main.h"

int main()
{
 srand(time(NULL));
 std::mt19937 gen(time(NULL));
// std::random_device rd{};
// std::mt19937 gen{rd()};
 std::normal_distribution<> distr(0,5);
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

 //{--------------------------------------------------------Main loop---------------------------------------------------------
 for(auto t=0; t<plate->t; t++){
    std::cout<<t<<std::endl;
    cells.clear_positions();
    for(auto it=cells.agent_list.begin();it!=cells.agent_list.end();++it){ //go through all agents
        feed(plate, nutrients, it);

        it->decide_state();
        it->can_divide(cells);

        if(it->state && it->divide){cell_division(plate, cells, it, gen, distr, 1000);}

        if(it->energy<(it->species->death_threshold)){ //see if agent should die
            cell_death(plate, nutrients, cells, it);
        }
        else{cells.copy_positions(it);}
    }
    draw->link_agent_position(cells.agent_positions);
 //}

 //{-------------------------------------------------------Diffusion----------------------------------------------------------
    update_laplace(plate, nutrients);
    diffusion(plate->diff_cnt, nutrients, gr);
    if(t%draw->steps==0){gr.Update();}
 }
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
