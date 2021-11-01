#include "master.h"
#include "main.h"

int main()
{
 time_t start = time(0);

 srand(time(NULL));
 std::mt19937 gen(time(NULL));
// std::random_device rd{};
// std::mt19937 gen{rd()};
 std::normal_distribution<> distr(0,5);
 //{--------------------------------------------------Create initial objects--------------------------------------------------
 Plate* plate = nullptr;
 yeastDraw* draw = nullptr;
 mglFLTK* gr = nullptr;
 std::vector<Nutrient> nutrients;
 std::vector<Cells::Species> species;
 //}

 //{--------------------------------------------------Read in and initialize--------------------------------------------------
 read_in(plate, draw, nutrients, species, "3D_parameters.txt");
 Cells cells(plate, species);
 plate->set_plate_dimensions(cells.gridcell_size);
 for(auto &nutrient:nutrients){nutrient.set_matrices(plate);}
 cells.init_agents(plate, species);
 if(draw->is_draw){gr = new mglFLTK(draw, "MathGL run test");}
 //}

 //{--------------------------------------------------------Main loop---------------------------------------------------------
 for(auto t=0; t<plate->t; t++){
    std::cout<<t<<std::endl;
    cells.clear_positions();

////life-cycle functions
//for(auto it=cells.agent_list.begin();it!=cells.agent_list.end();++it){
//    feed(plate, nutrients, it);
//
//    it->decide_state();
//    it->can_divide(cells);
//
//    if(it->state && it->divide){cell_division(plate, cells, it, gen, distr, 1000);}
//
//    if(it->energy<(it->species->death_threshold)){ //see if agent should die
//        cell_death(plate, nutrients, cells, it);
//    }
//    else if(draw->is_draw){cells.copy_positions(it);}
//}

    //life-cycle functions
    cells.feed(plate, nutrients);
    cells.decide_state();
    cells.can_divide();
    cells.cell_division(plate, gen, distr, 500, 0.5, draw->is_draw);
    cells.cell_death(plate, nutrients, draw->is_draw);

    if(draw->is_draw){draw->link_agent_position(cells.agent_positions);}
 //}

 //{-------------------------------------------------------Diffusion----------------------------------------------------------
    update_laplace(plate, nutrients);
    diffusion(plate->diff_cnt, nutrients);
    if(draw->is_draw && t%draw->steps==0){gr->Update();}
 }
 //}

 //{----------------------------------------------------Run FLTK window-------------------------------------------------------
 if(draw->is_draw){gr->Run();}// This goes at the end of main!

 int sum = 0;
 for(auto &i:cells.agent_gridmap){
    int cnt = cells.agent_gridmap.count(i.first);
    if(cnt!=0){sum += cnt;}
 }
 float mean = (float)sum/(float)cells.agent_gridmap.size();
 std::cout<<"Mean grid density in hash table: "<<mean<<std::endl;
 std::cout<<"Number of agents: "<<cells.agent_list.size()<<std::endl;

 double seconds_since_start = difftime( time(0), start);
 std::cout<<"Total runtime in seconds: "<<seconds_since_start<<std::endl;
 //}

 //{--------------------------------------------------------Cleanup-----------------------------------------------------------
 delete draw;
 delete gr;
 //clean up agents
 delete plate;
 nutrients.clear();
 species.clear();
 //}

 return 0;
}
