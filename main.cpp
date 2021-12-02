#include "master.h"
#include "main.h"

int main()
{
 time_t start = time(0);

 srand(time(NULL));
 std::mt19937 gen(/*time(NULL)*/42);
// std::random_device rd{};
// std::mt19937 gen{rd()};
 std::normal_distribution<> distr(0,5);
 //{--------------------------------------------------Create initial objects--------------------------------------------------
 Plate* plate = nullptr;
 yeastDraw* draw = nullptr;
 mglFLTK* runtime_window = nullptr;
 Statistics* stats = nullptr;
 mglFLTK* stats_window = nullptr;
 std::vector<Nutrient> nutrients;
 std::vector<Cells::Species> species;
 //}

 //{--------------------------------------------------Read in and initialize--------------------------------------------------
 read_in(plate, draw, nutrients, species, "3D_parameters.txt");
 Cells cells(plate, species);
 plate->set_plate_dimensions(cells.gridcell_size);
 for(auto &nutrient:nutrients){nutrient.set_matrices(plate);}
 cells.init_agents(plate, gen, distr, species);
 if(draw->is_draw){
    runtime_window = new mglFLTK(draw, "Runtime window");
    stats = new Statistics(plate->t, species);
 }
 //}

 //{--------------------------------------------------------Main loop---------------------------------------------------------
 for(auto t=0; t<=plate->t; t++){
    std::cout<<t<<std::endl;
    bool draw_condition = draw->is_draw && t%draw->steps==0;
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
    cells.cell_division(plate, gen, distr, 500, 0.1, draw_condition);
    cells.cell_death(plate, nutrients, draw_condition);

    if(draw_condition){draw->link_agent_position(cells.agent_positions);}
    if(draw->is_draw){stats->fill_growth_curve(t, cells.agent_list);}
 //}

 //{-------------------------------------------------------Diffusion----------------------------------------------------------
    update_laplace(plate, nutrients);
    diffusion(plate->diff_cnt, nutrients);
    if(draw_condition){runtime_window->Update();}
 }
 //}

 //{---------------------------------------------Extract remaining statistics-------------------------------------------------
 if(draw->is_draw){
    stats->fill_heatmaps(cells.agent_positions);
    stats_window = new mglFLTK(stats, "Statistics");
 }
 //}

 //{----------------------------------------------------Run FLTK window-------------------------------------------------------
 double seconds_since_start = difftime( time(0), start);
 std::cout<<"Total runtime in seconds: "<<seconds_since_start<<std::endl;

 if(draw->is_draw){runtime_window->Run();}
 //}

 //{--------------------------------------------------------Cleanup-----------------------------------------------------------
 delete draw;
 delete runtime_window;
 delete stats;
 delete stats_window;
 //clean up agents
 delete plate;
 nutrients.clear();
 species.clear();
 //}

 return 0;
}
