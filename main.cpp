#include "master.h"
#include "main.h"

int main()
{
 //{--------------------------------------------------Create initial objects--------------------------------------------------
 Plate* plate = nullptr;
 yeastDraw* draw = nullptr;
 Cells cells;
 std::vector<Nutrient> nutrients;
 std::vector<Cells::Species> species;
 //}

 //{--------------------------------------------------Read in and initialize--------------------------------------------------
 read_in(&plate, &draw, nutrients, species, "3D_parameters.txt");
 mglFLTK gr(draw, "MathGL run test");
 //}

 //{--------------------------------------------------------Testing-----------------------------------------------------------
 for(auto x=3;x<6;x++){
    for(auto y=3;y<6;y++){
        for(auto z=plate->agar_height;z<plate->agar_height+3;z++){
            plate->occupancy_space.flip(x+(plate->x)*(y+(plate->y)*z));
        }
    }
 }

  for(auto x=3;x<6;x++){
    for(auto y=3;y<6;y++){
        for(auto z=plate->agar_height;z<plate->agar_height+3;z++){
            for(auto i : plate->neighbours[plate->diff_type]){
                if(plate->occupancy_space[x+(plate->x)*(y+(plate->y)*z)] && plate->occupancy_space[x+i[0]+(plate->x)*(y+i[1]+(plate->y)*(z+i[2]))]){
                    nutrients[0].laplace_mat.coeffRef(x+(plate->x)*(y+(plate->y)*z), x+i[0]+(plate->x)*(y+i[1]+(plate->y)*(z+i[2]))) += -1*nutrients[0].diff_const_cell;
                }
            }
        }
    }
 }

 for(int i=0; i<75; i++){
    nutrients[0].density_space -= 0.01*nutrients[0].laplace_mat*nutrients[0].density_space;
    gr.Update();
 }
 //}

 //{----------------------------------------------------Run FLTK window-------------------------------------------------------
 gr.Run();// This goes at the end of main!
 //}

 //{--------------------------------------------------------Cleanup-----------------------------------------------------------
 delete draw;
 delete plate;
 nutrients.erase(nutrients.begin());
 species.erase(species.begin());
 //}

 return 0;
}
