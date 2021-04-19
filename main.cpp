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
 read_in(plate, draw, nutrients, species, "3D_parameters.txt");
 mglFLTK gr(draw, "MathGL run test");
 //}

 //{--------------------------------------------------------Testing-----------------------------------------------------------
 for(auto x=3;x<6;x++){
    for(auto y=3;y<6;y++){
        for(auto z=plate->agar_height;z<plate->z-2;z++){
            plate->occupancy_space.flip(x+(plate->x)*(y+(plate->y)*z));
        }
    }
 }
 //}

 //{-------------------------------------------------------Diffusion----------------------------------------------------------
 update_laplace(plate, nutrients);
 diffusion(plate->diff_cnt, nutrients);
 gr.Update();
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
