#include "master.h"
#include "main.h"

int main()
{
 //{---------------------------------------Create initial objects--------------------------------------
 Plate* plate = nullptr;
 Draw* draw = nullptr;
 Cells cells;
 std::vector<Nutrient> nutrients;
 std::vector<Cells::Species> species;
 //}

 //{-----------------------------------------Read in variables-----------------------------------------
 read_in(&plate, &draw, nutrients, species, "3D_parameters.txt");
 //}

 //{------------------------------------------Initialization-------------------------------------------
 run_mathgl_test(&plate, &nutrients[0].density_space);

 //}

 //{---------------------------------------------Cleanup----------------------------------------------
 delete draw;
 delete plate;
 nutrients.erase(nutrients.begin());
 species.erase(species.begin());
 //}

 return 0;
}
