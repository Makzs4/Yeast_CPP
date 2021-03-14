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

 //}

 //std::cout << species[1].init_pos_z << std::endl;
 //std::cout << plate->dt << std::endl;

 //valarray_handling_test();
 eigen_handling_test();

 //{---------------------------------------------Cleanup----------------------------------------------
 delete draw;
 delete plate;
 nutrients.erase(nutrients.begin());
 species.erase(species.begin());
 //}

 return 0;
}
