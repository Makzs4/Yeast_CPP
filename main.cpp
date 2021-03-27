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

//std::cout << nutrients[1].laplace_mat.size() << std::endl;
//std::cout << nutrients[0].name << std::endl;
//std::cout << nutrients[0].laplace_mat.coeff(0,1) << std::endl;

//eigen_handling_test();

 //{---------------------------------------------Cleanup----------------------------------------------
 delete draw;
 delete plate;
 nutrients.erase(nutrients.begin());
 species.erase(species.begin());
 //}

 return 0;
}
