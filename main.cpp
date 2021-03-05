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

 Eigen::MatrixXf m(2,2);
 m(0,0) = 3;
 m(1,0) = 2.5;
 m(0,1) = -1;
 m(1,1) = m(1,0) + m(0,1);

 std::cout << m << std::endl;
 std::cout << species[1].init_pos_z << std::endl;
 std::cout << plate->agar_height << std::endl;

 //{---------------------------------------------Cleanup----------------------------------------------
 delete draw;
 delete plate;
 nutrients.erase(nutrients.begin());
 species.erase(species.begin());
 //}

 return 0;
}
