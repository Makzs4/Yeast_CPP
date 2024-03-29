#include "master.h"

void read_in (Plate*& plate, yeastDraw*& draw, std::vector<Nutrient>& nutrients, std::vector<Cells::Species>& species, std::string filename)
{
 std::ifstream fin;
 std::string line;
 fin.open(filename, std::ios::in);
 if (fin.is_open())
  {
    while (getline(fin,line,'\n'))
     {
        if(line == "VISUALISATION"){
           draw = new yeastDraw(fin, line, plate, nutrients, species);
        }
        else if(line == "PLATE"){
           plate = new Plate(fin, line);
        }
        else if(line == "NUTRIENT"){
           nutrients.push_back(Nutrient(fin, line, plate));
        }
        else if(line == "CELLS"){
           species.push_back(Cells::Species(fin, line, plate, nutrients.size()));
        }
     }
    }
 else std::cout << "Error in read_in.cpp: unable to open file";
 fin.close();
}
