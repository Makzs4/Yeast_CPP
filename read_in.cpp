#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

std::vector<std::vector<float>> read_in (std::string filename)
{
 std::ifstream fin;
 std::string line;
 std::string var;
 std::vector<float> var_category;
 std::vector<std::vector<float>> var_table;

 fin.open(filename, std::ios::in);
 if (fin.is_open())
  {
    while (fin.good())
    {
      // ignore header lines
      fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      // read in and put variables into a table
      getline(fin,line,'\n');
      while(!line.empty())
      {
       std::stringstream s(line);
       getline(s,var,';'); //ignore first segment up to ';'
       while(getline(s,var,';'))
       {
        var_category.push_back(stof(var));
       }
       var_table.push_back(var_category);
       var_category.clear();
       getline(fin,line,'\n');
      }
    }
    fin.close();
  }

  else std::cout << "Error in read_in.cpp: unable to open file";
  return var_table;
}
