#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
using namespace std;

vector<vector<float>> read_in (string filename)
{
 ifstream fin;
 string line;
 string var;
 vector<float> var_category;
 vector<vector<float>> var_table;

 fin.open(filename, ios::in);
 if (fin.is_open())
  {
    while (fin.good())
    {
      // ignore header lines
      fin.ignore(numeric_limits<streamsize>::max(), '\n');
      // read in and put variables into a table
      getline(fin,line,'\n');
      while(!line.empty())
      {
       stringstream s(line);
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

  else cout << "Error in read_in.cpp: unable to open file";
  return var_table;
}
