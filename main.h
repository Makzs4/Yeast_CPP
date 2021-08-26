#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

void read_in (Plate*& plate, yeastDraw*& draw, std::vector<Nutrient>& nutrients, std::vector<Cells::Species>& species, std::string filename);
void diffusion(int& diff_cnt, std::vector<Nutrient>& nutrients, mglFLTK& gr);
void update_laplace(Plate*& plate, std::vector<Nutrient>& nutrients);
void feed(Plate*& plate, std::vector<Nutrient>& nutrients, std::list<Cells::Agent>::iterator &agent);
void cell_death(Plate*& plate, std::vector<Nutrient>& nutrients, Cells &cells,std::list<Cells::Agent>::iterator &agent);
void cell_division(Plate*& plate, Cells &cells, std::list<Cells::Agent>::iterator &agent, std::mt19937 &gen, std::normal_distribution<> &distr, int tries);

//might delete later
void piramid_diffusion(Plate*& plate);
void dome_diffusion(Plate*& plate);
void random_walk_diffusion(Plate*& plate);

// delete later!
void valarray_handling_test();
void eigen_handling_test();
int run_mathgl_test(Plate** plate, Eigen::SparseVector<float> *d);
void hashtable_test(Plate*& plate);

#endif // MAIN_H_INCLUDED
