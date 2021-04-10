#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

void read_in (Plate** plate, yeastDraw** draw, std::vector<Nutrient>& nutrients, std::vector<Cells::Species>& species, std::string filename);

// delete later!
void valarray_handling_test();
void eigen_handling_test();
int run_mathgl_test(Plate** plate, Eigen::SparseVector<float> *d);

#endif // MAIN_H_INCLUDED
