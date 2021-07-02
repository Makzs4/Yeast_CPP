#include "master.h"

void diffusion(int& diff_cnt, std::vector<Nutrient>& nutrients, mglFLTK& gr){
    for(auto& n:nutrients){
        for(int t=0; t<diff_cnt; t++){
            n.density_space -= n.laplace_mat*n.density_space;
            //gr.Update();
        }
   }
}

void update_laplace(Plate*& plate, std::vector<Nutrient>& nutrients){
    for(auto& n:nutrients){
        n.laplace_triplet_list.erase(n.laplace_triplet_list.begin()+n.laplace_nutrient_size,n.laplace_triplet_list.end());
        for(auto x=0;x<plate->x;x++){
            for(auto y=0;y<plate->y;y++){
                for(auto z=0;z<plate->z;z++){
                    float counter = 0;
                    for(auto i : plate->neighbours[plate->diff_type]){
                        if(is_on_plate(x,y,z,plate->x,plate->y,plate->z,i) && plate->occupancy_space[x+(plate->x)*(y+(plate->y)*z)] && plate->occupancy_space[x+i[0]+(plate->x)*(y+i[1]+(plate->y)*(z+i[2]))]){
                        n.laplace_triplet_list.push_back(Eigen::Triplet<float>(x+(plate->x)*(y+(plate->y)*z),x+i[0]+(plate->x)*(y+i[1]+(plate->y)*(z+i[2])),-1*n.diff_const_cell));
                        counter += n.diff_const_cell;
                        }
                    }
                    n.laplace_triplet_list.push_back(Eigen::Triplet<float>(x+(plate->x)*(y+(plate->y)*z),x+(plate->x)*(y+(plate->y)*z),counter));
                }
            }
        }
        n.laplace_mat.setFromTriplets(n.laplace_triplet_list.begin(),n.laplace_triplet_list.end());
    }
}
