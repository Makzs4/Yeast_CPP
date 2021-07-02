#include "master.h"
#include "main.h"

int main()
{
 //{--------------------------------------------------Create initial objects--------------------------------------------------
 Plate* plate = nullptr;
 yeastDraw* draw = nullptr;
 std::vector<Nutrient> nutrients;
 std::vector<Cells::Species> species;
 //}

 //{--------------------------------------------------Read in and initialize--------------------------------------------------
 read_in(plate, draw, nutrients, species, "3D_parameters.txt");
 Cells cells(plate, species);
 cells.init_agents(plate, species);

 std::cout<<cells.agent_gridmap.size()<<std::endl;
 for(auto i:cells.agent_gridmap){
    std::cout<<i.first<<", "<<i.second->pos[0]<<std::endl;
 }

 std::cout<<cells.agent_list.size()<<std::endl;
 for(auto i:cells.agent_list){
    std::cout<<i.pos[0]<<std::endl;
 }

// Cells::Agent agent(plate, cells, &species[1], {15.6,29.2,9.3});
// mglFLTK gr(draw, "MathGL run test");
 //}

 //{--------------------------------------------------------Testing-----------------------------------------------------------
// int n = 5;
// for(auto x=(int)plate->x/2-n;x<(int)plate->x/2+n;x++){
//    for(auto y=(int)plate->y/2-n;y<(int)plate->y/2+n;y++){
//        for(auto z=plate->agar_height;z<plate->z-2;z++){
//            plate->occupancy_space.flip(x+(plate->x)*(y+(plate->y)*z));
//        }
//    }
// }

/// Random walk
// int x = (int)plate->x/2;
// int y = (int)plate->y/2;
// int z = plate->agar_height;
// int n = 5000;
// srand( /*(unsigned)time( NULL )*/1 );
// for(auto i=0;i<n;i++){
//    float val = (float) rand()/RAND_MAX;
//    //std::cout<<x<<" "<<y<<" "<<z<<" "<<val<<std::endl;
//    for(auto i=0;i<3;i++){
//        x = (val<0.16)*1 + (val>=0.16 && val<0.33)*-1 + x;
//        x = (x<0)*1 + (x>=plate->x)*-1 + x;
//        y = (val>=0.33 && val<0.5)*1 + (val>=0.5 && val<0.66)*-1 + y;
//        y = (y<0)*1 + (y>=plate->y)*-1 + y;
//        z = (val>=0.66 && val<0.83)*1 + (val>=0.83 && val<=1)*-1 + z;
//        z = (z<plate->agar_height)*1 + (z>=plate->z)*-1 + z;
//        int cell = x+(plate->x)*(y+(plate->y)*z);
//        plate->occupancy_space.set(cell);
//    }
// }

/// Piramid
// int n = 3;
// int m = 0;
// for(auto i=3; i>0; i--){
//    for(auto x=(int)plate->x/2-(i*n);x<(int)plate->x/2+(i*n);x++){
//        for(auto y=(int)plate->y/2-(i*n);y<(int)plate->y/2+(i*n);y++){
//            for(auto z=plate->agar_height+(m*n);z<plate->agar_height+((m+1)*n);z++){
//                plate->occupancy_space.flip(x+(plate->x)*(y+(plate->y)*z));
//            }
//        }
//    }
//    m++;
// }

/// Dome
// for(auto r=0; r<10; r++){
//    for(auto f=0; f<360; f++){
//        for(auto p=0; p<90 ;p++){
//            int x = (int)plate->x/2+r*cos(f*M_PI/180)*sin(p*M_PI/180);
//            int y = (int)plate->y/2+r*sin(f*M_PI/180)*sin(p*M_PI/180);
//            int z = plate->agar_height+r*cos(p*M_PI/180);
//            plate->occupancy_space.set(x+(plate->x)*(y+(plate->y)*z));
//        }
//    }
// }
 //}

// hashtable_test(plate);

 //{-------------------------------------------------------Diffusion----------------------------------------------------------
// update_laplace(plate, nutrients);
// diffusion(plate->diff_cnt, nutrients, gr);
// gr.Update();
 //}

 //{----------------------------------------------------Run FLTK window-------------------------------------------------------
// gr.Run();// This goes at the end of main!
 //}

 //{--------------------------------------------------------Cleanup-----------------------------------------------------------
 delete draw;
 delete plate;
 nutrients.clear();
 species.clear();
 //}

 return 0;
}
