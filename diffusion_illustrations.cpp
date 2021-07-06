#include "master.h"

void piramid_diffusion(Plate*& plate){
    int n = 3;
    int m = 0;
    for(auto i=3; i>0; i--){
        for(auto x=(int)plate->x/2-(i*n);x<(int)plate->x/2+(i*n);x++){
            for(auto y=(int)plate->y/2-(i*n);y<(int)plate->y/2+(i*n);y++){
                for(auto z=plate->agar_height+(m*n);z<plate->agar_height+((m+1)*n);z++){
                    plate->occupancy_space.flip(x+(plate->x)*(y+(plate->y)*z));
                }
            }
        }
        m++;
    }
}

void dome_diffusion(Plate*& plate){
    for(auto r=0; r<10; r++){
        for(auto f=0; f<360; f++){
            for(auto p=0; p<90 ;p++){
                int x = (int)plate->x/2+r*cos(f*M_PI/180)*sin(p*M_PI/180);
                int y = (int)plate->y/2+r*sin(f*M_PI/180)*sin(p*M_PI/180);
                int z = plate->agar_height+r*cos(p*M_PI/180);
                plate->occupancy_space.set(x+(plate->x)*(y+(plate->y)*z));
            }
        }
    }
}

void random_walk_diffusion(Plate*& plate){
    int x = (int)plate->x/2;
    int y = (int)plate->y/2;
    int z = plate->agar_height;
    int n = 5000;
    //srand( /*(unsigned)time( NULL )*/1 );
    for(auto i=0;i<n;i++){
        float val = (float) rand()/RAND_MAX;
        for(auto i=0;i<3;i++){
            x = (val<0.16)*1 + (val>=0.16 && val<0.33)*-1 + x;
            x = (x<0)*1 + (x>=plate->x)*-1 + x;
            y = (val>=0.33 && val<0.5)*1 + (val>=0.5 && val<0.66)*-1 + y;
            y = (y<0)*1 + (y>=plate->y)*-1 + y;
            z = (val>=0.66 && val<0.83)*1 + (val>=0.83 && val<=1)*-1 + z;
            z = (z<plate->agar_height)*1 + (z>=plate->z)*-1 + z;
            int cell = x+(plate->x)*(y+(plate->y)*z);
            plate->occupancy_space.set(cell);
        }
    }
}
