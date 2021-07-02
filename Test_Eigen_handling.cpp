#include "master.h"

void eigen_handling_test(){

int width = 2, depth = 3, height = 2;
int mat_size = width*depth*height;

std::vector<Eigen::Triplet<float>> laplace_triplet_list;
laplace_triplet_list.reserve(7);
Eigen::SparseMatrix<float> laplace_mat(mat_size,mat_size);

Eigen::SparseVector<float> U(mat_size);
for(int i=0; i<4; i++){
    U.insert(i) = 1;
}

//for(int i=0; i<mat_size; i++){
//    if(i%width!=(width-1)){laplace_triplet_list.push_back(Eigen::Triplet<double>(i,i+1,1));}
//    if(i%width!=0){laplace_triplet_list.push_back(Eigen::Triplet<double>(i,i-1,1));}
//    if((i/width)%depth!=(depth-1)){laplace_triplet_list.push_back(Eigen::Triplet<double>(i,i+width,1));}
//    if((i/width)%depth!=0){laplace_triplet_list.push_back(Eigen::Triplet<double>(i,i-width,1));}
//    if(i/(width*depth)!=(height-1)){laplace_triplet_list.push_back(Eigen::Triplet<double>(i,i+(width*depth),1));}
//    if(i/(width*depth)!=0){laplace_triplet_list.push_back(Eigen::Triplet<double>(i,i-(width*depth),1));}
//}

for(int x=0; x<width; x++){
    for(int y=0; y<depth; y++){
        for(int z=0; z<height; z++){
                int counter = 0;
                //Adjacency matrix
                if(x!=(width-1)){laplace_triplet_list.push_back(Eigen::Triplet<float>(x+width*(y+depth*z),(x+1)+width*(y+depth*z),-1));counter++;}
                if(x!=0){laplace_triplet_list.push_back(Eigen::Triplet<float>(x+width*(y+depth*z),(x-1)+width*(y+depth*z),-1));counter++;}
                if(y!=(depth-1)){laplace_triplet_list.push_back(Eigen::Triplet<float>(x+width*(y+depth*z),x+width*((y+1)+depth*z),-1));counter++;}
                if(y!=0){laplace_triplet_list.push_back(Eigen::Triplet<float>(x+width*(y+depth*z),x+width*((y-1)+depth*z),-1));counter++;}
                if(z!=(height-1)){laplace_triplet_list.push_back(Eigen::Triplet<float>(x+width*(y+depth*z),x+width*(y+depth*(z+1)),-1));counter++;}
                if(z!=0){laplace_triplet_list.push_back(Eigen::Triplet<float>(x+width*(y+depth*z),x+width*(y+depth*(z-1)),-1));counter++;}

                laplace_triplet_list.push_back(Eigen::Triplet<float>(x+width*(y+depth*z),x+width*(y+depth*z),counter));
        }
    }
}

laplace_mat.setFromTriplets(laplace_triplet_list.begin(),laplace_triplet_list.end()); //Laplace operator matrix

std::cout << "Laplacian operator matrix" << std::endl;
for(int i=0;i<mat_size;i++){
    for(int j=0;j<mat_size;j++){
        std::cout << laplace_mat.coeff(i,j) << ' ';
    }
    std::cout << std::endl;
}

float c = 0.1;
for(int i=0; i<50; i++){
    U -= c*laplace_mat*U; //FUCK YEAH!
}
for(int i=0; i<mat_size; i++){
    std::cout<<U.coeff(i)<<std::endl;
}

}
