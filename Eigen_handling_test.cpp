#include "master.h"

void eigen_handling_test(){

int width = 2, depth = 3, height = 2;
int mat_size = width*depth*height;

std::vector<Eigen::Triplet<double>> adj_triplet_list;
adj_triplet_list.reserve(6);
Eigen::SparseMatrix<double> adj_mat(mat_size,mat_size);

std::vector<Eigen::Triplet<double>> deg_triplet_list;
deg_triplet_list.reserve(1);
Eigen::SparseMatrix<double> deg_mat(mat_size,mat_size);

//for(int i=0; i<mat_size; i++){
//    if(i%width!=(width-1)){adj_triplet_list.push_back(Eigen::Triplet<double>(i,i+1,1));}
//    if(i%width!=0){adj_triplet_list.push_back(Eigen::Triplet<double>(i,i-1,1));}
//    if((i/width)%depth!=(depth-1)){adj_triplet_list.push_back(Eigen::Triplet<double>(i,i+width,1));}
//    if((i/width)%depth!=0){adj_triplet_list.push_back(Eigen::Triplet<double>(i,i-width,1));}
//    if(i/(width*depth)!=(height-1)){adj_triplet_list.push_back(Eigen::Triplet<double>(i,i+(width*depth),1));}
//    if(i/(width*depth)!=0){adj_triplet_list.push_back(Eigen::Triplet<double>(i,i-(width*depth),1));}
//}

for(int i=0; i<width; i++){
    for(int j=0; j<depth; j++){
        for(int k=0; k<height; k++){
                int counter = 0;
                //Adjacency matrix
                if(i!=(width-1)){adj_triplet_list.push_back(Eigen::Triplet<double>(i+width*(j+depth*k),(i+1)+width*(j+depth*k),1));counter++;}
                if(i!=0){adj_triplet_list.push_back(Eigen::Triplet<double>(i+width*(j+depth*k),(i-1)+width*(j+depth*k),1));counter++;}
                if(j!=(depth-1)){adj_triplet_list.push_back(Eigen::Triplet<double>(i+width*(j+depth*k),i+width*((j+1)+depth*k),1));counter++;}
                if(j!=0){adj_triplet_list.push_back(Eigen::Triplet<double>(i+width*(j+depth*k),i+width*((j-1)+depth*k),1));counter++;}
                if(k!=(height-1)){adj_triplet_list.push_back(Eigen::Triplet<double>(i+width*(j+depth*k),i+width*(j+depth*(k+1)),1));counter++;}
                if(k!=0){adj_triplet_list.push_back(Eigen::Triplet<double>(i+width*(j+depth*k),i+width*(j+depth*(k-1)),1));counter++;}

                //Degree matrix
                deg_triplet_list.push_back(Eigen::Triplet<double>(i+width*(j+depth*k),i+width*(j+depth*k),counter));
        }
    }
}

adj_mat.setFromTriplets(adj_triplet_list.begin(),adj_triplet_list.end()); //Adjacency matrix
deg_mat.setFromTriplets(deg_triplet_list.begin(),deg_triplet_list.end()); //Degree matrix

Eigen::SparseMatrix<double> lap_mat(mat_size,mat_size);
lap_mat = deg_mat - adj_mat;

std::cout << "Adjacency matrix" << "                "<< "Degree matrix" << "                   "<< "Laplacian matrix" << std::endl;
for(int i=0;i<mat_size;i++){
    for(int j=0;j<mat_size;j++){
        std::cout << adj_mat.coeff(i,j) << ' ';
    }

    std::cout << "        ";

    for(int j=0;j<mat_size;j++){
        std::cout << deg_mat.coeff(i,j) << ' ';
    }

    std::cout << "        ";

    for(int j=0;j<mat_size;j++){
        std::cout << lap_mat.coeff(i,j) << ' ';
    }
    std::cout << std::endl;
}
}
