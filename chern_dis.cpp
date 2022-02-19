#include <iostream>
#include <cmath>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <time.h>
#include <fstream>
#include <chrono>
using namespace std;
using namespace Eigen;
using namespace std::chrono;
ofstream myfile;



//FUNCTION DEFINITION
void Hamiltonian(MatrixXcf& Hpr, int size, float B, VectorXf Disorder,float W);
void Couple_matrix(MatrixXcf& Cqqprime, Vector2cf q, Vector2cf qPrime, MatrixXcf vecs, int size_sys);
int modulo(int x, int y){return (x%y + y)%y;}
int right(int x, int y, int N) {return N*y + modulo(x+1,N);}
int left(int x, int y, int N) {return N*y + modulo(x-1,N);}
int up(int x, int y, int N){return N*modulo(y+1, N) + x;}
int down(int x, int y, int N){return N*modulo(y-1, N) + x;}

int right_up(int x, int y, int N) {return N*modulo(y+1,N) + modulo(x+1,N);}
int right_down(int x, int y, int N) {return N*modulo(y-1,N) + modulo(x+1,N);}
int left_up(int x, int y, int N) {return N*modulo(y+1,N) + modulo(x-1,N);}
int left_down(int x, int y, int N) {return N*modulo(y-1,N) + modulo(x-1,N);}

//MAIN FUNCTION
int main(int argc, char* argv[]){
srand (static_cast <unsigned> (time(0)));
auto start = high_resolution_clock::now();
    int N = atoi(argv[1]);
    float B = atof(argv[2]);
    float W = atof(argv[3]);

   // Create The Hamiltonian
    VectorXf Dis(2*N*N); 
    Dis.setRandom();

    MatrixXcf H(2*N*N, 2*N*N);
    H.setZero();
    Hamiltonian(H, N, B, Dis, W);

   
   //find the eigenvector
    SelfAdjointEigenSolver<MatrixXcf> es(2*N*N);
    es.compute(H);
    
// define the vectors and coupling matrix

    Vector2cf q0,q1,q2,q3;
    q0 << 0.0,0.0;  
    q1 << 2*M_PI/(float) N,0.0;
    q2 << 2*M_PI/(float) N,2*M_PI/(float) N;
    q3 << 0.0,2*M_PI/(float) N;

    MatrixXcf C01(N*N,N*N), C12(N*N,N*N), C23(N*N,N*N), C30(N*N,N*N);
    C01.setZero();
    C12.setZero();
    C23.setZero();
    C30.setZero();


// Compute the coupling matrices elements
Couple_matrix(C01, q0, q1,es.eigenvectors(),N);
Couple_matrix(C12, q1, q2,es.eigenvectors(),N);
Couple_matrix(C23, q2, q3,es.eigenvectors(),N);
Couple_matrix(C30, q3, q0,es.eigenvectors(),N);


//Do the Dot Product
MatrixXcf DotProduct(N*N, N*N);
DotProduct = C01*C12*C23*C30;

//find the eigenvalues
    ComplexEigenSolver<MatrixXcf> es2(N*N);
    es2.compute(DotProduct);
//cout the sum of the phases divided by 2pi
   float Chern_number = 0.0;
    Chern_number = es2.eigenvalues().cwiseArg().sum();
   cout << -1.0*Chern_number/2/M_PI;
   auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start)/1000;

  cout << "Done! - " << duration.count() << " ms" <<endl;
}

void Couple_matrix(MatrixXcf& Cqqprime, Vector2cf q, Vector2cf qPrime, MatrixXcf vecs, int size_sys){
    for (int m=0; m <size_sys*size_sys;m++){
        
        for (int n=0; n <size_sys*size_sys;n++){
        
            for(int y=0; y<size_sys; y++){
                for( int x=0; x<size_sys;x++){
                    complex <float> qqpdot = (q(0)-qPrime(0))*(float)x + (q(1)-qPrime(1))*(float)y;
                    Cqqprime(m,n) += vecs.conjugate().col(m)(size_sys*y+x)*vecs.col(n)(size_sys*y+x)*exp(qqpdot*(complex <float>)1.0i) + vecs.conjugate().col(m)(size_sys*size_sys + size_sys*y+x)*vecs.col(n)(size_sys*size_sys + size_sys*y+x)*exp(qqpdot*(complex <float>)1.0i);
                }
            }
        }
    }
}


void Hamiltonian(MatrixXcf& Hpr, int size, float B, VectorXf Disorder, float W){
    for(int y=0; y<size; y++){
        for(int x=0;x<size;x++){
            int site = size*y+x;
            int r = right(x,y,size);
            int l = left(x,y,size);
           
            int u = up(x,y,size);
            int d = down(x,y,size);

            int ru = right_up(x,y,size);
            int rd = right_down(x,y,size);
            int lu = left_up(x,y,size);
            int ld = left_down(x,y,size);
            Hpr(site,site) = Disorder(site)*W/2.0;
            Hpr(size*size + site,size*size + site) = Disorder(size*size + site)*W/2.0;

            Hpr(r, site) = 1.0;
            Hpr(l, site) = 1.0;
            Hpr(u, site) = -1.0;
            Hpr(d, site) = -1.0;

            Hpr(size*size + r,size*size +  site) = -1.0;
            Hpr(size*size + l, size*size + site) = -1.0;
            Hpr(size*size + u, size*size + site) = 1.0;
            Hpr(size*size + d, size*size + site) = 1.0;

            Hpr(ru,size*size + site) = 0.5;
            Hpr(rd,size*size + site) = -0.5;
            Hpr(lu,size*size + site) = -0.5;
            Hpr(ld,size*size + site) = 0.5;

            Hpr(size*size +ru, site) = 0.5;
            Hpr(size*size +rd, site) = -0.5;
            Hpr(size*size +lu, site) = -0.5;
            Hpr(size*size +ld, site) = 0.5;
            
            Hpr(size*size + site, site) = -1.0i;
            Hpr(site, size*size + site) = 1.0i;

            Hpr(size*size + r, site) = -0.25*(B+1.0)*1.0i;
            Hpr(size*size + l, site) = -0.25*(B+1.0)*1.0i;
            Hpr(size*size + u, site) = -0.25*(B+1.0)*1.0i;
            Hpr(size*size + d, site) = -0.25*(B+1.0)*1.0i;

            Hpr(r,size*size + site) = 0.25*(B+1.0)*1.0i;      
            Hpr(l,size*size + site) = 0.25*(B+1.0)*1.0i;
            Hpr(u,size*size + site) = 0.25*(B+1.0)*1.0i;
            Hpr(d,size*size + site) = 0.25*(B+1.0)*1.0i;
        
        }
    }
}