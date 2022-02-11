#include <iostream>
#include <cmath>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
using namespace std;
using namespace Eigen;





//FUNCTION DEFINITION
void Hamiltonian(MatrixXcf& Hpr, int size);
void Couple_matrix(MatrixXcf& Cqqprime, Vector2f q, Vector2f qPrime, MatrixXf vecs, int size_sys);
int modulo(int x, int y){return (x%y + y)%y;}
int right(int x, int y, int N) {return N*y + modulo(x+1,N);}
int left(int x, int y, int N) {return N*y + modulo(x-1,N);}
int up(int x, int y, int N){return N*modulo(y+1, N) + x;}
int down(int x, int y, int N){return N*modulo(y-1, N) + x;}



int main(int argc, char* argv[]){

    int N = atoi(argv[1]);

   // Create The Hamiltonian

// BLA BLA BLA

   //find the eigenvector

// BLA BLA BLA

// define the vectors and coupling matrix

    Vector2f q0,q1,q2,q3;
    q0 << 0.0,0.0;
    q1 << 2*M_PI/(float) N,0.0;
    q2 << 2*M_PI/(float) N,2*M_PI/(float) N;
    q3 << 0.0,2*M_PI/(float) N;

    MatrixXcf C01(N*N), C12(N*N), C23(N*N), C30(N*N);
    C01.setZero();
    C12.setZero();
    C23.setZero();
    C30.setZero();
cout<<exp(2*M_PI/(float)N*1.0i);
// Compute the coupling matrices elements



//Do the Dot Product

//find the eigenvalues

//cout the sum of the phases divided by 2pi

}

void Couple_matrix(MatrixXcf& Cqqprime, Vector2f q, Vector2f qPrime, MatrixXf vecs, int size_sys){
    for (int m=0; m <pow(size_sys,2);m++){
        for (int n=0; m <pow(size_sys,2);n++){
            for(int y=0; y<size_sys; y++){
                for( int x=0; x<size_sys;x++){
                    float qqpdot = (q(0)-qPrime(0))*x+(q(1)-qPrime(1))*y;
                    Cqqprime(m,n) += conj(vecs.col(m)(size_sys*y+x))*vecs.col(m)(size_sys*y+x);
                }
            }
        }
    }
}