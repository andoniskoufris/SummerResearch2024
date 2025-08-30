#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <immintrin.h>

#define ALIGN 64

// this function adds a number to a (flattened) matrix of size N by 2, such as the array that will store positions. This will be used in the integration step of Runge-Kutta when I will need to calculate \dot{r}_i at the current positions plus some additional number, for the calculation of each f_0, f_1, etc.
void matrix_add_number(double* array, double constant, int N) {
    
    for (int i=0; i<2*N; i++) {
        array[i] = array[i] + constant;
    }

}


// this function will allow me to add together two (flattened) arrays of size N by 2, such as when I need to do some things
void matrix_add_array(double* array_out, double* array_add, int N) {

}


// -------------------------------------- CALCULATING NUMBERS FOR SIMULATION -------------------------------------- //

// this function calculates x_{ij} and stores it in a matrix called X_ij 
void x_ij(double* X_ij, double* pos_Real, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {
        X_ij[NumVortices*i + i]=0;
        for (int j=0; j<i; j++) {
            X_ij[NumVortices*i + j] = pos_Real[2*i + 0] - pos_Real[2*j + 0];
            X_ij[NumVortices*j + i] = X_ij[NumVortices*i + j];
        }
    }

}


// this function calculates y_{ij} and stores it in a matrix called Y_ij 
void y_ij(double* Y_ij, double* pos_Real, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {
        Y_ij[NumVortices*i + i]=0;
        for (int j=0; j<i; j++) {
            Y_ij[NumVortices*i + j] = pos_Real[2*i + 1] - pos_Real[2*j + 1];
            Y_ij[NumVortices*j + i] = Y_ij[NumVortices*i + j];
        }
    }

}


// this function calculates r_{ij}^2 and stores it in an array called R_ij_squared
void r_ij_squared(double* X_ij, double* Y_ij, double* R_ij_squared, double* pos_Real, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {
        for (int j=0; j<i; j++) {
            R_ij_squared[NumVortices*i + j] = X_ij[NumVortices*i + j]*X_ij[NumVortices*i + j] + Y_ij[NumVortices*i + j]*Y_ij[NumVortices*i + j];
            R_ij_squared[NumVortices*j + i] = R_ij_squared[NumVortices*i + j];
        }
    }

}


// this function calculates and returns an array that contains the position of each image charge based on the locations of the real charges inside the circular domain of radius R
void posImage(double* pos_Image, double* pos_Real, double R, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {
        for (int j=0; j<2; j++) {
            pos_Image[2*i + j] = R*R*pos_Real[2*i + j]/(pos_Real[2*i + 0]*pos_Real[2*i + 0] + pos_Real[2*i + 1]*pos_Real[2*i + 1]);
        }
    }

}


// this function calculates x_bar_{ij}
void x_bar_ij(double* X_bar_ij, double* pos_Real, double* pos_Image, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {

        X_bar_ij[NumVortices*i + i] = pos_Real[2*i + 0] - pos_Image[2*i + 0]; // sets diagonal element equal to x_i - \bar{x}_i

        for (int j=0; j<NumVortices; j++) {
            X_bar_ij[NumVortices*i + j] = pos_Real[2*i + 0] - pos_Image[2*j + 0];
        }
    }

}


// this function calculates y_bar_{ij}
void y_bar_ij(double* Y_bar_ij, double* pos_Real, double* pos_Image, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {

        Y_bar_ij[NumVortices*i + i] = pos_Real[2*i + 0] - pos_Image[2*i + 0]; // sets diagonal element equal to x_i - \bar{x}_i

        for (int j=0; j<NumVortices; j++) {
            Y_bar_ij[NumVortices*i + j] = pos_Real[2*i + 0] - pos_Image[2*j + 0];
        }
    }

}


// this function calculates \bar{r}_{ij}^2
void r_bar_ij_squared(double* R_bar_ij_squared, double* pos_Real, double* pos_Image, double* X_bar_ij, double* Y_bar_ij, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {
        for (int j=0; j<NumVortices; j++) {
            R_bar_ij_squared[NumVortices*i + j] = X_bar_ij[NumVortices*i + j]*X_bar_ij[NumVortices*i + j] + Y_bar_ij[NumVortices*i + j]*Y_bar_ij[NumVortices*i + j];
        }
    }

}


// ----------------------------------------------- RUNGE-KUTTA ----------------------------------------------- //

// this function calculates v_x of each point vortex as a vector
void v_x(double* Vx, double* X_ij, double* Y_ij, double* R_ij_squared, double* X_bar_ij, double* Y_bar_ij, double* R_bar_ij_squared, double* pos_Real, double* pos_Image, double* gamma, int NumVortices) {

    // this for loop completes the sum
    for (int i=0; i<NumVortices; i++) { // outer loop goes through each particle according to the index i

        Vx[i] = 0; // clears the previous value of Vx that was in the ith spot from the previous integration step of the simulation
        
        for (int j=0; j<NumVortices; j++) { // inner summation loop goes through each term in the sum according to the index j
            
            if (j==i) {
                Vx[i] = Vx[i] + (1.0/2*M_PI)*gamma[j]*Y_bar_ij[NumVortices*i + j]/R_bar_ij_squared[NumVortices*i + j];
            } else {
                Vx[i] = Vx[i] - (1.0/2*M_PI)*gamma[j]*Y_ij[NumVortices*i + j]/(R_ij_squared[NumVortices*i + j]) + (1.0/2*M_PI)*gamma[j]*Y_bar_ij[NumVortices*i + j]/R_bar_ij_squared[NumVortices*i + j];
            }

        }

    }

}



// this function calculates v_y of each point vortex as a vector
void v_y(double* Vy, double* X_ij, double* Y_ij, double* R_ij_squared, double* X_bar_ij, double* Y_bar_ij, double* R_bar_ij_squared, double* pos_Real, double* pos_Image, double* gamma, int NumVortices) {

    for (int i=0; i<NumVortices; i++) {

        Vy[i] = 0; // clears the previous value of Vy that was in the ith spot from the previous integration step of the simulation

        for (int j=0; j<NumVortices; j++) {

            if (j==i) {
                Vy[i] = Vy[i] - (1.0/2*M_PI)*gamma[j]*X_bar_ij[NumVortices*i + j]/R_bar_ij_squared[NumVortices*i + j];
            }
            else {
                Vy[i] = Vy[i] + (1.0/2*M_PI)*gamma[j]*X_ij[NumVortices*i + j]/R_ij_squared[NumVortices*i + j] - (1.0/2*M_PI)*gamma[j]*X_bar_ij[NumVortices*i + j]/R_bar_ij_squared[NumVortices*i + j];
            }

        }

    }

}


// calculates the correction to the x- and y-coordinates of each point vortex and stores them each in a vector and puts them in a tuple according to the Runge-Kutta method
std::tuple<double*, double*> RK(double* Vx, double* Vy, double* X_ij, double* Y_ij, double* R_ij_squared, double* X_bar_ij, double* Y_bar_ij, double* R_bar_ij_squared, double* pos_Real, double* pos_Image, double* gamma, int NumVortices, double dt) {
    
    // calculates vx and vy and stores the results in Vx and Vy, respectively
    v_x(Vx, X_ij, Y_ij, R_ij_squared, X_bar_ij, Y_bar_ij, R_bar_ij_squared, pos_Real, pos_Image, gamma, NumVortices);
    v_y(Vy, X_ij, Y_ij, R_ij_squared, X_bar_ij, Y_bar_ij, R_bar_ij_squared, pos_Real, pos_Image, gamma, NumVortices);

    // calculates the position vector at the current step, plus 0.5*dt*f0, to be used in Runge-Kutta to calculate f1
    double* vec1_Real = pos_Real;
    double* vec1_Image = pos_Image;

    for (int i=0; i<NumVortices; i++) {
        vec1_Real[2*i + 0] = vec1_Real[2*i + 0] + gamma[i]*pos_Real[2*i + 1];
        vec1_Real[2*i + 1] = vec1_Real[2*i + 1] - gamma[i]*pos_Real[2*i + 0];

        vec1_Real[2*i + 0] = vec1_Real[2*i + 0] + 0.5*dt*Vx[i];
        vec1_Image[2*i + 0] = vec1_Image[2*i + 0] + 0.5*dt*Vx[i];

        vec1_Real[2*i + 1] = vec1_Real[2*i + 1] +  0.5*dt*Vx[i];
        vec1_Image[2*i + 1] = vec1_Image[2*i + 1] +  0.5*dt*Vy[i];
    }

    // AVX: 
    // make a new array of type __m256 (SIMD vectors)
    // _mm256_loadu_pd all the data from vec of doubles to SIMD vector
    // _mm256_add_pd whatever numbern in SIMD vector form to the vector 
    // _mm256_add_pd(some number, your vector)

    double* f1_x[NumVortices];
    double* f1_y[NumVortices];
    v_x(vec1_Real, vec1_Image, gamma);
    v_y(vec1_Real, vec1_Image, gamma);

    // calculates the position vector at the current step plus 0.5*dt*f1, to be used in Runge-Kutta
    double* vec2_Real = pos_Real;
    double* vec2_Image = pos_Image;

    for (int i=0; i<NumVortices; i++) {


        vec2_Real[i][0] = vec2_Real[i][0] + 0.5*dt*f1_x[i];
        vec2_Image[i][0] = vec2_Image[i][0] + 0.5*dt*f1_x[i];

        vec2_Real[i][1] = vec2_Real[i][1] + 0.5*dt*f1_y[i];
        vec2_Image[i][1] = vec2_Image[i][1] + 0.5*dt*f1_y[i];
    }

    std::vector<double> f2_x = xDeriv(vec2_Real, vec2_Image, gamma);
    std::vector<double> f2_y = yDeriv(vec2_Real, vec2_Image, gamma);

    // calculates the position vector at the current step plus 0.5*dt*f1, to be used in Runge-Kutta
    std::vector<std::vector<double>> vec3_Real = vec2_Real;
    std::vector<std::vector<double>> vec3_Image = vec2_Image;

    for (int i=0; i<pos_Real.size(); i++) {
        vec3_Real[i][0] = vec3_Real[i][0] + 0.5*dt*f2_x[i];
        vec3_Image[i][0] = vec3_Image[i][0] + 0.5*dt*f2_x[i];
        
        vec3_Real[i][1] = vec3_Real[i][1] + 0.5*dt*f2_y[i];
        vec3_Image[i][1] = vec3_Image[i][1] + 0.5*dt*f2_y[i];
    }

    std::vector<double> f3_x = xDeriv(vec3_Real, vec3_Image, gamma);
    std::vector<double> f3_y = yDeriv(vec3_Real, vec3_Image, gamma);

    std::vector<double> vec_x(pos_Real.size());
    std::vector<double> vec_y(pos_Real.size());

    for (int i=0; i<pos_Real.size(); i++) {
        vec_x[i] = dt*(f0_x[i] + 2*f1_x[i] + 2*f2_x[i] + f3_x[i])/6;
        vec_y[i] = dt*(f0_y[i] + 2*f1_y[i] + 2*f2_y[i] + f3_y[i])/6;
    }

    return std::make_tuple(vec_x, vec_y);

}


// --------------------------------------------------------- END OF FUNCTIONS --------------------------------------------------------- //


int main() {

    // ------------------------------------ DEFINING CONSTANTS ------------------------------------ //

    int N = 10000; // number of time steps
    double tmax = 10; // simulation time
    double dt = tmax/N; // time step

    double R = 10; // radius of superfluid


    // ------------------------------------ INITIAL CONDITION ------------------------------------ //

    std::vector<std::vector<double>> vortexPosReal = { {-1, 0} , {1, 0} };
    std::vector<std::vector<double>> vortexPosImage = posImage(vortexPosReal, R);

    std::vector<double> gamma = {-8, 8}; // vector of circulations of the point vortices
    double damping = 0.01; // damping coefficient
    
    std::ofstream file("test.dat");

    file << vortexPosReal[0][0] << " " << vortexPosReal[0][1] << " " << vortexPosReal[1][0] << " " << vortexPosReal[1][1] << std::endl;


    // ------------------------------------ SIMULATION LOOP ------------------------------------ //

    

    for (int i=1; i<N; i++) {

        std::vector<double> x_correction = std::get<0>(RK(vortexPosReal, vortexPosImage, gamma, dt)); 
        std::vector<double> y_correction = std::get<1>(RK(vortexPosReal, vortexPosImage, gamma, dt));
        
        for (int j=0; j<vortexPosReal.size(); j++) {
            vortexPosReal[j][0] = vortexPosReal[j][0] + x_correction[j];
            vortexPosReal[j][1] = vortexPosReal[j][1] + y_correction[j];

            file << vortexPosReal[j][0] << " " << vortexPosReal[j][1] << " ";

        }

        file << std::endl;

        vortexPosImage = posImage(vortexPosReal, R);
        
    }

   

    return 0;
}