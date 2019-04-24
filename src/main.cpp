/*
 * main.cpp
 *
 *  Created on: Apr 2, 2019
 *      Author: d-w-h
 *
 *      Implementation of the LBM.
 *      Simulation of Couette flow.
 */

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <vector>
#include <string>

#include "../inc/memory_operations.hpp"
#include "../inc/user.hpp"

int main(int argc, char* argv[])
{
    /*Parameter declaration*/
    int Nx, Ny, Nd;
    int max_timesteps, timestep;
    double dt, tau, dt_tau, dx;
    double u_lid;
    std::map<int, int> index_map;

    Nx = 20;
    Ny = 10;
    Nd = 9;

    max_timesteps = 500;
    dt = 1.0;
    dx = 1.0;
    tau = 0.8;
    u_lid = 0.1;

    dt_tau = dt/tau;

    /*Memory allocation*/
    double*** fi     = allocate_discrete_velocity_2D(Nx, Ny, Nd);
    double*** fieq   = allocate_discrete_velocity_2D(Nx, Ny, Nd);
    double*** fiprop = allocate_discrete_velocity_2D(Nx, Ny, Nd);
    double**     rho = allocate_2D_array_doubles(Nx, Ny);
    double**      Ux = allocate_2D_array_doubles(Nx, Ny);
    double**      Uy = allocate_2D_array_doubles(Nx, Ny);
    double** vort_field = allocate_2D_array_doubles(Nx, Ny);
    double**       X = allocate_2D_array_doubles(Nx, Ny);
    double**       Y = allocate_2D_array_doubles(Nx, Ny);
    double*       wi = allocate_1D_array_doubles(Nx);
    double*      cix = allocate_1D_array_doubles(Nx);
    double*      ciy = allocate_1D_array_doubles(Nx);

    /*Initialization*/
    /*Initialize index map*/
    index_map.insert(std::pair <int, int> (0, 0));
    index_map.insert(std::pair <int, int> (1, 3));
    index_map.insert(std::pair <int, int> (2, 4));
    index_map.insert(std::pair <int, int> (3, 1));
    index_map.insert(std::pair <int, int> (4, 2));
    index_map.insert(std::pair <int, int> (5, 7));
    index_map.insert(std::pair <int, int> (6, 8));
    index_map.insert(std::pair <int, int> (7, 5));
    index_map.insert(std::pair <int, int> (8, 6));

    /*Initializing cij and wi*/
    for(int i = 0; i < 27; ++i){
        cix[0] = 0; ciy[0] = 0; wi[0] = 4.0/9.0;
        cix[1] = 1; ciy[1] = 0; wi[1] = 1.0/9.0;
        cix[2] = 0; ciy[2] = 1; wi[2] = 1.0/9.0;
        cix[3] = -1;ciy[3] = 0; wi[3] = 1.0/9.0;
        cix[4] = 0; ciy[4] = -1;wi[4] = 1.0/9.0;
        cix[5] = 1; ciy[5] = 1; wi[5] = 1.0/36.0;
        cix[6] = -1;ciy[6] = 1; wi[6] = 1.0/36.0;
        cix[7] = -1;ciy[7] = -1;wi[7] = 1.0/36.0;
        cix[8] = 1; ciy[8] = -1;wi[8] = 1.0/36.0;
    }

    /*Initializing density and velocity*/
    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j) {
            Ux[i][j] = 0.0;
            Uy[i][j] = 0.0;
            vort_field[i][j] = 0.0;
            X[i][j] = i*dx + 0.5*dx;
            Y[i][j] = j*dx + 0.5*dx;
            rho[i][j] = 1000.0;
        }

    /*Initializing discrete velocity distributions*/
    for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
            for(int d = 0; d < Nd; ++d) {
                fi[i][j][d] = wi[d]*rho[i][j];
                fieq[i][j][d] = wi[d]*rho[i][j];
                fiprop[i][j][d] = wi[d]*rho[i][j];
            }


    /*LBM algorithm*/
    timestep = 0;
    double max_v_y = 0.0;
    double min_v_y = 1000.0;
    double max_v_x = 0.0;
    double min_v_x = 1000.0;
    double max_speed = 0.0;
    double max_vorticity = 0.0;
    do
    {
        /*Calculate density and velocity*/
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j) {
                rho[i][j] = 0.0;
                Ux[i][j] = 0.0;
                Uy[i][j] = 0.0;
                for(int d = 0; d < Nd; ++d)
                {
                    rho[i][j] += fiprop[i][j][d];
                    Ux[i][j] += cix[d] * fiprop[i][j][d];
                    Uy[i][j] += ciy[d] * fiprop[i][j][d];
                }
                Ux[i][j] = Ux[i][j]/rho[i][j];
                Uy[i][j] = Uy[i][j]/rho[i][j];
            }

        /*Calculate vorticity*/
        for(int i = 1; i < Nx - 1; ++i)
            for(int j = 1; j < Ny - 1; ++j) {
                vort_field[i][j] = 0.5 * (Uy[i+1][j] - Uy[i-1][j] - Ux[i][j+1] + Ux[i][j-1]) * 1.0 / dx;
            }

        /*Calculate equilibrium distribution*/
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j) {
                for(int d = 0; d < Nd; ++d) {
                    double uci = cix[d]*Ux[i][j] + ciy[d]*Uy[i][j];
                    double uu = Ux[i][j]*Ux[i][j] + Uy[i][j]*Uy[i][j];
                    fieq[i][j][d] = wi[d]*rho[i][j]*(1.0 + 3.0*uci + 4.5*uci*uci - 1.5*uu);
                }
            }

        /*Collision*/
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j) {
                for(int d = 0; d < Nd; ++d) {
                    fi[i][j][d] = fiprop[i][j][d]*(1 - dt_tau) + fieq[i][j][d]*dt_tau;
                }
            }

        /*Streaming internal nodes*/
        for(int i = 1; i < Nx - 1; ++i)
            for(int j = 1; j < Ny - 1; ++j) {
                for(int d = 0; d < Nd; ++d) {
                    int ir, jr;
                    ir = i - cix[d];
                    jr = j - ciy[d];
                    fiprop[i][j][d] = fi[ir][jr][d];
                }
            }

        /*Streaming south side*/
        for(int i = 0; i < Nx; ++i)
            for(int d = 0; d < Nd; ++d) {
                int ir, jr;
                ir = i - cix[d];
                jr = 0 - ciy[d];
                if(ciy[d] == 1) {
                    fiprop[i][0][d] = fi[i][0][index_map[d]];
                }
                if(ir == -1)
                    ir = Nx - 1;
                if(ir == Nx)
                    ir = 0;
                if(ciy[d] != 1) {
                    fiprop[i][0][d] = fi[ir][jr][d];
                }
            }

        /*Streaming north side*/
        for(int i = 0; i < Nx; ++i)
            for(int d = 0; d < Nd; ++d) {
                int ir, jr;
                ir = i - cix[d];
                jr = Ny - 1 - ciy[d];
                if(ciy[d] == -1) {
                    fiprop[i][Ny-1][d] = fi[i][Ny-1][index_map[d]] - 2*wi[index_map[d]]*rho[i][Ny-1]*3.0*(u_lid*cix[index_map[d]]);
                }
                if(ir == -1)
                    ir = Nx - 1;
                if(ir == Nx)
                    ir = 0;
                if(ciy[d] != -1) {
                    fiprop[i][Ny-1][d] = fi[ir][jr][d];
                }
            }

        /*Streaming west side*/
        for(int j = 1; j < Ny - 1; ++j)
            for(int d = 0; d < Nd; ++d) {
                int ir, jr;
                ir = 0 - cix[d];
                jr = j - ciy[d];
                if(ir == -1)
                    ir = Nx - 1;
                fiprop[0][j][d] = fi[ir][jr][d];
            }

        /*Streaming east side*/
        for(int j = 1; j < Ny - 1; ++j)
            for(int d = 0; d < Nd; ++d) {
                int ir, jr;
                ir = Nx - 1 - cix[d];
                jr = j - ciy[d];
                if(ir == Nx)
                    ir = 0;
                fiprop[Nx-1][j][d] = fi[ir][jr][d];
            }

        /*Calculate max and min velocity, speed and vorticity*/
        double dummy_vel = 0.0;
        double dummy_vort = 0.0;
        for(int i = 0; i < Nx; ++i)
            for(int j = 0; j < Ny; ++j) {
                dummy_vel  = Ux[i][j] * Ux[i][j] + Uy[i][j] * Uy[i][j];
                dummy_vort = vort_field[i][j] * vort_field[i][j];

                if(Uy[i][j] > max_v_y)
                    max_v_y = Uy[i][j];
                if(Uy[i][j] < min_v_y)
                    min_v_y = Uy[i][j];

                if(Ux[i][j] > max_v_x)
                    max_v_x = Ux[i][j];
                if(Ux[i][j] < min_v_x)
                    min_v_x = Ux[i][j];

                if(dummy_vel > max_speed)
                    max_speed = dummy_vel;
                if(dummy_vort > max_vorticity)
                    max_vorticity = dummy_vort;
            }
        max_speed  = sqrt(max_speed);
        max_vorticity = sqrt(max_vorticity);

        /*Write time dependent data*/
        FILE* data = NULL;
        char file_name_complete[50];
        char* file_name_suffix = (char*) ".txt";
        char* file_name = (char*) "./simdata/data";

        sprintf(file_name_complete, "%s_%d%s", file_name, timestep, file_name_suffix);
        data = fopen(file_name_complete, "w");
        if(data != NULL) {
            for(int j = Ny - 1;j >= 0 ; --j)
                for(int i = 0;i < Nx; ++i) {
                    fprintf(data,"%f    ", i*dx + 0.5*dx);
                    fprintf(data,"%f    ", j*dx + 0.5*dx);
                    fprintf(data,"%f    ", Ux[i][j]);
                    fprintf(data,"%f    ", Uy[i][j]);
                    fprintf(data,"%f\n", vort_field[i][j]);
                }
        }
        fclose(data);

        std::cout << "timestep: " << timestep << std::endl;

        timestep++;

    } while(timestep < max_timesteps);
    /*End LBM algorithm*/

    /*Write time independent data*/
    FILE* params = NULL;
    params = fopen("parameters.txt", "w");
    if(params != NULL) {
        fprintf(params, "%i\n", Nx);
        fprintf(params, "%i\n", Ny);
        fprintf(params, "%i\n", timestep);
        fprintf(params,"%f\n", max_speed * 1.15);
        fprintf(params,"%f\n", max_vorticity * 1.15);
        fprintf(params,"%f\n", min_v_y);
        fprintf(params,"%f\n", max_v_y);
        fprintf(params,"%f\n", min_v_x);
        fprintf(params,"%f\n", max_v_x);
    }
    fclose(params);

    std::cout << "done" << std::endl;

    return 0;
}
