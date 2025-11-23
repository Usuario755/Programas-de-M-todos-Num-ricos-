#include <stdio.h>
#include "metodos.h"
#include <stdlib.h>
// Prototipos
int metodosecante();
int paquete2();
int programa3();

// Prototipos de funciones auxiliares
void fn1(double x0, double x1, int kmax, double tol);
void fn2(double x0, double x1, int kmax, double tol);
void fn3(double x0, double x1, int kmax, double tol);
void fn4(double x0, double x1, int kmax, double tol);
void readb();
void readdata();
void printmat();
int isEDD();
int determinante();
double calculateNorm(double *v1, double *v2);
void gaussSeidel(double tolerancia, int max_iter);
void eliminacionGaussiana();
void readdata3();
void readvector();
void CholeskyMat();
void sustdelante();
void sustregresiva();


int main() {
    int opcion;
    do {
        mostrarMenu();
        printf("Opción: ");
        scanf("%d", &opcion);
        if (opcion == 1) {
            fn2(1.0, 2.0, 10, 0.001);
        } else if (opcion == 2) {
            int n = 2;
            double Matrix[MAX][MAX] = {{4,1},{2,3}};
            double b[MAX] = {1,2};
            double X[MAX];
            gaussSeidel(n, Matrix, b, X, 0.001, 10);
            imprimirVector(n, X);
        } else if (opcion == 3) {
            int n = 2;
            double Matrix[MAX][MAX] = {{4,2},{2,3}};
            double L[MAX][MAX] = {0};
            if (cholesky(n, Matrix, L) == 0) imprimirMatriz(n, L);
        }
    } while(opcion != 0);
    return 0;
}
