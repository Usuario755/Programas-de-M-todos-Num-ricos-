//Y esto es el código de las funciones
//Macros
#ifndef METODOS_H
#define METODOS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MAX 100

//Prototipos
void mostrarMenu();
void leerMatriz(int n, double Matrix[MAX][MAX], double b[MAX]);
void imprimirMatriz(int n, double Matrix[MAX][MAX]);
void imprimirVector(int n, double v[MAX]);
void gaussSeidel(int n, double Matrix[MAX][MAX], double b[MAX], double X[MAX], double tolerancia, int max_iter);

void fn1(double x0, double x1, int kmax, double tol);

void fn2(double x0, double x1, int kmax, double tol);
int cholesky(int n, double A[MAX][MAX], double L[MAX][MAX]);


#endif
