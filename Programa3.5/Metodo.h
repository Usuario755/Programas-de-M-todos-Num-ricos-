//Y esto es el código de las funciones
//Macros
#ifndef METODO_H
#define METODO_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MAX 100

//Prototipos
void mostrarMenu();
void sustAtras(int n, double L[MAX][MAX], double y[MAX], double x[MAX]);
void leerMatriz(int n, double Matrix[MAX][MAX], double b[MAX]);
void imprimirMatriz(int n, double Matrix[MAX][MAX]);
void leerVector(int n, double b[MAX]);
void sustDelante(int n, double L[MAX][MAX], double b[MAX], double Z[MAX]);
void sustregres(int n, double L[MAX][MAX], double y[MAX], double x[MAX]);
void imprimirVector(int n, double v[MAX]);
void gaussSeidel(int n, double Matrix[MAX][MAX], double b[MAX], double X[MAX], double tolerancia, int max_iter);

void fn1(double x0, double x1, int kmax, double tol);
//cholesky
void programa3();
void fn2(double x0, double x1, int kmax, double tol);
int cholesky(int n, double A[MAX][MAX], double L[MAX][MAX]);

//Prototipos para valores propios
void metodoPotencia(int n, double A[MAX][MAX], double *lambda, double v[MAX], double tol, int max_iter);
void metodoPotenciaInversa(int n, double A[MAX][MAX], double *lambda_min, double v_min[MAX], double tol, int max_iter);
void calcularValoresPropios(int n, double A[MAX][MAX], double tol, int max_iter);

#endif
