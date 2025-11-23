#ifndef METODOS_H
#define METODOS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX 100
void mostrarMenu() {
    printf("\nMenu:\n");
    printf("1. Método de la Secante\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Cholesky\n");
    printf("0. Salir\n");
}

void leerMatriz(int n, double Matrix[MAX][MAX], double b[MAX]) {
    printf("Ingrese la matriz (%dx%d):\n", n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &Matrix[i][j]);
        }
    }
    printf("Ingrese el vector b: ");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }
}

void imprimirMatriz(int n, double Matrix[MAX][MAX]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%8.3f ", Matrix[i][j]);
        }
        printf("\n");
    }
}

void imprimirVector(int n, double v[MAX]) {
    for (int i = 0; i < n; i++) {
        printf("%8.6f\n", v[i]);
    }
}

int cholesky(int n, double Matrix[MAX][MAX], double L[MAX][MAX]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }
            if (i == j) {
                double val = Matrix[i][i] - sum;
                if (val <= 0) {
                    printf("ERROR: La matriz no es definida positiva.\n");
                    return 1;
                }
                L[i][j] = sqrt(val);
            } else {
                L[i][j] = (Matrix[i][j] - sum) / L[j][j];
            }
        }
    }
    return 0;
}

void gaussSeidel(int n, double Matrix[MAX][MAX], double b[MAX], double X[MAX], double tolerancia, int max_iter) {
    for (int i = 0; i < n; i++) X[i] = 0.0;
    for (int iter = 1; iter <= max_iter; iter++) {
        double error = 0.0;
        for (int i = 0; i < n; i++) {
            double suma = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) suma -= Matrix[i][j] * X[j];
            }
            double nuevo = suma / Matrix[i][i];
            error += fabs(nuevo - X[i]);
            X[i] = nuevo;
        }
        if (error < tolerancia) {
            printf("Convergencia alcanzada en %d iteraciones\n", iter);
            return;
        }
    }
    printf("No se alcanzó la tolerancia tras %d iteraciones\n", max_iter);
}

static double f(double x) {
    return x * x - 4; // Ejemplo: raíz en x=2
}

void fn2(double x0, double x1, int kmax, double tol) {
    printf("\nMétodo de la Secante\n");
    printf("Iter\tx0\tx1\tx2\tf(x2)\terror\n");
    for (int k = 1; k <= kmax; k++) {
        double fx0 = f(x0);
        double fx1 = f(x1);
        if (fabs(fx1 - fx0) < 1e-12) {
            printf("División por cero\n");
            return;
        }
        double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        double error = fabs(x2 - x1);
        printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", k, x0, x1, x2, f(x2), error);
        if (error < tol) {
            printf("Convergencia en %d iteraciones. Raíz: %.6f\n", k, x2);
            return;
        }
        x0 = x1;
        x1 = x2;
    }
    printf("No se alcanzó la tolerancia tras %d iteraciones\n", kmax);
}

#endif
