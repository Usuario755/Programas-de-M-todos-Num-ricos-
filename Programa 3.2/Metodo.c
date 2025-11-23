#include "metodos.h"
//Cuerpos de fuciones
void mostrarMenu() {
    printf("\n Menu:\n");
    printf("1. Metodo de la Secante funcion 1\n");
    printf("1. Metodo de la Secante funcion 2\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Cholesky\n");
    printf("4. Eigenvalores\n");
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
    printf("No se alcanz贸 la tolerancia tras %d iteraciones\n", max_iter);
}

static double f(double x) {
    return x * x - 4; // Ejemplo: ra铆z en x=2
}
// Funci贸n 1: x^2 cos(x) - 2x
void fn1(double x0, double x1, int kmax, double tol) {
    double fx0, fx1, fx2, x2, Er;
    int k = 1;
    printf("\n |k|\t |x2|\t\t |fx2|\t\t |Er|\n ");
    do {
        fx0 = (x0 * x0) * cos(x0) - 2 * x0;
        fx1 = (x1 * x1) * cos(x1) - 2 * x1;
        x2 = (x1 - fx1) * ((x0 - x1) / (fx0 - fx1));
        fx2 = (x2 * x2) * cos(x2) - 2 * x2;
        if (k > 2)
            Er = fabs(x2 - x1) / fabs(x2);
        else
            Er = 0;
        printf("\n  |%d|\t |%lf|\t |%lf|\t |%lf|\n ", k, x2, fx2, Er);
        if (Er <= tol) {
            break;
        }
        x0 = x1;
        x1 = x2;
        k++;
    } while (k <= kmax);
    printf("La ra韟 es %lf\n", x2);
}

// Funci贸n 2: (6 - 2/x^3)(e^2 + x/4) + 1
void fn2(double x0, double x1, int kmax, double tol) {
    printf("\nM茅todo de la Secante\n");
    printf("Iter\tx0\tx1\tx2\tf(x2)\terror\n");
    for (int k = 1; k <= kmax; k++) {
        double fx0 = f(x0);
        double fx1 = f(x1);
        if (fabs(fx1 - fx0) < 1e-12) {
            printf("Divisi贸n por cero\n");
            return;
        }
        double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        double error = fabs(x2 - x1);
        printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", k, x0, x1, x2, f(x2), error);
        if (error < tol) {
            printf("Convergencia en %d iteraciones. Ra铆z: %.6f\n", k, x2);
            return;
        }
        x0 = x1;
        x1 = x2;
    }
    printf("No se alcanz贸 la tolerancia tras %d iteraciones\n", kmax);
}


int cholesky(int n, double A[MAX][MAX], double L[MAX][MAX]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++)
                sum += L[i][k] * L[j][k];

            if (i == j) {
                double val = A[i][i] - sum;
                if (val <= 0) {
                    printf("ERROR: La matriz no es definida positiva.\n");
                    return 0;
                }
                L[i][j] = sqrt(val);
            } else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    return 1;
}


void sustDelante(int n, double L[n][n], double b[n], double Z[n]) {
	for (int i = 0; i < n; i++) {
	double sum = 0.0;
	for (int j = 0; j < i; j++){
		sum += L[i][j] * Z[j];
		Z[i] = (b[i] - sum) / L[i][i];
		}	
	}
}

void sustregres(int n, double L[n][n], double y[n], double x[n]) {
	for (int i = n - 1; i >= 0; i--) {
	double sum = 0.0;
	for (int j = i + 1; j < n; j++){
		sum += L[j][i] * x[j];
		x[i] = (y[i] - sum) / L[i][i];
		}
	}
}

void programa3() {
	int n;
	printf("Introduzca el tama駉 de la matriz: ");
	scanf("%d", &n);
	double A[n][n], L[n][n], b[n], y[n], x[n];
	leerMatriz(n, A);
	leerVector(n, b);
	if (!cholesky(n, A, L)) return;
	sustDelante(n, L, b, y);
	sustAtras(n, L, y, x);


printf("\n Soluciones finales del sistema: \n");
	for (int i = 0; i < n; i++){
	printf("x[%d] = %.6lf", i, x[i]);
	}
}
