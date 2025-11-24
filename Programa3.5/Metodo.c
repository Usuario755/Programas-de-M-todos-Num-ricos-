#include "Metodo.h"
//Cuerpos de fuciones
void mostrarMenu() {
    printf("\n Menu:\n");
    printf("1. Metodo de la Secante funcion 1\n");
    printf("2. Metodo de la Secante funcion 2\n");
    printf("3. Gauss-Seidel\n");
    printf("4. Cholesky\n");
    printf("5. Valores Propios\n");
    printf("0. Salir\n");
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
    printf("No se alcanzÃ³ la tolerancia tras %d iteraciones\n", max_iter);
}

static double f(double x) {
    return x * x - 4; // Ejemplo: raÃ­z en x=2
}
// FunciÃ³n 1: x^2 cos(x) - 2x

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
    printf("La raíz es %lf\n", x2);
}

// FunciÃ³n 2: (6 - 2/x^3)(e^2 + x/4) + 1
void fn2(double x0, double x1, int kmax, double tol) {
    printf("\nMÃ©todo de la Secante\n");
    printf("Iter\tx0\tx1\tx2\tf(x2)\terror\n");
    for (int k = 1; k <= kmax; k++) {
        double fx0 = f(x0);
        double fx1 = f(x1);
        if (fabs(fx1 - fx0) < 1e-12) {
            printf("DivisiÃ³n por cero\n");
            return;
        }
        double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        double error = fabs(x2 - x1);
        printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", k, x0, x1, x2, f(x2), error);
        if (error < tol) {
            printf("Convergencia en %d iteraciones. RaÃ­z: %.6f\n", k, x2);
            return;
        }
        x0 = x1;
        x1 = x2;
    }
    printf("No se alcanzÃ³ la tolerancia tras %d iteraciones\n", kmax);
}

void leerMatriz(int n, double Matrix[MAX][MAX], double b[MAX]) {
    printf("Ingrese la matriz (%dx%d):\n", n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("[%d][%d]:",i,j);
            scanf("%lf", &Matrix[i][j]);
        }
        printf(" \n");
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
void leerVector(int n, double b[MAX]) {
    printf("Ingrese el vector b:\n");
    for (int i = 0; i < n; i++) {
        printf("[%d]:",i);
        scanf("%lf", &b[i]);
    }
}


void imprimirVector(int n, double v[MAX]) {
    for (int i = 0; i < n; i++) {
        printf("%8.6f\n", v[i]);
    }
}
//Cholesky
void programa3() {
    int n;
    printf("Introduzca el tamaño de la matriz: ");
    scanf("%d", &n);

    double A[MAX][MAX], L[MAX][MAX], U[MAX][MAX], b[MAX], y[MAX], x[MAX];

    // Inicializar L y U
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;
        }
    }

    // Leer matriz A
    printf("Ingrese la matriz A (%dx%d):\n", n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("[%d][%d]:",i,j);
            scanf("%lf", &A[i][j]);
        }
    }

    // Leer vector b
    printf("Ingrese el vector b:\n");
    for (int i = 0; i < n; i++) {
        printf("[%d]:",i);
        scanf("%lf", &b[i]);
    }

    // Descomposición LU
    for (int j = 0; j < n; j++) {
        // Calcular elementos de L
        for (int i = j; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - sum;
        }

        // Calcular elementos de U
        for (int k = j + 1; k < n; k++) {
            double sum = 0.0;
            for (int s = 0; s < j; s++) {
                sum += L[j][s] * U[s][k];
            }
            U[j][k] = (A[j][k] - sum) / L[j][j];
        }
        U[j][j] = 1.0; // Diagonal de U
    }

    // Sustitución hacia adelante (Ly = b)
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // Sustitución hacia atrás (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = y[i] - sum; // U tiene 1 en diagonal
    }

    // Mostrar resultados
    printf("\nSoluciones del sistema:\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.6lf\n", i, x[i]);
    }
}
//Valores propios 
double normaInfinita(int n, double v[MAX]) {
    double max = 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(v[i]) > max) {
            max = fabs(v[i]);
        }
    }
    return max;
}

void metodoPotencia(int n, double A[MAX][MAX], double *lambda_max, double v[MAX], double tol, int max_iter) {
    double v_old[MAX];
    double z[MAX]; 
    double error = 1.0;
    int iter = 0;

    for (int i = 0; i < n; i++) v[i] = 1.0;

    printf("\n--- Metodo de la Potencia (MAXIMO) ---\n");
    printf("Iter\tLambda Max\tError\n");

    while (error > tol && iter < max_iter) {
        for(int i = 0; i < n; i++) v_old[i] = v[i];
        
        for (int i = 0; i < n; i++) {
            z[i] = 0.0;
            for (int j = 0; j < n; j++) {
                z[i] += A[i][j] * v[j];
            }
        }

        double mu_max = 0.0;
        double max_abs = 0.0;
        
        for(int i=0; i<n; i++) {
            if(fabs(z[i]) > max_abs) {
                max_abs = fabs(z[i]);
                mu_max = z[i];
            }
        }
        
        *lambda_max = mu_max; 

        if (fabs(*lambda_max) < 1e-9) { printf("Lambda cercano a 0.\n"); break; }

        for (int i = 0; i < n; i++) {
            v[i] = z[i] / *lambda_max;
        }

        error = 0.0;
        for (int i = 0; i < n; i++) {
            error += fabs(v[i] - v_old[i]);
        }

        iter++;
        printf("%d\t%.6f\t%.6f\n", iter, *lambda_max, error);
    }

    if (iter >= max_iter) {
        printf("Advertencia: Maximo de iteraciones alcanzado.\n");
    } else {
        printf("Convergencia alcanzada.\n");
    }
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
void sustAtras(int n, double L[MAX][MAX], double y[MAX], double x[MAX]) {
    for (int i = n - 1; i >= 0; i--) {
        double suma = 0.0;
        for (int j = i + 1; j < n; j++) {
            suma += L[j][i] * x[j];
        }
        x[i] = (y[i] - suma) / L[i][i];
    }
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
void metodoPotenciaInversa(int n, double A[MAX][MAX], double *lambda_min, double v_min[MAX], double tol, int max_iter) {
    double v_old[MAX];
    double z[MAX]; 
    double y[MAX]; 
    double L[MAX][MAX]; 
    double error = 1.0;
    int iter = 0;
    
    if (!cholesky(n, A, L)) {
        printf("ERROR (Potencia Inversa): La matriz no es S.D.P. No se pudo factorizar.\n");
        *lambda_min = 0.0;
        return;
    }

    for (int i = 0; i < n; i++) v_min[i] = 1.0;

    printf("\n--- Metodo de la Potencia Inversa (MINIMO) ---\n");
    printf("Iter\t1/Lambda (mu)\tLambda Min\tError\n");

    while (error > tol && iter < max_iter) {
        for(int i = 0; i < n; i++) v_old[i] = v_min[i];
        
        sustDelante(n, L, v_min, y); 
        sustAtras(n, L, y, z);
        
        double mu = 0.0; 
        double max_abs = 0.0;
        
        for(int i=0; i<n; i++) {
            if(fabs(z[i]) > max_abs) {
                max_abs = fabs(z[i]);
                mu = z[i];
            }
        }
        
        if (fabs(mu) < 1e-9) { printf("Mu cercano a 0.\n"); break; }

        for (int i = 0; i < n; i++) {
            v_min[i] = z[i] / mu;
        }

        error = 0.0;
        for (int i = 0; i < n; i++) {
            error += fabs(v_min[i] - v_old[i]);
        }

        *lambda_min = 1.0 / mu; 
        
        iter++;
        printf("%d\t%.6f\t\t%.6f\t%.6f\n", iter, mu, *lambda_min, error);
    }
    
    if (iter >= max_iter) {
        printf("Advertencia: Maximo de iteraciones alcanzado.\n");
    } else {
        printf("Convergencia alcanzada.\n");
    }
}

void calcularValoresPropios(int n, double A[MAX][MAX], double tol, int max_iter) {
    double v[MAX], v_min[MAX], lambda_max, lambda_min;

    // Inicializar vectores
    for (int i = 0; i < n; i++) {
        v[i] = 1.0;
        v_min[i] = 1.0;
    }

    // Método de Potencia
    metodoPotencia(n, A, &lambda_max, v, tol, max_iter);
    printf("\nValor propio máximo aproximado: %.6lf\n", lambda_max);
    printf("Vector propio asociado:\n");
    for (int i = 0; i < n; i++) {
        printf("%.6lf ", v[i]);
    }
    printf("\n");

    // Método de Potencia Inversa
    metodoPotenciaInversa(n, A, &lambda_min, v_min, tol, max_iter);
    printf("\nValor propio mínimo aproximado: %.6lf\n", lambda_min);
    printf("Vector propio asociado:\n");
    for (int i = 0; i < n; i++) {
        printf("%.6lf ", v_min[i]);
    }
    printf("\n");
}
