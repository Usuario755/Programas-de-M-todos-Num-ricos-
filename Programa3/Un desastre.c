#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#define MAX 100
using namespace std;

int main(){
	unsigned short int mmast; //Menu maestro
	char repeticion='s';
printf(" Paquete completo de programas solicitados a lo largo del curso de Metodos Numericos I\n ");
printf("A continuacion se mostrara un menu en donde puede elegir que programa revisar\n");
do{
printf("1. Metodo Abierto de la Secante\n");
printf("2. Resolucion de sistemas de ecuaciones Lineales (Gauss-Seidel)\n");
printf("3. Factorizacion LU de una matriz\n");
printf("4. Obtencion de Eigenvalues de una matriz\n");
scanf("%d",&mmast);
switch(mmast){
	case 1:
		printf("Metodo de la secante \n");
		// FunciÃ³n 1: x^2 cos(x) - 2x
void fn1(double x0, double x1, int kmax, double tol){
    double fx0, fx1, fx2, x2, Er;
    int k=1;
    printf("\n |k|\t |x2|\t\t |fx2|\t\t |Er|\n ");
    do{
        fx0 = (x0*x0)*cos(x0) - 2*x0;
        fx1 = (x1*x1)*cos(x1) - 2*x1;
        x2 = (x1 - fx1) * ((x0 - x1)/(fx0 - fx1));
        fx2 = (x2*x2)*cos(x2) - 2*x2;
        if(k > 2)
            Er = fabs(x2 - x1)/fabs(x2);
        else
            Er = 0;
        printf("\n  |%d|\t |%lf|\t |%lf|\t |%lf|\n ", k, x2, fx2, Er);
        if (Er <= tol){
            break;
        }
        x0 = x1;
        x1 = x2;
        k++;
    }while(k <= kmax);
    printf("La raÃ­z es %lf\n", x2);
}

// FunciÃ³n 2: (6 - 2/x^3)(e^2 + x/4) + 1
void fn2(double x0, double x1, int kmax, double tol){
    double fx0, fx1, fx2, x2, Er;
    int k=1;
    printf("\n |k|\t |x2|\t\t |fx2|\t\t |Er|\n ");
    do{
        fx0 = (6 - 2/pow(x0,3)) * (exp(2) + x0/4) + 1;
        fx1 = (6 - 2/pow(x1,3)) * (exp(2) + x1/4) + 1;
        x2 = (x1 - fx1) * ((x0 - x1)/(fx0 - fx1));
        fx2 = (6 - 2/pow(x2,3)) * (exp(2) + x2/4) + 1;
        if(k > 2)
            Er = fabs(x2 - x1)/fabs(x2);
        else
            Er = 0;
        printf("\n  |%d|\t |%lf|\t |%lf|\t |%lf|\n ", k, x2, fx2, Er);
        if (Er <= tol){
            break;
        }
        x0 = x1;
        x1 = x2;
        k++;
    }while(k <= kmax);
    printf("La raÃ­z es %lf\n", x2);
}

// FunciÃ³n 3: x^3 - 3sen(x^2) + 1
void fn3(double x0, double x1, int kmax, double tol){
    double fx0, fx1, fx2, x2, Er;
    int k=1;
    printf("\n |k|\t |x2|\t\t |fx2|\t\t |Er|\n ");
    do{
        fx0 = pow(x0,3) - 3*sin(pow(x0,2)) + 1;
        fx1 = pow(x1,3) - 3*sin(pow(x1,2)) + 1;
        x2 = (x1 - fx1) * ((x0 - x1)/(fx0 - fx1));
        fx2 = pow(x2,3) - 3*sin(pow(x2,2)) + 1;
        if(k > 2)
            Er = fabs(x2 - x1)/fabs(x2);
        else
            Er = 0;
        printf("\n  |%d|\t |%lf|\t |%lf|\t |%lf|\n ", k, x2, fx2, Er);
        if (Er <= tol){
            break;
        }
        x0 = x1;
        x1 = x2;
        k++;
    }while(k <= kmax);
    printf("La raÃ­z es %lf\n", x2);
}

// FunciÃ³n 4: x^3 + 6x^2 + 9.4x + 2.5
void fn4(double x0, double x1, int kmax, double tol){
    double fx0, fx1, fx2, x2, Er;
    int k=1;
    printf("\n |k|\t |x2|\t\t |fx2|\t\t |Er|\n ");
    do{
        fx0 = pow(x0,3) + 6*pow(x0,2) + 9.4*x0 + 2.5;
        fx1 = pow(x1,3) + 6*pow(x1,2) + 9.4*x1 + 2.5;
        x2 = (x1 - fx1) * ((x0 - x1)/(fx0 - fx1));
        fx2 = pow(x2,3) + 6*pow(x2,2) + 9.4*x2 + 2.5;
        if(k > 2)
            Er = fabs(x2 - x1)/fabs(x2);
        else
            Er = 0;
        printf("\n  |%d|\t |%lf|\t |%lf|\t |%lf|\n ", k, x2, fx2, Er);
        if (Er <= tol){
            break;
        }
        x0 = x1;
        x1 = x2;
        k++;
    }while(k <= kmax);
    printf("La raÃ­z es %lf\n", x2);
}

int metodosecante(){
    int kmax, menu;
    char rep;
    double x0, x1, tol;

    printf("Bienvenido, aquÃ­ la guÃ­a del programa\n");
    printf("Usted puede teclear la funciÃ³n que quiere analizar\n");
    printf("Presione enter para continuar...\n");
    cin.get();

    do{
        printf("\nTeclee la tolerancia de error deseada (ej. 0.0001):\n");
        scanf("%lf", &tol);
        printf("\nIntroduzca el nÃºmero mÃ¡ximo de iteraciones:\n");
        scanf("%d", &kmax);
        printf("\nElija una funciÃ³n para analizar:\n");
        printf("1. x^2 cos(x) - 2x\n");
        printf("2. (6 - 2/x^3)(e^2 + x/4) + 1\n");
        printf("3. x^3 - 3sen(x^2) + 1\n");
        printf("4. x^3 + 6x^2 + 9.4x + 2.5\n");
        scanf("%d", &menu);

        printf("\nIntroduzca el valor x0 para el intervalo inicial:\n");
        scanf("%lf", &x0);
        printf("\nIntroduzca el valor x1 para el intervalo inicial:\n");
        scanf("%lf", &x1);

        switch(menu){
            case 1:
                fn1(x0, x1, kmax, tol);
                break;
            case 2:
                fn2(x0, x1, kmax, tol);
                break;
            case 3:
                fn3(x0, x1, kmax, tol);
                break;
            case 4:
                fn4(x0, x1, kmax, tol);
                break;
            default:
                printf("\nOpciÃ³n invÃ¡lida, intente nuevamente.\n");
                break;
        }

        printf("\nÂ¿Desea repetir el programa? Teclee 's' para sÃ­:\n");
        scanf(" %c", &rep);
    }while(rep == 's');

    return 0;
}
		break;
	case 2:
		printf("Solucion de sistemas de ecuaciones lineales (Gauss-Seidel)\n");
		int n;
int Matrix[MAX][MAX];
double b[MAX]; 
double X[MAX]; 

void readb() {
    printf("\n--- Lectura del Vector de Terminos Independientes (b) ---\n");
    for (int i = 0; i < n; i++) {
        printf("Introduzca el valor de b[%d]: ", i);
        scanf("%lf", &b[i]); 
    }
}

void readdata() {
    printf("\n--- Lectura de la Matriz de Coeficientes (A) ---\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("Introduzca el valor de A[%d][%d]: ", i, j);
            scanf("%d", &Matrix[i][j]);
        }
    }
}

void printmat() {
    printf("\n--- Matriz Aumentada [A | b] ---\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("|%d| ", Matrix[i][j]);
        }
        printf("| %.2lf |\n", b[i]); 
    }
}


int isEDD(){
    printf("\n--- Verificacion EDD ---\n");
    int esEDD = 1;
    for (int i = 0; i < n; i++) {
        int suma = 0;
        for (int j = 0; j < n; j++) {
            if (i != j) suma += abs(Matrix[i][j]);
        }
        if (abs(Matrix[i][i]) <= suma) {
            esEDD = 0;
            break;
        }
    }
    if (esEDD){
        printf("\nLa matriz es Estrictamente Dominante Diagonalmente (EDD).");
    }
    else{
        printf("\nLa matriz NO es Estrictamente Dominante Diagonalmente.");
    }
    return esEDD;
}


int determinante(){
    if(n!=2){
        printf("\nEsta funcion solo calcula determinantes de matrices 2x2. (IMPLEMENTAR n x n)\n");
        return 1;
    }
    return Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];

}


double calculateNorm(double *v1, double *v2) {
    double max_error = 0.0;
    for (int i = 0; i < n; i++) {
        double error = fabs(v1[i] - v2[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}


void gaussSeidel(double tolerancia, int max_iter) {
    printf("\n--- Iniciando Metodo de Gauss-Seidel ---\n");
    printf("Tolerancia: %lf | Maximo de Iteraciones: %d\n", tolerancia, max_iter);

    double X_prev[MAX]; 
    int iter = 0;

    for (int i = 0; i < n; i++) {
        X_prev[i] = X[i];
    }

    while (iter < max_iter) {
        
        for (int i = 0; i < n; i++) {
            X_prev[i] = X[i];
        }

        for (int i = 0; i < n; i++) {
            double sum1 = 0.0; 
            double sum2 = 0.0; 

            for (int j = 0; j < n; j++) {
                if (i != j) {
                    if (j < i) {
                        sum1 += Matrix[i][j] * X[j]; 
                    } else {
                        sum2 += Matrix[i][j] * X_prev[j]; 
                    }
                }
            }

            if (Matrix[i][i] != 0) {
                 X[i] = (b[i] - sum1 - sum2) / Matrix[i][i];
            } else {
                printf("Error: Elemento diagonal Matrix[%d][%d] es cero. El metodo no puede continuar.\n", i, i);
                return;
            }
        }
        
        iter++;
        
        double error = calculateNorm(X, X_prev);

        printf("\nIteracion %d, Error: %.6lf\n", iter, error);
        printf("Vector X: (");
        for (int i = 0; i < n; i++) {
            printf("%.6lf%s", X[i], (i == n - 1) ? "" : ", ");
        }
        printf(")\n");
        
        if (error < tolerancia) {
            printf("\nCONVERGENCIA ALCANZADA.\n");
            break;
        }
        
    }

    if (iter >= max_iter) {
        printf("\nMAXIMO DE ITERACIONES ALCANZADO. El metodo termino sin converger al error deseado.\n");
    }

    printf("\n** SOLUCION FINAL OBTENIDA **\n");
    printf("La solucion X es: (");
    for (int i = 0; i < n; i++) {
        printf("%.6lf%s", X[i], (i == n - 1) ? "" : ", ");
    }
    printf(")\n");
}

/* ----------------------------------------------------------
   NUEVA FUNCION: Eliminacion Gaussiana con Sustitucion Regresiva
-------------------------------------------------------------*/
void eliminacionGaussiana() {
    printf("\n--- Iniciando Metodo de Eliminacion Gaussiana ---\n");

    double A[MAX][MAX];
    double B[MAX];
    double factor;

    // Copiar matriz y vector a variables locales tipo double
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = (double)Matrix[i][j];
        }
        B[i] = b[i];
    }

    // EliminaciÃ³n progresiva
    for (int k = 0; k < n - 1; k++) {
        if (fabs(A[k][k]) < 1e-12) {
            printf("Error: pivote nulo en la fila %d.\n", k);
            return;
        }
        for (int i = k + 1; i < n; i++) {
            factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            B[i] -= factor * B[k];
        }
    }

    // SustituciÃ³n regresiva
    for (int i = n - 1; i >= 0; i--) {
        double suma = 0.0;
        for (int j = i + 1; j < n; j++) {
            suma += A[i][j] * X[j];
        }
        if (fabs(A[i][i]) < 1e-12) {
            printf("Error: Division por cero detectada.\n");
            return;
        }
        X[i] = (B[i] - suma) / A[i][i];
    }

    // Mostrar resultados
    printf("\n** SOLUCION POR ELIMINACION GAUSSIANA **\n");
    for (int i = 0; i < n; i++) {
        printf("X[%d] = %.6lf\n", i, X[i]);
    }
}
 int paquete2(){
 int determ;
    int max_iter = 0;
    double tolerancia = 0.0;
    int edd_flag = 0;
    
    printf("\t\t\t\t\tBienvenido al programa 2.\n\tÂ°Solucion de Sistemas de Ecuaciones\n");
    printf("- Este programa leera un sistema de ecuaciones y lo resolvera usando metodos numericos\n");

    printf("\nIntroduzca el valor de la dimension de la matriz (n x n): ");
    scanf("%d", &n);

    if (n > MAX) {
        printf("Dimension demasiado grande. Maximo permitido: %d\n", MAX);
        return 1;
    }

    readdata();
    readb(); 
    printmat();


    printf("\n--- VERIFICACION DE SOLUCION (Paso 1) ---\n");
    determ = determinante(); 
    if (n == 2 && determ == 0) {
        printf("Det = 0. El sistema de ecuaciones NO tiene solucion.\n");
        return 1;
    }
    printf("Sistema con solucion (asumiendo det != 0).\n");


    printf("\n--- VERIFICACION EDD (Paso 2) ---\n");
    edd_flag = isEDD();
    
    if (edd_flag) {
        printf("\nLa matriz es EDD. Se aplica Gauss-Seidel directamente.\n");
        
        printf("\n--- Parametros para el Metodo Iterativo (Gauss-Seidel) ---\n");
        printf("Introduzca el vector inicial (X0) por componente:\n");
        for(int i = 0; i < n; i++) {
            printf("X0[%d]: ", i);
            scanf("%lf", &X[i]); 
        }
        printf("Introduzca la tolerancia (ej. 0.0001): ");
        scanf("%lf", &tolerancia);
        printf("Introduzca el maximo de iteraciones: ");
        scanf("%d", &max_iter);
        
        gaussSeidel(tolerancia, max_iter);
        
    } else {
        printf("\nADVERTENCIA: La matriz NO es EDD. La convergencia NO se garantiza.\n");
        printf("\n--- Se procede con el metodo de Eliminacion Gaussiana ---\n");
        eliminacionGaussiana();
    }

    return 0;
}
		break;
	case 3:
		printf("Factorizacion LU\n");
		break;
	case 4:
		printf("Obtencion de Eigenvalues\n");
		break;
	default:
		printf("Opcion incorrecta, vuelva a ingresar un numero por favor \n");
	break;

}
printf("Presione 's' para repetir el programa \n");
scanf("%s",&repeticion);
}while(repeticion=='s'||repeticion=='S');
return 1;
}
