
#include "Metodo.h"

int main() {
    int opcion, n;
    double Matrix[MAX][MAX], b[MAX], X[MAX], L[MAX][MAX];

    do {
        mostrarMenu();
        printf("Seleccione una opción: ");
        scanf("%d", &opcion);

        switch (opcion) {
        	
        	case 1:
		fn1(1.0, 2.0, 50, 0.0001); 
                break;

            case 2:
                fn2(1.0, 2.0, 50, 0.0001);
                break;
            case 3:
                printf("Ingrese tamaño de la matriz: ");
                scanf("%d", &n);
                leerMatriz(n, Matrix, b);
                gaussSeidel(n, Matrix, b, X, 0.0001, 100);
                imprimirVector(n, X);
                break;
            case 4:
                programa3();
                break;
            case 5: // Método de Potencia
                
                printf("Ingrese tamaño de la matriz: ");
                scanf("%d", &n);
                
                double tol;
                int max_iter;

    
                // Leer matriz A
                leerMatriz(n, Matrix, b);
                // Leer tolerancia y máximo de iteraciones
                printf("Ingrese la tolerancia: ");
                scanf("%lf", &tol);
                printf("Ingrese el máximo de iteraciones: ");
                scanf("%d", &max_iter);

                
                // Llamar a la función unificada
                calcularValoresPropios(n, Matrix, tol, max_iter);
                break;
            case 0:
                printf("Saliendo...\n");
                break;
            default:
                printf("Opción inválida.\n");
        }
    } while (opcion != 0);

    return 0;
}

