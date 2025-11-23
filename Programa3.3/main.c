
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
                printf("Ingrese tamaño de la matriz: ");
                scanf("%d", &n);
                leerMatriz(n, Matrix, b);
                if (cholesky(n, Matrix, L)) {
                    printf("Factorización Cholesky exitosa.\n");
                }
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

