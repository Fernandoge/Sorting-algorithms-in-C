/*
*Fernando García Polgatti
*Felipe Parra
*
* Entorno de prueba:
* Ubuntu 18.04
* Compilador: GCC 7.2.1
* C11 (C standard revision)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> 



// Funciones varias
void myprintf(double x, const int is_e);
void info(char *x);
void imprimirArreglo(double *x, const unsigned int n, const int is_e);
void copiarArreglo(double *x, double *y, const unsigned int largo);

// MergeSort
void MergeSort(double *A, const int n);
void Merge(double *A, double *L, const int leftCount, double *R, const int rightCount);


// MAIN
int main (int argc, char *argv[]) {  //argv[0] = ".\lab2" argv[1] = nombre archivo entrada }
	double *input, *ms, *qs, *hs, *cs, *tmp, x;
	double min, max, A, mediana;
	int cant_datos=0, cant_distintos, k; // k -> cantidad de clases
	char check;
	int is_e = 0, full_out = 0;
  	
  	if( argc<=1 )
  		info(argv[0]);
  	else {
  		FILE *file;
	  	file = fopen(argv[1], "r");

	  	if(!file) {
	  		printf("Archivo %s no existe!\n", argv[1]);
	  		exit(0);
	  	}
	  	do{
	  		check = fgetc(file);
	  		if (check == 'e'){
	  			is_e = 1;
	  			break;
	  		}

	  	} while ( check != '\n');

	  	rewind(file);
	  	input = malloc(sizeof(double));	

	  	//Lee los datos y los guarda en el arreglo input[]
	  	printf("Leyendo: %s\n", argv[1]);

	  	while (fscanf(file, "%lf", &x) != EOF) {	  		
	  		input[cant_datos] = x;	
	  		cant_datos++;
	  		tmp = realloc(input, (cant_datos+1)*sizeof(double));
	  		if (!tmp)
			{
			  printf("No se pudo asignar memoria\n");
			  exit(0);
			}
			else
			{
				input = tmp;
			}
	  	}

	  	printf("\n---------------------------------------------------------------------------------------\n");
  		printf("Entrada: \n\n");
  		imprimirArreglo(input, cant_datos, is_e);
	  	

	  	printf("\n---------------------------------------------------------------------------------------\n");

	  	

	  	printf("MergeSort: \n\n");
	  	ms = malloc(cant_datos*sizeof(double));
	  	copiarArreglo(input, ms, cant_datos);
	  	MergeSort(ms, cant_datos);
  		imprimirArreglo(ms, cant_datos, is_e);

		min = ms[0];
		max = ms[cant_datos-1];

		free(ms);
	  	printf("\n---------------------------------------------------------------------------------------\n");

	  	
	  }
}

// END MAIN

////////////////////////////////////////////

/////	MergeSort	////////////////////////

void Merge(double *A, double *L, const int leftCount, double *R, const int rightCount) {
	// i - indice subarreglo izquierdo (L)
	// j - indice subarreglo derecho(R)
	// k - indice subarreglo merge (A)

	int i = 0, j = 0, k =0;

	while(i<leftCount && j< rightCount) {
		if(L[i]  < R[j])
			A[k++] = L[i++];
		else
			A[k++] = R[j++];
	}
	while(i < leftCount)
		A[k++] = L[i++];
	while(j < rightCount)
		A[k++] = R[j++];
}

void MergeSort(double *A, const int n) {
	double *L, *R;
	int i, mid;

	if(n < 2)
		return; 

	mid = n/2;  // punto medio del arreglo
	
	L = (double*)malloc(mid*sizeof(double)); // crea subarreglo izquierdo
	R = (double*)malloc((n - mid)*sizeof(double)); // crea subarreglo derecho
	
	for(i = 0;i<mid;i++) 
		L[i] = A[i]; 
	for(i = mid;i<n;i++)
		R[i-mid] = A[i]; 

	MergeSort(L,mid);  // ordena lado izquierdo
	MergeSort(R,n-mid);  // ordena lado derecho
	Merge(A, L, mid, R, n-mid);  // realiza el merge
    free(L);
    free(R);
}
////////////////////////////////////////////

/////	Funciones Auxiliares	////////////

void myprintf(double x, const int is_e) {
	if(is_e)
			printf("%e ", x);
		else
			printf("%.2lf ", x);
}

void info(char *x) {
	printf("Archivo de entrada no especificado, ejecutar de la siguiente manera (sin comillas):\n");
	printf("%s \"nombre_archivo\"\n", x);
	printf("Ejemplo:\n");
	printf("%s datos.tex\n", x);
}

void imprimirArreglo(double *x, const unsigned int n, const int is_e) {
	unsigned int i=0;
	
	for(i=0; i<n; i++)
		myprintf(x[i], is_e);
	printf("\n");
}

void copiarArreglo(double *x, double *y, const unsigned int largo) {
	unsigned int i=0;

	for(i=0; i<largo; i++)
		y[i]=x[i];
}


///////////////////////////////////////////////////
