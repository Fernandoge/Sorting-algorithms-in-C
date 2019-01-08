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
#include <time.h>
#include <math.h>
// limits.h necesaria para conocer el maximo valor de int que el computador puede representar (necesario para Ordenamiento por busqueda) Ver linea 25, 95 y 237
#include <limits.h> 

#define IZQ 1
#define BAL 0
#define DER -1
#define TRUE 1
#define FALSE 0
#define VACIO -1
#define maxvalue(type) _Generic(type, int: INT_MAX)

// Estructuras utilizadas

typedef struct NodoArbol23 {
	double raiz1, raiz2;
	struct NodoArbol23 *hijo1, *hijo2, *hijo3;
} TArbol23, *Arbol23;


typedef struct NodoAVL
{
	double info;
	struct NodoAVL *izq, *der;
	int balan;
} TAVL, *AVL;

// Funciones varias
void myprintf(double x, const int is_e);
void info(char *x);
void imprimirArreglo(double *x, const unsigned int n, const int is_e);
void copiarArreglo(double *x, double *y, const unsigned int largo);
double moda(const int lon, double *v, int *frecuencia);

// MergeSort
void MergeSort(double *A, const int n);
void Merge(double *A, double *L, const int leftCount, double *R, const int rightCount);
//HeapSort
void heapify (double *Arreglo, unsigned int i, unsigned int heapsize );
void buildheap (double *Arreglo, int largo);
void Heapsort(double *Arreglo, unsigned int largo);

// AVL tree
AVL roteIzq(AVL a);
AVL roteDer(AVL a);
AVL roteDerIzq(AVL a);
AVL roteIzqDer(AVL a);
AVL balanceaDer(AVL a);
AVL balanceaIZQ(AVL a);
AVL insertar(AVL a, TAVL *p, int *masAlto);
AVL insAVL(AVL a, const double elem);
void destruirAVL(AVL a);
void mostrarAVL(AVL a, const int is_e);

// 2-3 Tree
int subirInfo1 (Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq);
int subirInfo2 (Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq);
int subirInfo3 (Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq);
int insHoja (Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq);
int insertar2 (Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq);
Arbol23 insArbol23 (Arbol23 a, double elem);
void mostrar23 (Arbol23 a,const int is_e);
void destruir23 (Arbol23 a);

// QuickSort
void swap(double* a, double* b);
void QuickSort(double a[], const int start, const int end);
double HoarePartition (double *a, const int p, const int r);

// Ordenamiento por Conteo
void counting_sort_mm(double *array, int n, double min, double max);

// MAIN
int main (int argc, char *argv[]) {  //argv[0] = ".\lab2" argv[1] = nombre archivo entrada }
	double *input, *ms, *qs, *hs, *cs, *tmp, x;
	double min, max, A, mediana;
	int cant_datos=0, cant_distintos, k; // k -> cantidad de clases
    int MAXINT = maxvalue(MAXINT);
	char check;
	int is_e = 0, mostrarTodo = 0;
	clock_t begin, end;	
	double time_spent;
	AVL avl = NULL;
	Arbol23 dostres = NULL;
  	
  	if( argc<=1 )
  		info(argv[0]);
  	else {
  		FILE *file;
	  	file = fopen(argv[1], "r");

	  	if(!file) {
	  		printf("Archivo %s no existe!\n", argv[1]);
	  		exit(0);
	  	}
	  	// chequear si la entrada es en notacion cientifica o notacion, si es la primera linea hay una 'e' entonces los es y is_check valdra 1
	  	do{
	  		check = fgetc(file);
	  		if (check == 'e'){
	  			is_e = 1;
	  			break;
	  		}

	  	} while ( check != '\n');

	  	rewind(file);
	  	input = malloc(sizeof(double));	

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

        printf("Se recomienda no mostrar los arreglos por pantalla debido al tiempo que demora en mostrarlos\n");
	  	printf("Desea mostrar los arreglos por pantalla? Si: 1 No: 0\n");
	  	do{
	  		scanf("%d", &mostrarTodo);
	  	} while (mostrarTodo!=0 && mostrarTodo!=1);
	  	
	  	if(mostrarTodo){
	  		printf("Entrada: \n");
	  		imprimirArreglo(input, cant_datos, is_e);
	  	}
	  	

	  	printf("\n////////////////////////////////////////////////////////////////////////\n");

	  	printf("QuickSort: \n");
	  	qs = malloc(cant_datos*sizeof(double));
	  	copiarArreglo(input, qs, cant_datos);
	  	begin = clock();
	  	QuickSort(qs, 0, cant_datos);
	  	end = clock();
	  	if(mostrarTodo)
	  		imprimirArreglo(qs, cant_datos, is_e);
	  	time_spent = (double)(end - begin)/CLOCKS_PER_SEC; 
		printf("\nTiempo empleado:\t%f\n", time_spent);		

	  	
	  	printf("\n////////////////////////////////////////////////////////////////////////\n");

	  	
        printf("MergeSort: \n");
	  	ms = malloc(cant_datos*sizeof(double));
	  	copiarArreglo(input, ms, cant_datos);
	  	begin = clock();
	  	MergeSort(ms, cant_datos);
	  	end = clock();
	  	if(mostrarTodo)
	  		imprimirArreglo(ms, cant_datos, is_e);
	  	time_spent = (double)(end - begin)/CLOCKS_PER_SEC; 
		printf("\nTiempo empleado:\t%f\n", time_spent);

		min = ms[0];
		max = ms[cant_datos-1];

		free(ms);

	  	printf("\n////////////////////////////////////////////////////////////////////////\n");	  	

	
	  	printf("HeapSort: \n");
	  	hs = malloc(cant_datos*sizeof(double));
	  	copiarArreglo(input, hs, cant_datos);
	  	begin = clock();
	  	Heapsort(hs, cant_datos); 
	  	end = clock();
	  	if(mostrarTodo)
	  		imprimirArreglo(hs, cant_datos, is_e);
	  	time_spent = (double)(end - begin)/CLOCKS_PER_SEC; 
		printf("\nTiempo empleado:\t%f\n", time_spent);	

		free(hs);
		rewind(file);
	  	printf("\n////////////////////////////////////////////////////////////////////////\n");	
	  	
	  	
	  	printf("AVL Tree: \n");	  	
	  	begin = clock();
	  	while (fscanf(file, "%lf", &x) != EOF)	  		
	  		avl = insAVL(avl, x);
	  	end = clock();
	  	if(mostrarTodo)
	  		mostrarAVL(avl, is_e);
	  	time_spent = (double)(end - begin)/CLOCKS_PER_SEC; 			
		printf("\nTiempo empleado:\t%f\n", time_spent);	
		
		destruirAVL(avl);
		rewind(file);
		printf("\n////////////////////////////////////////////////////////////////////////\n");
		

		printf("2-3 Tree: \n");
	  	begin = clock();
	  	while (fscanf(file, "%lf", &x) != EOF)	  		
	  		dostres = insArbol23(dostres, x);
	  	end = clock();
	  	if(mostrarTodo)	
	  		mostrar23(dostres, is_e);
	  	time_spent = (double)(end - begin)/CLOCKS_PER_SEC; 
	  	printf("\nTiempo empleado:\t%f\n", time_spent);

		destruir23(dostres);		
		printf("\n////////////////////////////////////////////////////////////////////////\n");

        if( max > MAXINT) //Verifica el maximo valor de int para no salir del rango
		{
			printf("No es posible realizar Ordenamiento por Conteo dado que el elemento maximo ");
			myprintf(max, is_e);
			printf("\nes mayor al maximo entero representable (%d) en este computador\n", MAXINT);
		}

		else {
			printf("Counting Sort: \n");
			cs = malloc(cant_datos*sizeof(double));
		  	copiarArreglo(input, cs, cant_datos);
		  	begin = clock();
		  	counting_sort_mm(cs, cant_datos, min, max);
		  	end = clock();
		  	if(mostrarTodo)
		  		imprimirArreglo(cs, cant_datos, is_e);
		  	time_spent = (double)(end - begin)/CLOCKS_PER_SEC; 
			printf("\nTiempo empleado:\t%f\n", time_spent);
			free(cs);
		}

		free(input);
		printf("\n////////////////////////////////////////////////////////////////////////\n");

        printf("\n");
		printf("Cantidad de datos: %d\n", cant_datos);

		printf("Minimo: ");
		myprintf(min, is_e);
		printf("\n");

		printf("Maximo: ");
		myprintf(max, is_e);
		printf("\n");

		//Elementos distintos
		cant_distintos = 0;
		for(int i=0; i<cant_datos; i++) {
			if(qs[i]==qs[i+1])
				continue;
			else
				cant_distintos++;
		}
		printf("Cantidad de elementos distintos: %d\n", cant_distintos);

		// Calculo de la moda
		int frec = 1; //para obtener la moda
		printf("Moda: %lf ", moda(cant_datos, qs, &frec));
		printf("Frecuencia: %d\n", frec);

		// Calculo de la mediana, se usa el arreglo ya ordenado qs[]
		printf("Mediana: ");
		if(!cant_datos%2) // cant_datos es par
			mediana = (qs[(cant_datos/2)-1]+qs[cant_datos/2])/2;
		else // cant_datos es impar
			mediana = qs[cant_datos/2];
		myprintf(mediana, is_e);
		printf("\n");


		// Calculo de cuartiles
		printf("Cuartiles: \n");	
		printf("Q1: ");
		myprintf(qs[(cant_datos/4)-1], is_e);
		printf("\n");
		printf("Q2: ");
		myprintf(mediana, is_e);
		printf("\n");
		printf("Q3: ");
		myprintf(qs[((3*cant_datos)/4)-1], is_e);
		printf("\n");


        printf("\n\n");

	  	free(qs);	  	
	  	fclose(file);
	  	return EXIT_SUCCESS;
	

  	} 
}

// END MAIN

////////////////////////////////////////////

double moda(const int lon, double *v, int *frecuencia)
{
	int i=0, cont=1;
	float moda;

	moda=v[i];

	while(i<lon) {		
		if(v[i]==v[i+1]) {
			cont++;
			i++;
		}
		else {
			if(cont>(*frecuencia)) {
				moda=v[i];
				*frecuencia=cont;
				cont=1;
				i++;
			} else {
				cont=1;
				i++;
			}
		}
	}
	return(moda);
}

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

/////	QuickSort	////////////////////////

void swap(double* a, double* b) {
    double t = *a;

    *a = *b;
    *b = t;
}

double HoarePartition (double *a, const int p, const int r) {
    double x=a[p];
    int i = p-1, j=r;

    while (1) {
        do
        	j--;
        while(a[j] > x);

        do
        	i++;
        while(a[i] < x);

        if(i < j)
        	swap(&a[i],&a[j]);
        else 
            return j+1;
    }
}

void QuickSort(double a[], const int start, const int end) {
    double q;

    if (end-start<2)
    	return;
    q = HoarePartition(a, start, end);
    QuickSort(a, start, q);
    QuickSort(a, q, end);
}

////////////////////////////////////////////

/////	HeapSort	////////////////////////

void heapify (double *Arreglo, const unsigned int i, const unsigned int heapsize ) {
	unsigned int derecha, izquierda, mayor;
	double temp;
   
	izquierda = 2*i+1;
	derecha = 2*i+2;

	if (izquierda<=heapsize && *(Arreglo+izquierda)>*(Arreglo+i)) // *(Arreglo+x) = Arreglo[x]
		mayor=izquierda;
	else
		mayor=i;
	if (derecha<=heapsize && *(Arreglo+derecha)>*(Arreglo+mayor))
	 	mayor=derecha;
	if (mayor != i) {
		temp = *(Arreglo+i);
		*(Arreglo+i) = *(Arreglo+mayor);
		*(Arreglo+mayor) = temp;
		heapify (Arreglo, mayor, heapsize);
	}
	return;
}
 
void buildheap (double *Arreglo, const int largo) {
	int i;

	for (i=(int) ((largo/2)-1); i>=0; i--)
		heapify(Arreglo, i, largo);
   	return;
}
 
void Heapsort(double *Arreglo, unsigned int largo) {
	unsigned int i, heapsize;              
	double temp;
         
	largo=largo-1;
	heapsize=largo;
	     
	buildheap(Arreglo, largo);
	for (i = largo;i>=1;i--) {
		temp=*(Arreglo+0);
		*(Arreglo+0)=*(Arreglo+i);
		*(Arreglo+i)=temp;
		heapsize=heapsize-1;
		heapify(Arreglo,0,heapsize);
	}
	return;
}

////////////////////////////////////////////

/////	 AVL	////////////////////////////

AVL roteIzq(AVL a) {
	AVL temp = a->der;
	a->der = temp->izq;
	temp->izq = a;
	return temp;

}

AVL roteDer(AVL a) {
	AVL temp = a->izq;
	a->izq = temp->der;
	temp->der = a;
	return temp;

}

AVL roteDerIzq(AVL a) {
	a->der = roteDer(a->der);
	return roteIzq(a);
}

AVL roteIzqDer(AVL a) {
	a->izq = roteIzq(a->izq);
	return roteDer(a);
}

AVL balanceaDer(AVL a) {
	if( a->der->balan == DER){
		a->balan = a->der->balan = BAL;
		a = roteIzq(a);
	}
	else{
		switch(a->der->izq->balan){
			case IZQ:
				a->balan = BAL;
				a->der->balan = DER;
				break;
			case BAL:
				a->balan = a->der->balan = BAL;
				break;
			case DER:
				a->balan = IZQ;
				a->der->balan = BAL;
				break;
		}
		a->der->izq->balan = BAL;
		a = roteDerIzq(a);
	}
	return a;
}

AVL balanceaIZQ(AVL a) {
	if( a->izq->balan == IZQ){
		a->balan = a->izq->balan = BAL;
		a = roteDer(a);
	}
	else{
		switch(a->izq->der->balan){
			case DER:
				a->balan = BAL;
				a->izq->balan = IZQ;
				break;
			case BAL:
				a->balan = a->izq->balan = BAL;
				break;
			case IZQ:
				a->balan = DER;
				a->izq->balan = BAL;
				break;
		}
		a->izq->der->balan = BAL;
		a = roteIzqDer(a);
	}
	return a;
}

AVL insertar(AVL a, TAVL *p, int *masAlto) {
	if(a == NULL){
		*masAlto = TRUE;
		a = p;
	}
	else{
		if (a->info > p->info){
			a->izq = insertar(a->izq, p, masAlto);
			if(*masAlto)
				switch(a->balan){
					case IZQ:
						*masAlto = FALSE;
						a = balanceaIZQ(a);
						break;
					case BAL:
						a->balan = IZQ;
						break;
					case DER:
						*masAlto = FALSE;
						a->balan = BAL;
						break;
				}
			/* code */
		} else{
			a->der = insertar(a->der, p, masAlto);
			if(*masAlto)
				switch(a->balan){
					case IZQ:
						*masAlto = FALSE;
						a->balan = BAL;
						break;
					case BAL:
						a->balan = DER;
						break;
					case DER:
						*masAlto = FALSE;
						a = balanceaDer(a);
				}
		}
	}//else principal
		
	return a;
}

AVL insAVL(AVL a, const double elem) {
	AVL p = (AVL)malloc(sizeof(TAVL));
	int masAlto;
	p->izq = NULL;
	p->der = NULL;
	p->info = elem;
	p->balan = BAL;

	return insertar(a, p, &masAlto);
}

void destruirAVL(AVL a) {
	if( a != NULL){
		destruirAVL(a->izq);
		destruirAVL(a->der);
		free(a);
	}
}

void mostrarAVL(AVL a, const int is_e) {
	/** recorrer inorden: izq - raiz - der **/
	if( a != NULL){
		mostrarAVL(a->izq, is_e);
		myprintf(a->info, is_e);
		mostrarAVL(a->der, is_e);
	}

}

////////////////////////////////////////////

/////	2-3						////////////
//static functions are functions that are only visible to other functions in the same file (more precisely the same translation unit).

int subirInfo1(Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq) {
  double temp;
  Arbol23 nuevo;
  if (a->raiz2 == VACIO)
/* hay campo en ese nodo: reorganizar */
    {
      a->raiz2 = a->raiz1;
      a->hijo3 = a->hijo2;
      a->raiz1 = *elem;
      a->hijo1 = *arbolIzq;
      a->hijo2 = *arbolDer;
      return FALSE;
    }
  else
/* no hay campo en el nodo: partir y subir */
    {
      nuevo = (Arbol23) malloc (sizeof (TArbol23));
      nuevo->hijo3 = NULL;
      nuevo->raiz2 = VACIO;
      nuevo->raiz1 = a->raiz2;
      nuevo->hijo1 = a->hijo2;
      nuevo->hijo2 = a->hijo3;
      temp = *elem;
      *elem = a->raiz1;
      a->raiz1 = temp;
      a->raiz2 = VACIO;
      a->hijo1 = *arbolIzq;
      a->hijo2 = *arbolDer;
      a->hijo3 = NULL;
      *arbolIzq = a;
      *arbolDer = nuevo;
      return TRUE;
    }
}

int subirInfo2(Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq) {
  Arbol23 nuevo;
  if (a->raiz2 == VACIO)
/* hay campo en ese nodo: reorganizar */
    {
      a->raiz2 = *elem;
      a->hijo2 = *arbolIzq;
      a->hijo3 = *arbolDer;
      return FALSE;
    }
  else
/* no hay campo en el nodo: partir y subir */
    {
      nuevo = (Arbol23) malloc (sizeof (TArbol23));
      nuevo->hijo3 = NULL;
      nuevo->raiz2 = VACIO;
      nuevo->raiz1 = a->raiz2;
      nuevo->hijo1 = *arbolDer;
      nuevo->hijo2 = a->hijo3;
      a->hijo2 = *arbolIzq;
      a->hijo3 = NULL;
      a->raiz2 = VACIO;
      *arbolIzq = a;
      *arbolDer = nuevo;
      return TRUE;
    }
}

int subirInfo3(Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq) {
  Arbol23 nuevo;
  nuevo = (Arbol23) malloc (sizeof (TArbol23));
  nuevo->hijo3 = NULL;
  nuevo->raiz2 = VACIO;
/* No hay campo en el nodo: tiene que seguir subiendo */
  nuevo->raiz1 = *elem;
  nuevo->hijo1 = *arbolIzq;
  nuevo->hijo2 = *arbolDer;
  *elem = a->raiz2;
  a->raiz2 = VACIO;
  a->hijo3 = NULL;
  *arbolIzq = a;
  *arbolDer = nuevo;
  return TRUE;
}

int insHoja(Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq) {
  Arbol23 nuevo;
  if (a->raiz2 == VACIO)
/* caso 1: hay espacio en el nodo */
    {
      if (*elem < a->raiz1)
  {
    a->raiz2 = a->raiz1;
    a->raiz1 = *elem;
  }
      else
  a->raiz2 = *elem;
      return FALSE;
    }
  else
/*caso 2: no hay espacio en el nodo y se debe partir */
    {
      nuevo = (Arbol23) malloc (sizeof (TArbol23));
      nuevo->hijo1 = nuevo->hijo2 = nuevo->hijo3 = NULL;
      nuevo->raiz2 = VACIO;
      if (*elem < a->raiz1)
/* sube la raC-z 1 */
  {
    nuevo->raiz1 = *elem;
    *elem = a->raiz1;
    a->raiz1 = a->raiz2;
    a->raiz2 = VACIO;
  }
      else if (*elem < a->raiz2)
/* sube el elemento nuevo */
  {
    nuevo->raiz1 = a->raiz1;
    a->raiz1 = a->raiz2;
    a->raiz2 = VACIO;
  }
      else
/* sube la raC-z 2 */
  {
    nuevo->raiz1 = a->raiz1;
    a->raiz1 = *elem;
    *elem = a->raiz2;
    a->raiz2 = VACIO;
  }
      *arbolIzq = nuevo;
      *arbolDer = a;
      return TRUE;
    }
}

int insertar2(Arbol23 a, double *elem, Arbol23 *arbolDer, Arbol23 *arbolIzq) {
  if (a->hijo1 == NULL)
    return insHoja (a, elem, arbolDer, arbolIzq);
  else if (*elem < a->raiz1)
    return insertar2 (a->hijo1, elem, arbolDer, arbolIzq) ? subirInfo1 (a, elem, arbolDer, arbolIzq) : FALSE;
  else if (a->raiz2 == VACIO || *elem < a->raiz2)
    return insertar2 (a->hijo2, elem, arbolDer, arbolIzq) ? subirInfo2 (a, elem, arbolDer, arbolIzq) : FALSE;
  else
    return insertar2 (a->hijo3, elem, arbolDer, arbolIzq) ? subirInfo3 (a, elem, arbolDer, arbolIzq) : FALSE;
}


Arbol23 insArbol23(Arbol23 a, double elem) {
  Arbol23 derecho = NULL, izquierdo = NULL, nuevo;
  if (a == NULL)
    {
      a = (Arbol23) malloc (sizeof (TArbol23));
      a->raiz1 = elem;
      a->raiz2 = VACIO;
      a->hijo1 = a->hijo2 = a->hijo3 = NULL;
      return a;
    }
  else if (insertar2 (a, &elem, &derecho, &izquierdo))
    {
      nuevo = (Arbol23) malloc (sizeof (TArbol23));
      nuevo->raiz1 = elem;
      nuevo->raiz2 = VACIO;
      nuevo->hijo1 = izquierdo;
      nuevo->hijo2 = derecho;
      nuevo->hijo3 = NULL;
      return nuevo;
    }
  else
    return a;
}

void mostrar23(Arbol23 a, const int is_e) {

  if ((a->hijo1 == NULL) && (a->hijo2 == NULL) && (a->hijo3 == NULL)) {

      /* Here print node's value/values */
      if (a->raiz1 != -1)
        myprintf(a->raiz1, is_e);
      if (a->raiz2 != -1)
      	myprintf(a->raiz2, is_e);

  } 
  else {
      if ((a->hijo1 != NULL) && (a->hijo2 != NULL) && (a->hijo3 != NULL)) {
        mostrar23 (a->hijo1, is_e);
        if (a->raiz1 != -1)
          myprintf(a->raiz1, is_e);
        mostrar23 (a->hijo2, is_e);
        if (a->raiz2 != -1)
          myprintf(a->raiz2, is_e);
        mostrar23 (a->hijo3, is_e);
      }
      else {     // 2 hijos
        mostrar23 (a->hijo1, is_e);
        if (a->raiz1 != -1)
          myprintf(a->raiz1, is_e);
        if (a->raiz2 != -1)
          myprintf(a->raiz2, is_e);
        mostrar23 (a->hijo2, is_e);
      }
  }
}

void destruir23(Arbol23 a) { 

  if ((a->hijo1 == NULL) && (a->hijo2 == NULL) && (a->hijo3 == NULL)) {
    free(a);
  } 
  else {
    if ((a->hijo1 != NULL) && (a->hijo2 != NULL) && (a->hijo3 != NULL)) {
      destruir23 (a->hijo1);
      destruir23 (a->hijo2);
      destruir23 (a->hijo3);
      free(a);
    }
    else {     // 2 hijos
      destruir23 (a->hijo1);
      destruir23 (a->hijo2);
      free(a);
    }
  }
}

////////////////////////////////////////////

/////	Counting Sort	////////////////////

void counting_sort_mm(double *array, int n, double min, double max)
{
  int i, j, z;
 
  double range = max - min + 1;
  double *count = malloc(range * sizeof(*array));
 
  for(i = 0; i < range; i++) count[i] = 0;
  for(i = 0; i < n; i++) count[ (int)array[i] - (int)min ]++; //se hace el casteo a int (trunca los decimales)
 
  for(i = min, z = 0; i <= max; i++) {
    for(j = 0; j < count[i - (int)min]; j++) {
      array[z++] = i;
    }
  } 
 
  free(count);
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
