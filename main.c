//Chivu Constantin Razvan , 332 CA
//Tema 2 ASC

#include <stdlib.h>
#include <cblas.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define TRI 1
#define NTRI -1
#define MAX 50

int n, m;


/* type = TRI pentru a genera A , type = NTRI pentru a genera B */
void genMat(int m, int n, double* mat, int type){

    for (int i = 0; i < m; ++i){

        for (int j = 0; j < n; ++j){

            if (type == TRI){
                if (j >= i){
                    mat[i*n + j] = 1;
                }
            }
            else{
                mat[i*n + j] = rand() % MAX;
            }

        }

    }

}


void print(int m, int n, double*mat){

    for (int i = 0; i < m; ++i){

        for (int j = 0; j < n; ++j){
            printf("%.2lf ", mat[i*n + j]);
        }
        printf("\n");
    }
    printf("\n");
}


void sst(int m, int n, double*a, double*b, double*x){


    double sum = 0;

    for (int j = 0; j < n; ++j){
        x[(m - 1)*n + j] = b[(m - 1)*n + j] / a[(m - 1) * m + (m - 1)];
    }


    for (int k = 0; k < n; ++k){

        /* for each column of Xs */
        /* apply single-col sst */

        for (int i = m - 2; i >= 0; --i){
            sum = 0;

            for (int j = i + 1; j < m; ++j){

                sum += a[i*m + j] * x[j*n + k];

            }

            x[i*n + k] = (b[i*n + k] - sum) / a[i*m + i];
        }

    }


}


void better_sst(int m, int n, double*a, double*b){


    double sum = 0;

    double *pa, *pb, *pk, *pa2;

    int cnt = m - 2;

    double* lim;

    for (pk = b; pk < b + n; pk++){

        *(pk + (m - 1)*n) = *(pk + (m - 1)*n) / *(a + (m - 1)*(m - 1));

        cnt = m - 2;

        for (pa = a + (m - 2)*m; pa >= a; pa -= m, cnt--){

            if (pa < a) break;

            sum = 0;

            lim = pk + m*n;

            for (pb = pk + (cnt + 1)*n; pa2 = pa + cnt + 1; pb += n, pa2++){

                if (pb == lim) break;

                sum += *pa2 * *pb;

            }

            *(pk + cnt*n) = (*pk + cnt*n - sum) / *(pa + cnt);

        }

    }

}


/* pentru a verifica corectitudinea */
void ver(int m, int n, double*l, double*r){

    for (int i = 0; i < m; ++i){

        for (int j = 0; j < n; ++j){

            if (abs(l[i*n + j] - r[i*n + j]) > 0.00001){
                printf("DIFF %.10lf \n", l[i*n + j] - r[i*n + j]);
                return;
            }

        }

    }

    printf("OK\n");

}


void cpy(int m, int n, double* from, double* to){

    for (int i = 0; i < m; ++i){

        for (int j = 0; j < n; ++j){

            to[i*n + j] = from[i*n + j];

        }

    }

}


void readSize(char*filen){

    scanf("%d %d", &m, &n);

}


void my_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N,
    const double alpha, const double  *A, const int lda, double  *B, const int ldb){

    double* x = (double*)calloc(M*N, sizeof(double));

    sst(M, N, A, B, x);

    for (int i = 0; i < M; ++i){

        for (int j = 0; j < N; ++j){
            B[i*N + j] = x[i*N + j];
        }

    }


}


void my_better_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N,
    const double alpha, const double  *A, const int lda, double  *B, const int ldb){


    better_sst(M, N, A, B);


}



int main(){

    srand(time(NULL));
    clock_t start, end;
    float seconds;
    double alf = 1.0;

    FILE *t1, *t2, *t3;

#ifdef OPTERON
    printf("OPTERON\n");
    t1 = fopen("time_normal_o.txt", "aw");
    t2 = fopen("time_opt_o.txt", "aw");
    t3 = fopen("time_blas_o.txt", "aw");
#endif
#ifdef NEHALEM
    printf("NEHALEM\n");
    t1 = fopen("time_normal.txt", "aw");
    t2 = fopen("time_opt.txt", "aw");
    t3 = fopen("time_blas.txt", "aw");
#endif

    printf("starting reading..\n");
    readSize("size.txt");
    printf("read m = %d , n = %d \n", m, n);

    printf("alloc..\n");


    /* alocari */
    double* mat = calloc(m*m, sizeof(double)); /* A matrix */
    double* mat2 = calloc(m*n, sizeof(double)); /* B matrix */
    double* res = calloc(m*n, sizeof(double)); /* result */
    double* res_cbls = calloc(m*n, sizeof(double));

    printf("ready..\n");


    /* genereaza A */
    genMat(m, m, mat, TRI);


    /* genereaza  B */
    genMat(m, n, mat2, NTRI);

    cpy(m, n, mat2, res);

    /* rezolva sistemul "de mana " si afiseaza rezultatul */
    start = clock();
    my_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, alf, mat, m, mat2, n);
    end = clock();
    seconds = (float)(end - start) / CLOCKS_PER_SEC;
    fprintf(t1, "%f %d normal\n", seconds, m*n);



    /* rezolva sistemul 'de mana imbunatatit' si afiseaza rezultatul */
    start = clock();
    my_better_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, alf, mat, m, res, n);
    end = clock();
    seconds = (float)(end - start) / CLOCKS_PER_SEC;
    fprintf(t2, "%f %d optimizat\n", seconds, m*n);


    /* rezolva sistemul cu BLAS si afiseaza rezultatul */
    start = clock();
    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, alf, mat, m, res_cbls, n);
    end = clock();
    seconds = (float)(end - start) / CLOCKS_PER_SEC;
    fprintf(t3, "%f %d blas\n", seconds, m*n);


    /* verifica sa fie rezultatele identice */
    ver(m, n, mat2, res);
    ver(m, n, res, res_cbls);


    free(mat);
    free(mat2);
    free(res);
    free(res_cbls);

    fclose(t1);
    fclose(t2);
    fclose(t3);

    return 0;
}
