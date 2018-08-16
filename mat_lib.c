#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>






typedef struct double_mat{
	double **vals;
	unsigned int axis0;
	unsigned int axis1;
	unsigned int length;
	double sum;
	double avg;
}DoubleMat;



void free_nm(double **vals, unsigned int n){
	unsigned int i;
	for (i = 0; i < n; i++)
		free(vals[i]);
	free(vals);
}

void freeMat(DoubleMat *mat){
	free_nm(mat->vals, mat->axis0);
	free(mat);
}



DoubleMat* createMat(unsigned int n, unsigned int m, double **values){
	DoubleMat *new_mat = (DoubleMat*)malloc(sizeof(DoubleMat));

	new_mat->vals = values;
	new_mat->axis0 = n ;
	new_mat->axis1 = m;
	return new_mat;

}

double max(DoubleMat *mat){
	unsigned int i, j, max_i = 0, max_j = 0;
	double curr, max = mat->vals[0][0];
	for (i = 0; i < mat->axis0; i++)
		for(j = 0; j < mat->axis1; j++){
			curr = mat->vals[i][j];
			if (curr > max){
				max = curr;
				max_i = i;
				max_j = j;
			}
		}
	return max;

}

void updateMat(DoubleMat *mat){
	unsigned int i, j;
	int sum = 0;

}

void printMat(DoubleMat *mat){
	unsigned int i, j;
	printf("\n");
	for(i=0; i<mat->axis0; i++){
		printf("\t");
		for (j = 0; j < mat->axis1; j++)
			printf("%.2lf ", mat->vals[i][j]);

		printf("\n");
	}
}



double dotProd(DoubleMat *arr1, DoubleMat *arr2){
	unsigned int i, j;
	double sum = 0;

	if (arr1->axis1 != arr2->axis1){
		printf("\ndimension mismatch");
		exit(1);
	}
	for(i=0; i<arr1->axis1; i++){
		sum += arr1->vals[0][i] * arr2->vals[0][i];
	}

	return sum;
}


DoubleMat *transpose(DoubleMat *mat){
	unsigned int i, j;
	DoubleMat *res = (DoubleMat *) malloc(sizeof(DoubleMat));
	res->axis0 = mat->axis1;
	res->axis1 = mat->axis0;
	double **res_vals;
	res_vals = (double **) malloc(sizeof(double) * mat->axis1);
	for (i = 0; i < res->axis0; i++){
		res_vals[i] = (double *) malloc(sizeof(double) * res->axis1);
		for (j = 0; j < res->axis1; j++){
			res_vals[i][j] = mat->vals[j][i];
		}
	}
	res->vals=res_vals;
	return res;

}


// NOT DONE
DoubleMat *slice(DoubleMat *mat, unsigned int s0, unsigned int e0, unsigned int s1, unsigned int e1){
	unsigned int i,j;
	double **res_vals;
	unsigned int new_size0 = e0-s0;
	unsigned int new_size1 = e1-s1;
	DoubleMat *res;
	res_vals = (double **) malloc(sizeof(double *) * new_size0);
	for (i = 0; i < new_size0; i++){
		res_vals[i] = (double *) malloc(sizeof(double) * new_size1);
		for (j = 0; j < new_size1; j++){
			res_vals[i][j]=mat->vals[i+s0][j+s1];
		}
	}
	res = createMat(new_size0, new_size1, res_vals);
	return res;

}


double **alloc_nm(unsigned int n, unsigned int m){
	unsigned int i;
	double **mat = (double**)malloc(sizeof(double*)*n);
	for (i = 0; i < n; i++)
		mat[i] = (double *) malloc(sizeof(double) * m);
	return mat;
}


DoubleMat *matrixMult(DoubleMat *mat1, DoubleMat *mat2){

	unsigned int i, j, k;
	DoubleMat *res;
	double **res_vals, **v1, **v2;
	v1 = mat1->vals;
	v2 = mat2->vals;
	if (mat1->axis1 != mat2->axis0){
		printf("dimension mismatch");
		return NULL;
	}
	
	res_vals = alloc_nm(mat1->axis0, mat2->axis1);
	for (i=0; i<mat1->axis0; i++){
		for (j = 0; j<mat2->axis1; j++){
			res_vals[i][j] = 0;
			for (k = 0; k<mat1->axis1; k++){
				res_vals[i][j]+=v1[i][k]*v2[k][j];
			}


			
		}
	}
	res = createMat(mat1->axis0, mat2->axis1, res_vals);
	return res;
}

DoubleMat *concatenate(DoubleMat *mat1, DoubleMat *mat2, unsigned int axis){
	unsigned int i, j;
	unsigned int new_size0, new_size1;
	DoubleMat *res = (DoubleMat *) malloc(sizeof(DoubleMat));
	double **res_vals;

	if (axis == 0){
		if (mat1->axis1 != mat2->axis1)
			printf("dimension mismatch\n");
		else{
			new_size0 = mat1->axis0 + mat2->axis0;
			new_size1 = mat1->axis1; 
		}
	}

	else{
		if (mat1->axis0 != mat2->axis0)
			printf("dimension mismatch\n");
		else{
			new_size0 = mat1->axis0;
			new_size1 = mat1->axis1 + mat2->axis1; 
		}
	}
	res_vals = alloc_nm(new_size0, new_size1);

	for (i = 0; i < new_size0; i++){
		for (j = 0; j < new_size1; j++){
			if (i >= mat1->axis0)
				res_vals[i][j] = mat2->vals[i - mat1->axis0][j];
			else
				res_vals[i][j] = mat1->vals[i][j];
		}
	}
	res->axis0 = new_size0;
	res->axis1 = new_size1;
	res->vals = res_vals;

	return res;	
}


DoubleMat *ID_n(unsigned int n){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(n, n);
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (j == i)
			{
				res_vals[i][j] = 1;
			}
			else
				res_vals[i][j] = 0;
		}
	}
	res = createMat(n, n, res_vals);
	return res;
}

DoubleMat *ones(unsigned int n, unsigned int m){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(n, m);
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++)
			res_vals[i][j] = 1;

	}
	res = createMat(n, m, res_vals);
	return res;

}

DoubleMat *zeros(unsigned int n, unsigned int m){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(n, m);
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++)
			res_vals[i][j] = 0;

	}
	res = createMat(n, m, res_vals);
	return res;

}

/*
	mat1 - mat2
*/
DoubleMat *matSub(DoubleMat *mat1, DoubleMat *mat2){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(mat1->axis0, mat1->axis1);
	res = createMat(mat1->axis0, mat1->axis1, res_vals);
	for (i = 0; i < res->axis0; i++){
		for (j = 0; j < res->axis1; j++)
			res_vals[i][j] = mat1->vals[i][j] - mat2->vals[i][j];

	}
	
	return res;

}


DoubleMat *matAdd(DoubleMat *mat1, DoubleMat *mat2){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(mat1->axis0, mat1->axis1);
	res = createMat(mat1->axis0, mat1->axis1, res_vals);
	for (i = 0; i < res->axis0; i++){
		for (j = 0; j < res->axis1; j++)
			res_vals[i][j] = mat1->vals[i][j] + mat2->vals[i][j];

	}
	return res;

}


DoubleMat *inverseMat(DoubleMat *mat){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(mat->axis0, mat->axis1);
	res = createMat(mat->axis0, mat->axis1, res_vals);
	for (i = 0; i < res->axis0; i++){
		for (j = 0; j < res->axis1; j++)
			res_vals[i][j] = 1 / mat->vals[i][j];

	}

	return res;
}

DoubleMat *scalarMult(DoubleMat *mat, double c){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(mat->axis0, mat->axis1);
	res = createMat(mat->axis0, mat->axis1, res_vals);
	for (i = 0; i < res->axis0; i++){
		for (j = 0; j < res->axis1; j++)
			res_vals[i][j] = mat->vals[i][j] * c;

	}
	return res;
}

DoubleMat *matExp(DoubleMat *mat){
	unsigned int i, j;
	DoubleMat *res;
	double **res_vals = alloc_nm(mat->axis0, mat->axis1);

	for (i = 0; i < mat->axis0; i++){
		for (j = 0; j < mat->axis1; j++)
			res_vals[i][j] = exp(mat->vals[i][j]);
	}

	res = createMat(mat->axis0, mat->axis1, res_vals);
	return res;
}

DoubleMat *hypLinear(DoubleMat *X1, DoubleMat *theta){
	DoubleMat *res = matrixMult(theta, X1);
	return res;
}

DoubleMat *hypLog(DoubleMat *X1, DoubleMat *theta){
	DoubleMat *step1 = matrixMult(theta, X1);
	DoubleMat *step2 = scalarMult(step1, -1);
	freeMat(step1);
	DoubleMat *step3 = matExp(step2);
	freeMat(step2);
	DoubleMat *step4 = matAdd(ones(1, step3->axis1), step3);
	freeMat(step3);
	DoubleMat *res = inverseMat(step4);
	freeMat(step4);
	return res;
}





void doIteration(DoubleMat *X1, DoubleMat *y, DoubleMat *theta, double alpha){
	// printMat(X1);
	// printMat(y);
	// printMat(theta);
	// printf("\n%lf\n", alpha);
	unsigned int i, j, p;
	double grad;
	
	for (j = 0; j < X1->axis1; j++){
		DoubleMat *h_th_x = hypLinear(X1, theta);
		// printf("h_th_x is: \n");
		// printMat(h_th_x);
		for(i = 0; i < theta->axis1; i++){
			grad = (y->vals[0][j] - h_th_x->vals[0][j]) * X1->vals[i][j];
			// printf("grad is: %lf\n", grad);
			theta->vals[0][i] = theta->vals[0][i] + alpha * grad;
			// printMat(X1);
			// printMat(y);
			// printMat(theta);
		}
		freeMat(h_th_x);
	}
	
	
}



void populate_arr_lin(DoubleMat *X, DoubleMat *y){
	unsigned int i;
    for (i = 0; i < X->axis1; i++){
        X->vals[0][i] = i;
        y->vals[0][i] = i * 2.542 + 7.65342;
    }
}

int main(){
	clock_t start, end;
	double cpu_time_used;

	
	srand(time(NULL));
	unsigned int i, j;
	DoubleMat *X, *y, *theta;
	unsigned int param_num, data_len;
	double alpha = 0.001;
	data_len = 3;
	param_num = 1;
	double **X_vals = alloc_nm(param_num, data_len);
	double **y_vals = alloc_nm(1,data_len);
	double **theta_vals = alloc_nm(param_num + 1, 1);

	X = createMat(param_num, data_len, X_vals);
	y = createMat(1, data_len, y_vals);
	theta = ones(1, param_num + 1);
	populate_arr_lin(X, y);
	printMat(X);
	printMat(y);
	printMat(theta);
	
	DoubleMat *X1 = concatenate(X, ones(1, X->axis1), 0);
	printf("hypLog is: \n");
	DoubleMat *theta_exp = hypLog(X1, theta);
	printMat(theta_exp);
	free(theta_exp);

	printf("hypLin is: \n");
	theta_exp = hypLinear(X1, theta);
	printMat(theta_exp);
	free(theta_exp);

	printf("\n");
	printMat(theta);
	start = clock();
	for (i = 0; i < 20000; i++){
		doIteration(X1, y, theta, alpha);
	}
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("cpu_time_used is: %lf\n", cpu_time_used);
	printMat(theta);
	freeMat(X1);
	freeMat(X);
	freeMat(theta);
	return 0;
}