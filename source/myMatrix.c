#include "myMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>



/*
* Description: init the matrix 
*/
void init_matrix(Matrix_t *mat)
{
     mat->column = 0;
     mat->row    = 0;
     mat->data   = NULL;
}


/*
create a matrix with 0
*/
Matrix_t create_mat(uint16_t row, uint16_t column, int8_t *ret_sts)
{
    Matrix_t mat;
    uint16_t i;

    //init the matrix data pointer
    mat.data = NULL;
 
    if (row > ROW_COL_MAX_NUM
      ||column > ROW_COL_MAX_NUM
      ||row == 0
      ||column == 0){
//     printf("error, in create_mat: row <= 0||column<=0\n");
        *ret_sts = ERROR_UNSPPT_ROW_COL_NUM;
        return mat;
    }
    if (row > 0 && column > 0){
        mat.row = row;
        mat.column = column;
        mat.data = (float **)malloc(row*sizeof(float *));//先指针的指针
        if (mat.data == NULL){
        //      printf("error, in create_mat: mat.data==NULL");
            *ret_sts = ERROR_MALLOC_FAILED;
            return mat;
        }

        for (i = 0; i < row; i++){
            *(mat.data + i) = (float *)malloc(column*sizeof(float));//再分配每行的指针
            if (mat.data[i] == NULL){
        //        printf("error, in create_mat: mat.data==NULL");
                *ret_sts = ERROR_MALLOC_FAILED;
                return mat;
            }
        }
        clear_mat(&mat, ret_sts);
    }
    return mat;
}
/*
free a matrix
*/
void free_mat(Matrix_t *mat)
{
    uint16_t i;

     //check input matrix 
    if(mat->data == NULL){
        return;
    }
    for (i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            return;
        }
    }

    for (i = 0; i < mat->row; i++){

        free(mat->data[i]);/*释放行*/
        mat->data[i] = NULL;/* we need to set it NULL after free*/
    }
    free(mat->data);/*释放头指针*/
    mat->data = NULL;/* we need to set it NULL after free*/
}

/*
display the matrix
*/
void show_mat_to_float(const char* name,const Matrix_t *mat, float *array_1D)
{
    uint16_t i, j;

     //check input matrix 
    if(mat->data == NULL){
        return;
    }
    for (i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            return;
        }
    }
 
    for (i = 0; i < mat->row; i++){
        for (j = 0; j < mat->column; j++){
//      printf("%.6f\t", mat->data[i][j]);
        array_1D[i*mat->column+j] = mat->data[i][j];
        }
//     printf("\n");
    }
}
/*
display the matrix
*/
void show_mat_to_double(const char* name,const Matrix_t *mat, double *array_1D)
{
    uint16_t i, j;

     //check input matrix 
    if(mat->data == NULL){
        return;
    }
    for (i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            return;
        }
    }
 
    for (i = 0; i < mat->row; i++){
        for (j = 0; j < mat->column; j++){
//      printf("%.6f\t", mat->data[i][j]);
        array_1D[i*mat->column+j] = mat->data[i][j];
        }
//     printf("\n");
    }
}
/*
set datas to the matrix
*/
void set_mat_data(Matrix_t* mat,const float *data, int8_t *ret_sts)
{
    uint16_t i, j;
 
    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return;
    }
    for (i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return;
        }
    }
    
    for (i = 0; i < mat->row; i++){
        for (j = 0; j < mat->column; j++){
            mat->data[i][j] = data[i*mat->column+j];
        }
    }
}
/*
mat=mat1+mat2
*/
Matrix_t add_mat(const Matrix_t* mat1, const Matrix_t* mat2, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;

    Matrix_t mat;
    uint16_t i, j;
   
    //init the matrix data pointer
    mat.data = NULL;
 
    //check input matrix 
    if(mat1->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat;
    }
    for(i = 0; i < mat1->row; i++){
        if (mat1->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat;
        }
    }
    //check input matrix 
    if(mat2->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat;
    }
    for(i = 0; i < mat2->row; i++){
        if (mat2->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat;
        }
    }

    if (mat1->row != mat2->row){
//     printf("error, in add_mat: mat1->row != mat2->row\n");
       *ret_sts = ERROR_MISMATCH_MATRIXS;
       return mat;
    }
    if (mat1->column != mat2->column){
//     printf("error, in add_mat: mat1->column != mat2->column\n");
       *ret_sts = ERROR_MISMATCH_MATRIXS;
       return mat;
    }

    mat = create_mat(mat1->row, mat1->column, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat;
    }
    
    for (i = 0; i < mat1->row; i++){

        for (j = 0; j < mat1->column; j++)
        mat.data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    }
    return mat;
}
/*
mat=mat1-mat2
*/
Matrix_t sub_mat(const Matrix_t* mat1, const Matrix_t* mat2, int8_t *ret_sts)//mat1-mat2;
{
    int8_t create_ret_sts = 0;
    Matrix_t mat;
    uint16_t i, j;

    //init the matrix data pointer
    mat.data = NULL;
 
    //check input matrix 
    if(mat1->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat;
    }
    for(i = 0; i < mat1->row; i++){
        if (mat1->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat;
        }
    }
    //check input matrix 
    if(mat2->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat;
    }
    for(i = 0; i < mat2->row; i++){
        if (mat2->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat;
        }
    }

    if (mat1->row != mat2->row){
//        printf("error, in sub_mat: mat1->row != mat2->row\n");
        *ret_sts = ERROR_MISMATCH_MATRIXS;
        return mat;
    }
    if (mat1->column != mat2->column){
//        printf("error, in sub_mat: mat1->column != mat2->column\n");
        *ret_sts = ERROR_MISMATCH_MATRIXS;
        return mat;
    }

    mat = create_mat(mat1->row, mat1->column, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat;
    }

    for (i = 0; i < mat1->row; i++){
        for (j = 0; j < mat1->column; j++)
        mat.data[i][j] = mat1->data[i][j] - mat2->data[i][j];
    }
    return mat;
}
/*
transpose the matrix, mat=mat'
*/
Matrix_t transpose_mat(const Matrix_t* mat, int8_t *ret_sts)//mat'
{
    int8_t create_ret_sts = 0;
    Matrix_t mat_T;
    uint16_t i, j;

    //init the matrix data pointer
    mat_T.data = NULL;
 
    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_T;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_T;
        }
    }

    mat_T = create_mat(mat->column, mat->row, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_T;
    }

    for (i = 0; i < mat_T.row; i++)
    {
        for (j = 0; j < mat_T.column; j++)
        {
            mat_T.data[i][j] = mat->data[j][i];
        }
    }
    return mat_T;
}
/*
mat=scaler*mat
every element in the matrix multiplys a scaler
*/
Matrix_t scale_mat(const Matrix_t* mat, const float scaler, int8_t *ret_sts)//scaler*Mat
{
    int8_t create_ret_sts = 0; 
    uint16_t i, j;
    Matrix_t mat1;

    //init the matrix data pointer
    mat1.data = NULL;
 
    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat1;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat1;
        }
    }
 
    mat1 = create_mat(mat->row, mat->column, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat1;
    }

    for (i = 0; i < mat->row; i++){
        for (j = 0; j < mat->column; j++){
            mat1.data[i][j] = mat->data[i][j]*scaler;
        }
    }
    return mat1;
}
/*
set all datas in matrix to zero
*/
void clear_mat(Matrix_t* mat, int8_t * ret_sts)
{
    uint16_t i, j;
 
    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return;
        }
    }
 
    for (i = 0; i < mat->row; i++)
    {
        for (j = 0; j < mat->column; j++)
        {
         mat->data[i][j] = 0;
        }
    }
}
/*
mat=mat1*mat2
*/
Matrix_t mult_mat(const Matrix_t *mat1, const Matrix_t* mat2, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    Matrix_t mat;
    uint16_t i, j, m;

    //init the matrix data pointer
    mat.data = NULL;
 
    //check input matrix 
    if(mat1->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat;
    }
    for(i = 0; i < mat1->row; i++){
        if (mat1->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat;
        }
    }

    //check input matrix 
    if(mat2->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat;
    }
    for(i = 0; i < mat2->row; i++){
        if (mat2->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat;
        }
    }

    if (mat1->column != mat2->row){
//     printf("error,In mult_mat: mat1->column != mat2->row\n");
        *ret_sts = ERROR_MISMATCH_MATRIXS;
        return mat;
    }else{
        mat = create_mat(mat1->row, mat2->column, &create_ret_sts);
        if(create_ret_sts != 0){
            *ret_sts = create_ret_sts;
            return mat;
        }
         clear_mat(&mat, ret_sts);

         for (i = 0; i < mat1->row; i++){
             for (j = 0; j < mat2->column; j++){
                 
                 for (m = 0; m < mat1->column; m++){
                     mat.data[i][j] += mat1->data[i][m] * mat2->data[m][j];
                 }
             }
        }
    }
    return mat;
}
/*
generate a I(nxn) matrix
*/
Matrix_t eye(uint16_t n, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i,j;
    Matrix_t mat;

    //init the matrix data pointer
    mat.data = NULL;
 
    mat = create_mat(n, n, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat;
    }

    if (n <= 0){
//        printf("error, in eye: n<0\n");
        *ret_sts = ERROR_UNSPPT_ROW_COL_NUM;
        return mat;
    }

    for (i = 0; i < n;i++){
        for (j = 0; j < n; j++){
            if (i == j)
             mat.data[i][j] = 1;
            else
             mat.data[i][j] = 0;
        }
    }

    return mat;
}
/*
generate a diagonal matrix with diag[n] as its main diagonal elements
*/
Matrix_t diag_mat(uint16_t n, float* diag, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat;

    //init the matrix data pointer
    mat.data = NULL;
 
    mat = create_mat(n, n, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat;
    }
    if (n <= 0){
//       printf("error: in diag_mat(n<0)\n");
        *ret_sts = ERROR_UNSPPT_ROW_COL_NUM;
        return mat;
    }

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i == j)
             mat.data[i][j] = diag[i];
            else
             mat.data[i][j] = 0;
        }
    }

    return mat;
}
/*
copy a matrix
*/
Matrix_t copy_mat(const Matrix_t *mat, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat_copy;

    //init the matrix data pointer
    mat_copy.data = NULL;

    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_copy;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_copy;
        }
    }

    mat_copy = create_mat(mat->row, mat->column, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_copy;
    }
 
    for (i = 0; i < mat->row; i++){
        for (j = 0; j < mat->column; j++){
            mat_copy.data[i][j] = mat->data[i][j];
        }
    }
    return mat_copy;
}

/*
* Description: copy a matrix columns and create a new matrix
* Input para:  1.column_table: from bit0 to bit15, means column 1 to column 16
               2.....
*/
Matrix_t copy_mat_columns(const Matrix_t *mat, uint16_t column_table, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat_copy;
    uint8_t  copy_column_amount = 0;
    uint8_t copy_column_id[16];
    uint8_t array_cursor = 0;

    //init the matrix data pointer
    mat_copy.data = NULL;

    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_copy;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_copy;
        }
    }

    /*get the column amount we need to copy, and the copy_column_id array*/
    for(i=0;i < mat->column;i++){
        if(column_table & (0x0001 << i)){
            copy_column_amount += 1;
            copy_column_id[array_cursor] = i;
            array_cursor++;
        }
    }

    mat_copy = create_mat(mat->row, copy_column_amount, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_copy;
    }

    for (i = 0; i < mat_copy.row; i++){
        for (j = 0; j < mat_copy.column; j++){
            mat_copy.data[i][j] = mat->data[i][ copy_column_id[j] ];
        }
    }
    return mat_copy;
}

/*
* Description: copy a matrix rows and create a new matrix
* Input para:  1.column_table: from bit0 to bit15, means row 1 to row 16
               2.....
*/
Matrix_t copy_mat_rows(const Matrix_t *mat, uint16_t row_table, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat_copy;
    uint8_t  copy_row_amount = 0;
    uint8_t copy_row_id[16];
    uint8_t array_cursor = 0;
   

    //init the matrix data pointer
    mat_copy.data = NULL;

    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_copy;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_copy;
        }
    }

    /*get the column amount we need to copy, and the copy_column_id array*/
    for(i=0;i < mat->row;i++){
        if(row_table & (0x0001 << i)){
            copy_row_amount += 1;
            copy_row_id[array_cursor] = i;
            array_cursor++;
        }
    }

    mat_copy = create_mat(copy_row_amount , mat->column, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_copy;
    }

    for (i = 0; i < mat_copy.row; i++){
        for (j = 0; j < mat_copy.column; j++){
            mat_copy.data[i][j] = mat->data[ copy_row_id[i] ][j];
        }
    }
    return mat_copy;
}


/*
copy matrix's data to another matrix
*/
void copy_mat_data(const Matrix_t* mat, Matrix_t* copy, int8_t *ret_sts)
{
    uint16_t i, j;

    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return;
        }
    }

    if (mat->row != copy->row || mat->column != copy->column){
//     printf("error, in copy_mat_data: mat->row != copy->row || mat->column != copy->column\n");
       *ret_sts = ERROR_MISMATCH_MATRIXS;
       return;
    }

    for (i = 0; i < mat->row;i++){
        for (j = 0; j < mat->column; j++){
            copy->data[i][j] = mat->data[i][j];
        }
    }
}
/*
get matrix's derterminent value
*/
float det_mat(Matrix_t *m, int8_t *ret_sts)
{
    int8_t copy_ret_sts = 0;
    uint16_t i, j, n, max_row;
    int16_t swap_f;
    float max, k;
    float det=1;
    swap_f = 0;
    Matrix_t mat;

    //init the matrix data pointer
    mat.data = NULL;
 
    //check input matrix 
    if(m->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return 1;
    }
    for(i = 0; i < m->row; i++){
        if (m->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return 1;
        }
    }
 
    if (m->column != m->row){
//     printf("error:In det_mat (m->column != m->row)\n");
       *ret_sts = ERROR_MISMATCH_MATRIXS;
       return det;
    }
    mat = copy_mat(m, &copy_ret_sts);
    if(copy_ret_sts != 0){
        *ret_sts = copy_ret_sts;
        //before return, we need to free space used
        if(mat.data != NULL){
            free_mat(&mat);
        }
        return det;
    }
    
    for (i = 0; i < mat.row-1; i++){
        max = fabs(mat.data[i][i]);
        max_row = i;
        for (j = i + 1; j < mat.row; j++){
            if (max < fabs(mat.data[j][i])){
                max = fabs(mat.data[j][i]);
                max_row = j;
            }
        }
        if (i != max_row){
            swap_row_mat(&mat, i, max_row, ret_sts);
            swap_f++;
        }
        for (j = i + 1; j < mat.row; j++){
            k = -mat.data[j][i]/mat.data[i][i];
            for (n= 0; n < mat.column; n++){
                mat.data[j][n] = mat.data[i][n] * k + mat.data[j][n];
            }
        }
    }
    if (swap_f%2==1)swap_f = -1;
    else swap_f = 1;
    det = 1;
    for (i = 0; i < mat.column; i++)
    det *= mat.data[i][i];
    det *= swap_f;

    free_mat(&mat);
    return det;
}
/*
get inverse matrix
use main column element of Gauss-Jordan algrithm: A|I  --->  I|A^(-1)
*/
Matrix_t inverse_mat(Matrix_t* m, int8_t *ret_sts)
{
    int8_t temp_ret_sts = 0;
    int i, j, n,max_row;
    int swap_f = 0;
    float max,k;
    Matrix_t inv_mat, mat;

    //init the matrix data pointer
    inv_mat.data = NULL;
    mat.data = NULL;

    //check input matrix 
    if(m->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return inv_mat;
    }
    for(i = 0; i < m->row; i++){
        if (m->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return inv_mat;
        }
    }
 
    inv_mat = eye(m->row, &temp_ret_sts);
    if(temp_ret_sts != 0){
        *ret_sts = temp_ret_sts;
        return inv_mat;
    }

    if (det_mat(m, &temp_ret_sts) == 0){
//       printf("error: In inverse_mat(det_mat(mat) == 0)\n");
        *ret_sts = ERROR_DET_VALUE_ZERO;
        return inv_mat;
    }
    if(temp_ret_sts != 0){
        *ret_sts = temp_ret_sts;
        return inv_mat;
    }

    mat = copy_mat(m, &temp_ret_sts);
    if(temp_ret_sts != 0){
        *ret_sts = temp_ret_sts;
        if(mat.data != NULL) free_mat(&mat);

        return inv_mat;
    }

    for (i = 0; i < mat.row - 1; i++){
        max = fabs(mat.data[i][i]);
        max_row = i;
        for (j = i + 1; j < mat.row; j++){
            if (max < fabs(mat.data[j][i])){
                max = fabs(mat.data[j][i]);
                max_row = j;
            }
        }
        if (i != max_row){
            swap_row_mat(&mat, i, max_row, ret_sts);
            swap_row_mat(&inv_mat, i, max_row, ret_sts);
            swap_f++;
        }
        for (j = i + 1; j < mat.row; j++){
            k = -mat.data[j][i] / mat.data[i][i];
            for (n = 0; n < mat.column; n++){
                mat.data[j][n] = mat.data[i][n] * k + mat.data[j][n];
                inv_mat.data[j][n] = inv_mat.data[i][n] * k + inv_mat.data[j][n];
            }
        }
    }

    for (i = 0; i < mat.row; i++){
        k = 1/mat.data[i][i];
        scale_row_mat(&mat,i, k, ret_sts);
        scale_row_mat(&inv_mat, i, k, ret_sts);
    }
    for (i = mat.row-1; i>0; i--){
        for (j = i - 1; j >=0; j--)
        {
            k = -mat.data[j][i] / mat.data[i][i];
            for (n = 0; n < mat.column; n++){
                mat.data[j][n] = mat.data[j][n] + k*mat.data[i][n];
                inv_mat.data[j][n] = inv_mat.data[j][n] + k*inv_mat.data[i][n];
            }
        }
    }

    free_mat(&mat);
    return inv_mat;
}
/*swap NO.m and NO.n row in mat*/
void swap_row_mat(Matrix_t* mat,uint16_t m,uint16_t n, int8_t *ret_sts)
{
    float temp;
    uint16_t i;
 
    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return;
        }
    }

    for (i = 0; i < mat->column; i++)
    {
        temp = mat->data[m][i];
        mat->data[m][i] = mat->data[n][i];
        mat->data[n][i] = temp;
    }
}
/*
NO.m row in matrix multiply a scaler
*/
void scale_row_mat(Matrix_t* mat, uint16_t m, float scaler, int8_t *ret_sts)
{
    uint16_t i;

    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return;
        }
    }

    for (i = 0; i < mat->column; i++)
     mat->data[m][i] *= scaler;
}

/*
assemble a matrix
*/
Matrix_t assemble_2mat_columns(const Matrix_t *mat1, const Matrix_t *mat2, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat_assemble;

    //init the matrix data pointer
    mat_assemble.data = NULL;

    //check input matrix 
    if(mat1->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_assemble;
    }
    for(i = 0; i < mat1->row; i++){
        if (mat1->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_assemble;
        }
    }

    //check input matrix 
    if(mat2->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_assemble;
    }
    for(i = 0; i < mat2->row; i++){
        if (mat2->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_assemble;
        }
    }
    
    //check if row_num equal
    if(mat1->row != mat2->row){
        *ret_sts = ERROR_MISMATCH_MATRIXS;
        return mat_assemble;
    }

    mat_assemble = create_mat( (mat1->row), (mat1->column + mat2->column), &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_assemble;
    }
 
    for (i = 0; i < mat_assemble.row; i++){
        for (j = 0; j < mat_assemble.column; j++){
            if(j < mat1->column){
                mat_assemble.data[i][j] = mat1->data[i][j];
            }else{
                mat_assemble.data[i][j] = mat2->data[i][j- mat1->column];
            }
        }
    }
    return mat_assemble;
}

/*
assemble  matrixs columns to one matrix
*/
Matrix_t assemble_mats_columns(Matrix_t *Matrix_array, uint8_t matrixs_amount , int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat_assemble;
 
    //init the matrix data pointer
    mat_assemble.data = NULL;

    
    //check input matrix 
    for(i = 0;i < matrixs_amount;i++){

        if(Matrix_array[i].data == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_assemble;
        }
        for(j = 0; j < Matrix_array[i].row; j++){
            if (Matrix_array[i].data[j] == NULL){
                *ret_sts = ERROR_NULL_INPUT_MATRIX;
                return mat_assemble;
            }
        }
    }

    
    //check if row_num equal
    uint16_t record_row_num = Matrix_array[0].row;
    for(i = 0;i < matrixs_amount;i++){
         if(Matrix_array[i].row != record_row_num){
             *ret_sts = ERROR_MISMATCH_MATRIXS;
             return mat_assemble;
         }
    }

    uint16_t total_columns_num = 0;
    for(i = 0;i < matrixs_amount;i++){
        total_columns_num += Matrix_array[i].column;
    }
    
    mat_assemble = create_mat(record_row_num, total_columns_num, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_assemble;
    }
 
    uint16_t k = 0;
    uint16_t ruler_cols_num = Matrix_array[0].column;
    uint16_t ruler_cols_bottom_cursor = 0;


    for (j = 0; j < mat_assemble.column; j++){
       for (i = 0; i < mat_assemble.row; i++){
                while(j >= ruler_cols_num){
                    if(k >= matrixs_amount -1)return mat_assemble;
                    k++;
                    ruler_cols_bottom_cursor = ruler_cols_num;
                    ruler_cols_num += Matrix_array[k].column;
                }

                if(j < ruler_cols_num
                && j >= ruler_cols_bottom_cursor){
                    mat_assemble.data[i][j] = Matrix_array[k].data[i][j - ruler_cols_bottom_cursor];
                }
        }
    }
    return mat_assemble;
}


/*
assemble two  matrixs to one matrix
*/
Matrix_t assemble_2mat_rows(const Matrix_t *mat1, const Matrix_t *mat2, int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat_assemble;

    //init the matrix data pointer
    mat_assemble.data = NULL;

    //check input matrix 
    if(mat1->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_assemble;
    }
    for(i = 0; i < mat1->row; i++){
        if (mat1->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_assemble;
        }
    }

    //check input matrix 
    if(mat2->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_assemble;
    }
    for(i = 0; i < mat2->row; i++){
        if (mat2->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_assemble;
        }
    }
    
    //check if row_num equal
    if(mat1->column != mat2->column){
        *ret_sts = ERROR_MISMATCH_MATRIXS;
        return mat_assemble;
    }

    mat_assemble = create_mat( (mat1->row + mat2->row), mat1->column, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_assemble;
    }
 
    for (i = 0; i < mat_assemble.row; i++){
        for (j = 0; j < mat_assemble.column; j++){
            if(i < mat1->row){
                mat_assemble.data[i][j] = mat1->data[i][j];
            }else{
                mat_assemble.data[i][j] = mat2->data[i - mat1->row][j];
            }
        }
    }
    return mat_assemble;
}

/*
assemble  matrixs rows to one matrix
*/
Matrix_t assemble_mats_rows(Matrix_t *Matrix_array, uint8_t matrixs_amount , int8_t *ret_sts)
{
    int8_t create_ret_sts = 0;
    uint16_t i, j;
    Matrix_t mat_assemble;
 
    //init the matrix data pointer
    mat_assemble.data = NULL;

    
    //check input matrix 
    for(i = 0;i < matrixs_amount;i++){

        if(Matrix_array[i].data == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_assemble;
        }
        for(j = 0; j < Matrix_array[i].row; j++){
            if (Matrix_array[i].data[j] == NULL){
                *ret_sts = ERROR_NULL_INPUT_MATRIX;
                return mat_assemble;
            }
        }
    }

    
    //check if row_num equal
    uint16_t record_col_num = Matrix_array[0].column;
    for(i = 0;i < matrixs_amount;i++){
         if(Matrix_array[i].column != record_col_num){
             *ret_sts = ERROR_MISMATCH_MATRIXS;
             return mat_assemble;
         }
    }

    uint16_t total_rows_num = 0;
    for(i = 0;i < matrixs_amount;i++){
        total_rows_num += Matrix_array[i].row;
    }
    
    mat_assemble = create_mat( total_rows_num, record_col_num, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_assemble;
    }
 
    uint16_t k = 0;
    uint16_t ruler_rows_num = Matrix_array[0].row;;
    uint16_t ruler_rows_bottom_cursor = 0;


    for (i = 0; i < mat_assemble.row; i++){
        for (j = 0; j < mat_assemble.column; j++){
                while(i >= ruler_rows_num){
                    if(k >= matrixs_amount -1)return mat_assemble;
                    k++;
                    ruler_rows_bottom_cursor = ruler_rows_num;
                    ruler_rows_num += Matrix_array[k].row;
                }

                if(i < ruler_rows_num
                && i >= ruler_rows_bottom_cursor){
                    mat_assemble.data[i][j] = Matrix_array[k].data[i - ruler_rows_bottom_cursor][j];
                }
        }
    }
    return mat_assemble;
}



/*
 get_mat_not_zero_columntable
*/
uint16_t get_mat_not_zero_columntable(const Matrix_t *mat, int8_t *ret_sts)
{
    uint16_t i, j;
    uint16_t ret_columntable = 0x0000;
    uint16_t zero_row_amount = 0;

    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return ret_columntable;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return ret_columntable;
        }
    }

    //init the ret_columntable 
    for(i = 0;i < mat->column;i++){
        ret_columntable |= (0x01 << i);
    }
    
    //get the column table here
    for (j = 0; j < mat->column; j++){
        if(mat->data[0][j] < (float)EPSILON){    //in this case, we think it is zero
            zero_row_amount = 0;
            for(i=0;i < mat->row;i++){
                if(mat->data[i][j] < (float)EPSILON){
                    zero_row_amount ++;
                }
            }
            //if the whole column is zero
            if(zero_row_amount == mat->row){
                ret_columntable &= ~(0x01 << j);
            }
        }
    }
    return ret_columntable;
}

/*
 get_mat_not_zero_rowtable
*/
uint16_t get_mat_not_zero_rowtable(const Matrix_t *mat, int8_t *ret_sts)
{
    uint16_t i, j;
    uint16_t ret_rowtable = 0x0000;
    uint16_t zero_column_amount = 0;

    //check input matrix 
    if(mat->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return ret_rowtable;
    }
    for(i = 0; i < mat->row; i++){
        if (mat->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return ret_rowtable;
        }
    }

    //init the ret_columntable 
    for(i = 0;i < mat->row;i++){
        ret_rowtable |= (0x01 << i);
    }
    
    //get the column table here
    for (j = 0; j < mat->row; j++){
        if(mat->data[j][0] < (float)EPSILON){    //in this case, we think it is zero
            zero_column_amount = 0;
            for(i=0;i < mat->column;i++){
                if(mat->data[j][i] < (float)EPSILON){
                    zero_column_amount ++;
                }
            }
            //if the whole column is zero
            if(zero_column_amount == mat->column){
                ret_rowtable &= ~(0x01 << j);
            }
        }
    }
    return ret_rowtable;
}


/*
description: vector multiply_cross
*/
Matrix_t mult_cross_3Dvector(const Matrix_t *mat1, const Matrix_t *mat2, int8_t *ret_sts){

    uint16_t i;
    float mat_cross[3];
    int8_t create_ret_sts = 0;
    Matrix_t mat_mult_cross;
    float array1[3], array2[3];
 
    //init the matrix data pointer
    mat_mult_cross.data = NULL;

    //check input matrix 
    if(mat1->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_mult_cross;
    }
    for(i = 0; i < mat1->row; i++){
        if (mat1->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_mult_cross;
        }
    }
    
    if(mat2->data == NULL){
        *ret_sts = ERROR_NULL_INPUT_MATRIX;
        return mat_mult_cross;
    }
    for(i = 0; i < mat2->row; i++){
        if (mat2->data[i] == NULL){
            *ret_sts = ERROR_NULL_INPUT_MATRIX;
            return mat_mult_cross;
        }
    }
    
    //check if input matrixs are vectors
    if(mat1->column != 1
    || mat2->column != 1
    || mat1->row != mat2->row
    || mat1->row != 3
    || mat2->row != 3){
        *ret_sts = ERROR_UNSPPT_ROW_COL_NUM;
        return mat_mult_cross;
    }
    
    mat_mult_cross = create_mat( 3, 1, &create_ret_sts);
    if(create_ret_sts != 0){
        *ret_sts = create_ret_sts;
        return mat_mult_cross;
    }
    
    show_mat_to_float("mat", mat1, array1);
    show_mat_to_float("mat", mat2, array2);
    
    mat_cross[0] = array1[1]*array2[2] - array1[2]*array2[1];
    mat_cross[1] = array1[2]*array2[0] - array1[0]*array2[2];
    mat_cross[2] = array1[0]*array2[1] - array1[1]*array2[0];
    set_mat_data(&mat_mult_cross, mat_cross, &create_ret_sts);
    return mat_mult_cross;

}



