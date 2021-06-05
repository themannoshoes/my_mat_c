#ifndef __MY_MATRIX_H
#define __MY_MATRIX_H

#include <stdint.h>

#define  EPSILON  1e-10

#define  ROW_COL_MAX_NUM     12

#define  ERROR_UNSPPT_ROW_COL_NUM  -2
#define  ERROR_MALLOC_FAILED       -3
#define  ERROR_MISMATCH_MATRIXS    -4
#define  ERROR_NULL_INPUT_MATRIX   -5
#define  ERROR_DET_VALUE_ZERO      -6


typedef struct
{
    uint16_t row;
    uint16_t column;
    float **data;
}Matrix_t;



    void init_matrix(Matrix_t *mat);
Matrix_t create_mat(uint16_t row, uint16_t column, int8_t *ret_sts);//create a matrix
Matrix_t add_mat(const Matrix_t* mat1, const Matrix_t* mat2, int8_t *ret_sts);//mat1+mat2;
Matrix_t sub_mat(const Matrix_t* mat1, const Matrix_t* mat2, int8_t *ret_sts);//mat1-mat2;
Matrix_t transpose_mat(const Matrix_t* mat, int8_t *ret_sts);//mat'
Matrix_t scale_mat(const Matrix_t* mat, const float scaler, int8_t *ret_sts);//scaler*Mat
Matrix_t mult_mat(const Matrix_t *mat1, const Matrix_t* mat2, int8_t *ret_sts);//mat1*mat2
Matrix_t assemble_2mat_columns(const Matrix_t *mat1, const Matrix_t *mat2, int8_t *ret_sts); //assemble 2 matrix into 1 matrix
Matrix_t assemble_2mat_rows(const Matrix_t *mat1, const Matrix_t *mat2, int8_t *ret_sts);
Matrix_t assemble_mats_columns(Matrix_t *Matrix_array, uint8_t matrixs_amount , int8_t *ret_sts);
Matrix_t assemble_mats_rows(Matrix_t *Matrix_array, uint8_t matrixs_amount , int8_t *ret_sts);
   float det_mat(Matrix_t *m, int8_t *ret_sts);//get matrix's derterminent value
Matrix_t inverse_mat(Matrix_t* m, int8_t *ret_sts);//get inverse matrix
Matrix_t eye(uint16_t n, int8_t *ret_sts);//generate I(nxn) matrix
Matrix_t diag_mat(uint16_t n, float* diag, int8_t *ret_sts);//generate diag matrix which is nxn matrix
Matrix_t copy_mat(const Matrix_t *mat, int8_t *ret_sts);//copy a matrix
Matrix_t copy_mat_columns(const Matrix_t *mat, uint16_t column_table, int8_t *ret_sts); //copy some columns of origin matrix and then create new matrix
Matrix_t copy_mat_rows(const Matrix_t *mat, uint16_t row_table, int8_t *ret_sts);
Matrix_t mult_cross_3Dvector(const Matrix_t *mat1, const Matrix_t *mat2, int8_t *ret_sts);
    void copy_mat_data(const Matrix_t* mat, Matrix_t* copy, int8_t *ret_sts);//copy matrix's data to another matrix
    void free_mat(Matrix_t *mat);//free a matrix
    void show_mat_to_float(const char* name,const Matrix_t *mat, float *array_1D);   //show the matrix float
    void show_mat_to_double(const char* name,const Matrix_t *mat, double *array_1D); //show the matrix double
    void set_mat_data(Matrix_t* mat,const float *data, int8_t *ret_sts);//set data to matrix
    void clear_mat(Matrix_t* mat, int8_t * ret_sts);//set all matrix's data to 0
uint16_t get_mat_not_zero_columntable(const Matrix_t *mat, int8_t *ret_sts);
uint16_t get_mat_not_zero_rowtable(const Matrix_t *mat, int8_t *ret_sts);

void swap_row_mat(Matrix_t* mat,uint16_t m,uint16_t n, int8_t *ret_sts);//swap NO.m and NO.n row in mat
void scale_row_mat(Matrix_t* mat, uint16_t m, float scaler, int8_t *ret_sts);//NO.m row in matrix multiply a scaler



int8_t test_mini_matrix(void);
int8_t test_mul_cross_matrix(void);
#endif



