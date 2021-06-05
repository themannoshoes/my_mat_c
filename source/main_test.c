#include "myMatrix.h"
#include "string.h"

int8_t test_mini_matrix()
{
    int8_t temp_ret_sts;
    Matrix_t mat,inv_mat;
    
    float test_data1[20];
    float test_data2[20];

    memset((void *)test_data1, 0, 16);
    memset((void *)test_data1, 0, 16);

    mat=create_mat(4,4, &temp_ret_sts);
    if(temp_ret_sts != 0){
        //before return, we need to free space used
        if(mat.data != NULL){
            free_mat(&mat);
        }
        return temp_ret_sts;
    }
 
    float data[16] = { 1, 2, 3, 4,
                       2, 1, 4, 5,
                       5, 4, 3, 5,
                       7, 6, 5, 4 };
    set_mat_data(&mat, data, &temp_ret_sts);
    show_mat_to_float("mat", &mat, test_data1);
    inv_mat=inverse_mat(&mat, &temp_ret_sts);
    if(temp_ret_sts != 0){
        //before return, we need to free space used
        if(mat.data != NULL)free_mat(&mat);
        if(inv_mat.data != NULL)free_mat(&inv_mat);

        return temp_ret_sts;
    }

    show_mat_to_float("inv_mat", &inv_mat, test_data2);
    free_mat(&mat);
    free_mat(&inv_mat);
    return 0;
}

int8_t test_mul_cross_matrix()
{
    int8_t temp_ret_sts;
    Matrix_t mat1,mat2,mat_cross;
    float m_temp_hold[3] = {0};


    mat1=create_mat(3,1, &temp_ret_sts);
    if(temp_ret_sts != 0){
        //before return, we need to free space used
        if(mat1.data != NULL)free_mat(&mat1);
        return temp_ret_sts;
    }
    float test_data1[3] = {3,4, 6};
    set_mat_data(&mat1, test_data1, &temp_ret_sts);

 
    mat2=create_mat(3,1, &temp_ret_sts);
    if(temp_ret_sts != 0){
        //before return, we need to free space used
        if(mat2.data != NULL)free_mat(&mat2);
        if(mat1.data != NULL)free_mat(&mat1);
        return temp_ret_sts;
    }
    float test_data2[3] = {1,-6, 0};
    set_mat_data(&mat2, test_data2, &temp_ret_sts);
    
    mat_cross = mult_cross_3Dvector(&mat1, &mat2, &temp_ret_sts);
    if(temp_ret_sts != 0){
        //before return, we need to free space used
        if(mat2.data != NULL)free_mat(&mat2);
        if(mat1.data != NULL)free_mat(&mat1);
        if(mat_cross.data != NULL)free_mat(&mat_cross);
        return temp_ret_sts;
    }
    show_mat_to_float("mat", &mat_cross, m_temp_hold);
    
    
    free_mat(&mat1);
    free_mat(&mat2);
    free_mat(&mat_cross);
    return 0;
}
