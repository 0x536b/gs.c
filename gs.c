/*
Author : Suresh (0x536b)
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>

// #define DEBUG_IDX 1579995
// #define DEBUG_IDX 1316
// #define DEBUG_IDX 285983
// #define DEBUG_IDX 0
#define DEBUG_IDX 9999999

#define MAX_LINE_LENGTH 500
#define SCALING_MOD 1

#define TEST_IMG_ID 10
#define TEST_IMG_SCALE 4

// #define IMAGE_W 1920
// #define IMAGE_H 1080

#define ZFAR 100.0
#define ZNEAR 0.01

#define BLOCK_X 16
#define BLOCK_Y 16
#define BLOCK_SIZE (BLOCK_X * BLOCK_Y)


#define MAX_SH_DEGREE 3
#define CHANNELS 3
// #define MAX_SH_COEFF 16


#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

typedef struct
{
    float x, y, z;       // 3
    float nx, ny, nz;    // 3
    float sh_coeff[48];  // 48
    float opacity;       // 1
    float scale[3];      // 3
    float quaternion[4]; // 4
} Gaussian;              // 62 or 59

// typedef struct
// {
//     int id;
//     char name;
//     int width;
//     int height;
//     float pos[3];
//     float rotation[9];
//     float focal_x;
//     float focal_y;
//     float fov_x;
//     float fov_y;
// } Viewport;

// typedef struct
// {
//     uint64_t g_key;
//     int g_idx;
// } g_kv; // gaussian_keys_values


const float SH_C0 = 0.28209479177387814;
const float SH_C1 = 0.4886025119029199;
const float SH_C2[5] = {
	1.0925484305920792,
	-1.0925484305920792,
	0.31539156525252005,
	-1.0925484305920792,
	0.5462742152960396
};
const float SH_C3[7] = {
	-0.5900435899266435,
	2.890611442640554,
	-0.4570457994644658,
	0.3731763325901154,
	-0.4570457994644658,
	1.445305721320277,
	-0.5900435899266435
};


float max(float a, float b) {
    if (a >= b) {
        return a;
    } else {
        return b;
    }
}

float min(float a, float b) {
    if (a <= b) {
        return a;
    } else {
        return b;
    }
}

float ndc2Pix(float v, int S) {
    // printf("ANS: %f\n", (v + 1.0));
    // printf("ANS: %f\n", ((v + 1.0)*S));
    // printf("ANS: %f\n", ((v + 1.0)*S) - 1.0);
    // printf("ANS: %f\n", (((v + 1.0)*S) - 1.0)*0.5 );
    float tmp = ((v + 1.0) * S - 1.0) * 0.5;
    // printf("TMP: %f\n", tmp);
    // printf("--------\n");
    return tmp;
}

// matrix multiply of two 3x3 matrices
// [(0,1,2),
//  (3,4,5),
//  (6,7,8)]
// [(00, 01, 02),
//  (10, 11, 12),
//  (20, 21, 22)]
void matmul3x3(float *a, float *b, float *result)
{
    result[0] = (a[0] * b[0]) + (a[1] * b[3]) + (a[2] * b[6]);
    result[1] = (a[0] * b[1]) + (a[1] * b[4]) + (a[2] * b[7]);
    result[2] = (a[0] * b[2]) + (a[1] * b[5]) + (a[2] * b[8]);
    result[3] = (a[3] * b[0]) + (a[4] * b[3]) + (a[5] * b[6]);
    result[4] = (a[3] * b[1]) + (a[4] * b[4]) + (a[5] * b[7]);
    result[5] = (a[3] * b[2]) + (a[4] * b[5]) + (a[5] * b[8]);
    result[6] = (a[6] * b[0]) + (a[7] * b[3]) + (a[8] * b[6]);
    result[7] = (a[6] * b[1]) + (a[7] * b[4]) + (a[8] * b[7]);
    result[8] = (a[6] * b[2]) + (a[7] * b[5]) + (a[8] * b[8]);
}

void matmul(float *a, int a_rows, int a_cols, float *b, int b_rows, int b_cols, float *result)
{
    if (a_cols != b_rows) {
        printf("Error: Matrices cannot be multiplied\n");
        exit(1);
    }

    int i,j, k;
    for (i = 0; i < a_rows; i++) {
        for (j = 0; j < b_cols; j++) {
            for (k = 0; k < a_cols; k++) {
                // printf("res[%d] += a[%d] * b[%d]\n", (i*a_cols + j), (i*a_cols + k), (k*b_cols + j) );
                result[(i*a_cols + j)] += a[(i*a_cols + k)] * b[(k*b_cols + j)];
            }
            // printf("------\n");
        }
    }
}

// multiply a 1x3 and 3x3 matrix. Result is 1x3
void matmul1x3(float *a, float *b, float *result) {
    result[0] = (a[0] * b[0]) + (a[1] * b[3]) + (a[2] * b[6]);
    result[1] = (a[0] * b[1]) + (a[1] * b[4]) + (a[2] * b[7]);
    result[2] = (a[0] * b[2]) + (a[1] * b[5]) + (a[2] * b[8]);
}

// transpose a 3x3 matrix
void transpose3x3(float *a, float *result)
{
    result[0] = a[0];
    result[1] = a[3];
    result[2] = a[6];
    result[3] = a[1];
    result[4] = a[4];
    result[5] = a[7];
    result[6] = a[2];
    result[7] = a[5];
    result[8] = a[8];
}

// transpose a 4x4 matrix
void transpose4x4(float *a, float *result)
{
    result[0] = a[0];
    result[1] = a[4];
    result[2] = a[8];
    result[3] = a[12];
    result[4] = a[1];
    result[5] = a[5];
    result[6] = a[9];
    result[7] = a[13];
    result[8] = a[2];
    result[9] = a[6];
    result[10] = a[10];
    result[11] = a[14];
    result[12] = a[3];
    result[13] = a[7];
    result[14] = a[11];
    result[15] = a[15];
}

// convert fov value to focal value
float fov2focal(float fov, int length) {
    float focal = 0.0;
    focal = length / (2 * tan((fov/2)));
    return focal;
}

// convert focal value to fov
float focal2fov(float focal, int length) {
    float fov = 0.0;
    fov = 2 * atan((length/(2 * focal)));
    return fov;
}

// find determinant of 3x3 matrix
void determinant3x3(float *a, float *ans)
{
    float det = 0.0;
    det = (a[0] * ((a[4] * a[8]) - (a[5] * a[7]))) - (a[3] * ((a[1] * a[8]) - (a[2] * a[7]))) + (a[6] * ((a[1] * a[5]) - (a[2] * a[4])));
    ans[0] = det;
}

// find cofactor of 3x3 matrix
void cofactor3x3(float *a, float *ans)
{
    ans[0] = ((a[4] * a[8]) - (a[5] * a[7]));
    ans[1] = -1 * ((a[3] * a[8]) - (a[6] * a[5]));
    ans[2] = ((a[3] * a[7]) - (a[6] * a[4]));
    ans[3] = -1 * ((a[1] * a[8]) - (a[2] * a[7]));
    ans[4] = ((a[0] * a[8]) - (a[2] * a[6]));
    ans[5] = -1 * ((a[0] * a[7]) - (a[1] * a[6]));
    ans[6] = ((a[1] * a[5]) - (a[2] * a[4]));
    ans[7] = -1 * ((a[0] * a[5]) - (a[3] * a[2]));
    ans[8] = ((a[0] * a[4]) - (a[1] * a[3]));
}

// multiply a constant to all elements
void multiply_const_3x3(float *a, float *constant, float *result)
{
    result[0] = a[0] * constant[0];
    result[1] = a[1] * constant[0];
    result[2] = a[2] * constant[0];
    result[3] = a[3] * constant[0];
    result[4] = a[4] * constant[0];
    result[5] = a[5] * constant[0];
    result[6] = a[6] * constant[0];
    result[7] = a[7] * constant[0];
    result[8] = a[8] * constant[0];
}

// inverse of 3x3 matrix
// using the adjoint method
void inverse3x3(float *a, float *result) {
    float *det = NULL;
    det = (float *)calloc(1, sizeof(float));
    float *cofactor = NULL;
    cofactor = (float *)calloc(9, sizeof(float));
    float *adj = NULL;
    adj = (float *)calloc(9, sizeof(float));
    float *one_by_det = NULL;
    one_by_det = (float *)calloc(1, sizeof(float));
    if ((det == NULL) || (cofactor == NULL) || (adj == NULL) || (one_by_det == NULL)) {
        printf("inv: Mem. Alloc. Failed");
        exit(1);
    }
    
    determinant3x3(a, det);
    if (det[0] != 0.0) {
        cofactor3x3(a, cofactor); // find the cofactors
        transpose3x3(cofactor, adj); // transpose to get adjoint
        one_by_det[0] = 1.0 / det[0];
        multiply_const_3x3(adj, one_by_det, result); // multiply with 1/det to get inverse
    } else {
        printf("inv: Determinant is zero\n");
        exit(1);
    }

    free(det);
    free(cofactor);
    free(adj);
}

void inversemat(float *matrix, int size, float *result) {
    int i, j, k, n;
    n = size;

    float *L = NULL;
    L = (float *)malloc(n * n * sizeof(float));
    float *I = NULL;
    I = (float *)malloc(n * n * sizeof(float));
    float *a = NULL;
    a = (float *)malloc(n * n * sizeof(float));
    float *Z = NULL;
    Z = (float *)malloc(n * n * sizeof(float));
    if ((L == NULL) || (I == NULL) || (a == NULL) || (Z == NULL)) {
        printf("inv2: Mem. alloc. failed\n");
        exit(1);
    }

    // init L and I to identity matrix
    // and also copy the original matrix
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                L[i * n + j] = 1;
                I[i * n + j] = 1;
            } else {
                L[i * n + j] = 0;
                I[i * n + j] = 0;
            }
            a[i * n + j] = matrix[i * n + j];
        }
    }

    // start LU decomposition
    // here a = U
    int norm, row, col;
    float multiplier;
    for (norm = 0; norm < n - 1; norm++) {
        for (row = norm + 1; row < n; row++) {
            // idx1 = row * n + norm;
            // idx2 = norm * n + norm;
            multiplier = a[row * n + norm] / a[norm * n + norm];
            for (col = norm; col < n; col++) {
                // idx3 = row * n + col;
                // idx4 = norm * n + col;
                a[row * n + col] -= a[norm * n + col] * multiplier;
            }
            L[row * n + norm] = multiplier;
        }
    }

    // LUz = B, where B = A^-1
    // Forward elimination
    for (k = 0; k < n; k++) {
        for (row = 0; row < n ; row++) {
            Z[row*n + k] = I[row*n + k];
            for (col = 0; col < row; col++) {
                Z[row*n + k] -=  L[row*n + col] * Z[col*n + k];
            }
            Z[row*n + k] /= L[row*n + row];
        }
    }

    // Backward elimination
    for (k = 0; k < n; k++) {
        for (row = n-1; row >= 0 ; row--) {
            result[row*n + k] = Z[row*n + k];
            for (col = n-1; col > row; col--) {
                result[row*n + k] -=  a[row*n + col] * result[col*n + k];
            }
            result[row*n + k] /= a[row*n + row];
        }
    }



    free(L);
    free(I);
    free(a);
    free(Z);
}

void transformPoint4x4(float*p, float* matrix, float* result){
    result[0] = matrix[0] * p[0] + matrix[4] * p[1] + matrix[8] * p[2] + matrix[12];
    result[1] = matrix[1] * p[0] + matrix[5] * p[1] + matrix[9] * p[2] + matrix[13];
    result[2] = matrix[2] * p[0] + matrix[6] * p[1] + matrix[10] * p[2] + matrix[14];
    result[3] = matrix[3] * p[0] + matrix[7] * p[1] + matrix[11] * p[2] + matrix[15];
}

void transformPoint4x3(float*p, float* matrix, float* result) {
    	result[0] = matrix[0] * p[0] + matrix[4] * p[1] + matrix[8] * p[2] + matrix[12];
		result[1] = matrix[1] * p[0] + matrix[5] * p[1] + matrix[9] * p[2] + matrix[13];
		result[2] = matrix[2] * p[0] + matrix[6] * p[1] + matrix[10] * p[2] + matrix[14];
}

void InclusiveSum(int *a, int *sum, int size) {
    float tmp = 0.0;
    for (int i = 0; i < size; i++) {
        tmp = tmp + a[i];
        sum[i] = tmp; 
    }
}

uint32_t getHigherMsb(uint32_t n)
{
	uint32_t msb = sizeof(n) * 4;
	uint32_t step = msb;
	while (step > 1)
	{
		step /= 2;
		if (n >> msb)
			msb += step;
		else
			msb -= step;
	}
	if (n >> msb)
		msb++;
	return msb;
}

// a - array to sort
// n - total elemetns in array
// bit - which bit to sort
void countingSortPairs(uint64_t *keys, int *values, int n, int bit) {
    uint64_t *out_keys = NULL;
    unsigned int *out_values = NULL;
    out_keys = (uint64_t*)malloc(n * sizeof(uint64_t));
    out_values = (unsigned int*)malloc(n * sizeof(unsigned int));

    if ((out_keys == NULL) || (out_values == NULL)) {
        printf("Mem alloc failed, countingsort!\n");
        exit(1);
    }

    int count[2] = {0, 0}; // count[0] - count of zero bits

    int i;
    for (i = 0; i < n; i++) { 
        int bitvalue = (keys[i] >> bit) & 1; // do an AND operation on the bit
        count[bitvalue]++;
    }
    // printf("c-0: %d\n", count[0]);
    // printf("c-1: %d\n", count[1]);

    // following this video https://www.youtube.com/watch?v=OKd534EWcdk
    // doing right shift of values manually
    count[1] = count[0];
    count[0] = 0;

    for (i = 0; i < n; i++) { 
        int bitvalue = (keys[i] >> bit) & 1; // do an AND operation on the bit
        out_keys[count[bitvalue]] = keys[i];
        out_values[count[bitvalue]] = values[i];
        count[bitvalue]++;
    }

    for (i = 0; i < n; i++) {
        keys[i] = out_keys[i];
        values[i] = out_values[i];
    }

    free(out_keys);
    free(out_values);
}

void radixSortPairs(uint64_t *keys, int *values, int n, int highestBit) {
    for (int i = 0; i < highestBit; i++) {
        countingSortPairs(keys, values, n, i);
    }
}

void uint64_to_bin(uint64_t number, char* bin) {
    uint64_t mask = 1 ;
    mask = mask << 63;
    for (int i = 0; i < 64; i++) {
        bin[i] = (number & mask) ? '1' : '0';
        mask = mask >> 1;
    }
    bin[64] = '\0';
}

void printProgress(int idx, int total, double val1) {
    float progress = (float)(idx) / (float)(total);
    int barwidth = 50;
    int pos = barwidth * progress;
    // printf("pos: %d\n", pos);

    printf("[");
    for (int i = 0; i < barwidth; i++) {
        if(i < pos) printf ("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %d%% (%d/%d) [%lf]\r", (int)(progress*100), idx, total, val1);
    fflush(stdout);
}

int main()
{
    clock_t tik1, tik2, tik3, tik4, tik5, tik6, tik7, tik8, tik9, tik10, tok;
    int i, j, k, l; // for for-loop use only
    tik1 = clock();

    FILE *file;
    file = fopen("test.ply", "r");

    if (file == NULL)
    {
        printf("File does not exist\n");
        return 1;
    }

    // Read the header
    char line[MAX_LINE_LENGTH];
    int num_gaussians = 0;

    while (fgets(line, sizeof(line), file) != NULL)
    {
        line[strcspn(line, "\n")] = 0; // put a NULL where the \n is
        if (strstr(line, "element vertex") != NULL)
        {
            sscanf(line, "element vertex %d", &num_gaussians);
            printf("N - %d\n", num_gaussians);
        }
        else if (strstr(line, "end_header") != NULL)
        {
            break;
        }
    }

    Gaussian *gaussians = NULL;
    gaussians = (Gaussian *)malloc(num_gaussians * sizeof(Gaussian));
    if (gaussians == NULL)
    {
        printf("Mem. alloc. failed\n");
        return 1;
    }

    fread(&gaussians[0], sizeof(Gaussian), num_gaussians, file); // read (num_gaussians * sizeof(gaussian struct)) bytes from file and store it at index 0

    //=== Pre process
    //=== Convert scales to exponent of scales
    //=== why?: Sec5.1: We use a sigmoid activation function for 𝛼 to constrain it in
    // the [0 − 1) range and obtain smooth gradients, and an exponential
    // activation function for the scale of the covariance for similar reassons.
    // To make all values positive
    // https://github.com/graphdeco-inria/gaussian-splatting/blob/8a70a8cd6f0d9c0a14f564844ead2d1147d5a7ac/scene/gaussian_model.py#L33
    for (i = 0; i < num_gaussians; i++)
    {
        for (j = 0; j < 3; j++)
        {
            gaussians[i].scale[j] = exp(gaussians[i].scale[j]);
        }
    }
    //== Implement scaling factor
    for (i = 0; i < num_gaussians; i++)
    {
        for (j = 0; j < 3; j++)
        {
            gaussians[i].scale[j] = gaussians[i].scale[j] * SCALING_MOD;
        }
    }

    //=== Normalize rotation matrix
    // https://github.com/graphdeco-inria/gaussian-splatting/blob/8a70a8cd6f0d9c0a14f564844ead2d1147d5a7ac/scene/gaussian_model.py#L41
    //=== Why ? :
    // L2 norm row-wise. I.e. individual rotation vector is normalized by itself and not the entire dataset
    for (i = 0; i < num_gaussians; i++)
    {
        float norm = 0.0;
        float sum = 0.0;
        // Calc L2 norm
        for (j = 0; j < 4; j++)
        {
            sum += pow(gaussians[i].quaternion[j], 2);
        }
        norm = sqrt(sum);
        if (norm > 0)
        {
            for (j = 0; j < 4; j++)
            {
                gaussians[i].quaternion[j] = gaussians[i].quaternion[j] / norm;
            }
        }
    }

    //=== calc sigmoid of opacities
    // sigmoid(x) = 1 / (1 - e^-x)
    // sigmoid(x) = 1 / (1 - exp(-x))
    // To map all values to [0,1)
    // Usefull because you want the weights of weighted sum to be positive and a probability of 0-1
    // https://github.com/graphdeco-inria/gaussian-splatting/blob/8a70a8cd6f0d9c0a14f564844ead2d1147d5a7ac/scene/gaussian_model.py#L38
    for (i = 0; i < num_gaussians; i++)
    {
        gaussians[i].opacity = (1 / (1 + exp(-1 * gaussians[i].opacity)));
    }

    //=== Calc quaternion
    // This converts a 4D vector to a 3x3 Matrix.
    // All quat to rot matrix funcs does normalization. I already did it up. But even if you do again, L2 norm of L2 norm gives the same answer.
    float *rotations = NULL;
    rotations = (float *)malloc(num_gaussians * 9 * sizeof(float)); // I know for sure Im setting all values of this 9-element array to something
    if (rotations == NULL)
    {
        printf("Mem. Alloc. failed\n");
        return 1;
    }
    for (i = 0; i < num_gaussians; i++)
    {
        float r = gaussians[i].quaternion[0];
        float x = gaussians[i].quaternion[1];
        float y = gaussians[i].quaternion[2];
        float z = gaussians[i].quaternion[3];

        rotations[(i * 9) + 0] = 1 - (2 * ((y * y) + (z * z)));
        rotations[(i * 9) + 1] = 2 * ((x * y) - (r * z));
        rotations[(i * 9) + 2] = 2 * ((x * z) + (r * y));
        rotations[(i * 9) + 3] = 2 * ((x * y) + (r * z));
        rotations[(i * 9) + 4] = 1 - (2 * ((x * x) + (z * z)));
        rotations[(i * 9) + 5] = 2 * ((y * z) - (r * x));
        rotations[(i * 9) + 6] = 2 * ((x * z) - (r * y));
        rotations[(i * 9) + 7] = 2 * ((y * z) + (r * x));
        rotations[(i * 9) + 8] = 1 - (2 * ((x * x) + (y * y)));
    }

    //== Create scale matrix
    // [(0,1,2),
    //  (3,4,5),
    //  (6,7,8)]
    // [(00, 01, 02),
    //  (10, 11, 12),
    //  (20, 21, 22)]
    float *scales = NULL;
    scales = (float *)calloc(num_gaussians * 9, sizeof(float));
    if (scales == NULL)
    {
        printf("Mem. Alloc. failed\n");
        return 1;
    }
    for (i = 0; i < num_gaussians; i++)
    {
        scales[(i * 9) + 0] = gaussians[i].scale[0];
        scales[(i * 9) + 4] = gaussians[i].scale[1];
        scales[(i * 9) + 8] = gaussians[i].scale[2];
    }

    // and calculate the cov3D
    float *cov3D = NULL;
    float *r_s = NULL;   // to store intermediate results
    float *r_s_t = NULL; // to store intermediate results
    float *sigma = NULL; // to store the final 3x3 matrix of E = R.S.S^T.R^T
    cov3D = (float *)malloc(num_gaussians * 6 * sizeof(float));
    r_s = (float *)malloc(9 * sizeof(float));
    r_s_t = (float *)malloc(9 * sizeof(float));
    sigma = (float *)malloc(9 * sizeof(float));
    if ((cov3D == NULL) || (r_s == NULL) || (r_s_t == NULL) || (sigma == NULL))
    {
        printf("Mem. Alloc. failed\n");
        return 1;
    }

    for (i = 0; i < num_gaussians; i++)
    {
        // R * S
        matmul3x3(rotations + (i * 9), scales + (i * 9), r_s);
        // R * S * S^T * R^T
        transpose3x3(r_s, r_s_t);
        matmul3x3(r_s, r_s_t, sigma);

        // store only the upper right diagnol elements as its symmetric
        // https://github.com/graphdeco-inria/diff-gaussian-rasterization/blob/59f5f77e3ddbac3ed9db93ec2cfe99ed6c5d121d/cuda_rasterizer/forward.cu#L145
        cov3D[(i * 6) + 0] = sigma[0];
        cov3D[(i * 6) + 1] = sigma[1];
        cov3D[(i * 6) + 2] = sigma[2];
        cov3D[(i * 6) + 3] = sigma[4];
        cov3D[(i * 6) + 4] = sigma[5];
        cov3D[(i * 6) + 5] = sigma[8];
    }
    free(r_s);
    free(r_s_t);
    free(sigma);

    // Read CSV file and get camera details
    FILE *csv_file = fopen("cameras.csv", "r");
    if (!csv_file)
    {
        printf("Unable to open CSV file!");
        return 1;
    }

    int img_id = 0;
    char img_name[20];
    int img_width = 0;
    int img_height = 0;
    float cam_pos[3] = {0.0, 0.0, 0.0};
    float cam_rotation[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    float cam_focal_x = 0.0;
    float cam_focal_y = 0.0;
    float cam_fov_x = 0.0;
    float cam_fov_y = 0.0;
    float tan_fov_x = 0.0;
    float tan_fov_y = 0.0;

    fgets(line, sizeof(line), csv_file); // reading to skip the header line
    while (fgets(line, sizeof(line), csv_file) != NULL)
    {
        line[strcspn(line, "\n")] = 0; // put a NULL where the \n is

        sscanf(line, "%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s",
               &img_id, &img_width, &img_height,
               &cam_pos[0], &cam_pos[1], &cam_pos[2],
               &cam_rotation[0], &cam_rotation[1], &cam_rotation[2], &cam_rotation[3], &cam_rotation[4], &cam_rotation[5], &cam_rotation[6], &cam_rotation[7], &cam_rotation[8],
               &cam_focal_x, &cam_focal_y,
               img_name);

        if (img_id == TEST_IMG_ID)
        {
            cam_fov_x = focal2fov(cam_focal_x, img_width);
            cam_fov_y = focal2fov(cam_focal_y, img_height);
            // printf("fov_x:\t%.10f\n", cam_fov_x);
            // printf("fov_y:\t%.10f\n", cam_fov_y);

            tan_fov_x = tan(cam_fov_x * 0.5);
            tan_fov_y = tan(cam_fov_y * 0.5);
            img_width = (int)ceil((float)img_width/TEST_IMG_SCALE);
            img_height = (int)ceil((float)img_height/TEST_IMG_SCALE);

            printf("GT Img:\t%s\n", img_name);
            printf("Img w:\t%d\n", img_width);
            printf("Img h:\t%d\n", img_height);
            printf("fov_x:\t%.10f\n", cam_fov_x);
            printf("fov_y:\t%.10f\n", cam_fov_y);
            printf("tan_x:\t%.10f\n", tan_fov_x);
            printf("tan_y:\t%.10f\n", tan_fov_y);
            // printf("focal_x:\t%.10f\n", cam_focal_x);
            // printf("focal_y:\t%.10f\n", cam_focal_y);
            cam_focal_x = img_width / (2.0 * tan_fov_x);
            cam_focal_y = img_height / (2.0 * tan_fov_y);

            printf("focal_x:\t%.10f\n", cam_focal_x);
            printf("focal_y:\t%.10f\n", cam_focal_y);
            printf("-----------\n");
            break;
        }
    }

    // Now calc rotation and translation of the camera, so that we can calc the world view transform matrix
    // COLMAP says it calculated "Position of Camera" by P = -R^T * T, were R^T is transpose of rotation matrix and T is translation
    // Now we have to find T. So we can do T = P * inv(-R^T), inv() is inverse of the matrix operation
    // https://colmap.github.io/format.html#sparse-reconstruction:~:text=The%20coordinates%20of%20the%20projection/camera%20center%20are%20given%20by%20%2DR%5Et%20*%20T%2C%20where%20R%5Et%20is%20the%20inverse/transpose%20of%20the%203x3%20rotation%20matrix%20composed%20from%20the%20quaternion%20and%20T%20is%20the%20translation%20vector.
    // P is 3x1
    // R is 3x3
    // T is supposed to be 3x1

    // T = P * inv(-R^T)
    float *cam_translation = NULL;
    cam_translation = (float *)malloc(3 * sizeof(float));
    float *cam_rotation_transpose = NULL;
    cam_rotation_transpose = (float *)malloc(9 * sizeof(float));
    float *cam_rotation_inverse = NULL;
    cam_rotation_inverse = (float *)malloc(9 * sizeof(float));
    if ((cam_translation == NULL) || (cam_rotation_transpose == NULL) || (cam_rotation_inverse == NULL)) {
        printf("Mem. Alloc. failed\n");
        return 1;
    }
    transpose3x3(cam_rotation, cam_rotation_transpose);
    float neg_one = -1.0;
    multiply_const_3x3(cam_rotation_transpose, &neg_one, cam_rotation_transpose);
    // inverse3x3(cam_rotation_transpose, cam_rotation_inverse);
    inversemat(cam_rotation_transpose, 3, cam_rotation_inverse);
    matmul1x3(cam_pos, cam_rotation_inverse, cam_translation);
    multiply_const_3x3(cam_rotation_transpose, &neg_one, cam_rotation_transpose); // reverse the sign back.


    // for (i = 0; i < 3; i++) {
    //     printf("%.15f\n", cam_translation[i]);
    // }

    // Get world2view matrix. Or world view transformation matrix
    // https://github.com/graphdeco-inria/gaussian-splatting/blob/8a70a8cd6f0d9c0a14f564844ead2d1147d5a7ac/utils/graphics_utils.py#L38
    float *Rt_temp = NULL;
    Rt_temp = (float *)calloc(16, sizeof(float));
    float *world_view_transform = NULL;
    world_view_transform = (float *)calloc(16, sizeof(float));
    if((Rt_temp == NULL) || (world_view_transform == NULL)) {
        printf("Mem. alloc. failed\n");
        return 1;
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Rt_temp[i*4 + j] = cam_rotation_transpose[i*3 + j];
        }
    }
    for (j = 0; j < 3; j++) {
        Rt_temp[j*4 + 3] = cam_translation[j];
    }
    Rt_temp[15] = 1.0;

    // it is done ! But Im going to add the feature to introduce translate and scale
    float translate[3] = {0.0, 0.0, 0.0};
    float scale = 1.0;
    float *C2W = NULL;
    C2W = (float *)calloc(16, sizeof(float));
    if(C2W == NULL) {
        printf("Mem. alloc. failed\n");
        return 1;
    }
    inversemat(Rt_temp, 4, C2W);
    for (j = 0; j < 3; j++) {
        C2W[j*4 + 3] = (C2W[j*4 + 3] + translate[j]) * scale;
    }
    inversemat(C2W, 4, Rt_temp);
    free(C2W);

    transpose4x4(Rt_temp, world_view_transform);
    free(Rt_temp);

    // calculate projection matrix
    // 00, 01, 02, 03 = 0, 1, 2, 3
    // 10, 11, 12, 13 = 4, 5, 6, 7
    // 20, 21, 22, 23 = 8, 9, 10, 11
    // 30, 31, 32, 33 = 12, 13, 14, 15
    float *projection_matrix_temp = NULL;
    projection_matrix_temp = (float *)calloc(16, sizeof(float));
    float *projection_matrix = NULL;
    projection_matrix = (float *)calloc(16, sizeof(float));
    if ((projection_matrix == NULL) || (projection_matrix_temp == NULL)) {
        printf("Mem. alloc. failed\n");
        return 1;
    }
    float top = tan_fov_y * ZNEAR;
    float bottom = -1 * top;
    float right = tan_fov_x * ZNEAR;
    float left = -1 * right;
    float z_sign = 1.0;
    projection_matrix_temp[0] = 2.0 * (ZNEAR / (right - left));
    projection_matrix_temp[5] = 2.0 * (ZNEAR / (top - bottom));
    projection_matrix_temp[2] = (right + left) / (right - left);
    projection_matrix_temp[6] = (top + bottom) / (top - bottom);
    projection_matrix_temp[10] = z_sign * ((ZFAR) / (ZFAR - ZNEAR));
    projection_matrix_temp[11] = (-1 * (ZFAR * ZNEAR)) / (ZFAR - ZNEAR);
    projection_matrix_temp[14] = z_sign;
    transpose4x4(projection_matrix_temp, projection_matrix);
    // printf("--------\n");
    free(projection_matrix_temp);
    // for (i = 0; i < 4; i++) {
    //     for (j = 0; j < 4; j++) {
    //         printf("%f\t", world_view_transform[i*4 + j]);
    //     }
    //     printf("\n");
    // }

    // for(i = 0; i < 4; i++) {
    //     printf(YEL "viewmatrix[%d]: ", i);
    //     for(j = 0; j < 4; j++) {
    //         printf("%f, ", world_view_transform[i*4 + j]);
    //     }
    //     printf(RESET "\n");
    // }

    // calc full_proj_transform
    float *full_proj_transform = NULL;
    full_proj_transform = (float *)calloc(16, sizeof(float));
    if (full_proj_transform == NULL) {
        printf("Mem. alloc. failed\n");
        return 1;
    }
    matmul(world_view_transform, 4, 4,projection_matrix, 4, 4, full_proj_transform);
    // for(i = 0; i < 4; i++) {
    //     printf(YEL "full_proj[%d]: ", i);
    //     for(j = 0; j < 4; j++) {
    //         printf("%f, ", full_proj_transform[i*4 + j]);
    //     }
    //     printf(RESET "\n");
    // }
    
    // check for in frustrum
    bool gaussian_validity[num_gaussians];
    for (i = 0; i < num_gaussians; i++){
        float p_hom[4];
        float p_orig[3];
        p_orig[0] = gaussians[i].x;
        p_orig[1] = gaussians[i].y;
        p_orig[2] = gaussians[i].z;
        transformPoint4x4(p_orig, full_proj_transform, p_hom);

        float p_w = 1.00 / (p_hom[3] + 0.0000001);
        float p_proj[3];
        p_proj[0] = p_hom[0] * p_w;
        p_proj[1] = p_hom[1] * p_w;
        p_proj[2] = p_hom[2] * p_w;

        float p_view[3];
        transformPoint4x3(p_orig, world_view_transform, p_view);

        if (i == DEBUG_IDX) {
            printf(RED "p_orig: %f\n" RESET, p_orig[0]);
            printf(RED "p_orig: %f\n" RESET, p_orig[1]);
            printf(RED "p_orig: %f\n" RESET, p_orig[2]);
            printf("--------\n");
            printf(RED "p_view: %f\n" RESET, p_view[0]);
            printf(RED "p_view: %f\n" RESET, p_view[1]);
            printf(RED "p_view: %f\n" RESET, p_view[2]);
        }

        if (p_view[2] <= 0.2) {
            // printf("False\n");
            gaussian_validity[i] = false;
            continue;
        }
        gaussian_validity[i] = true;
    }


    // calculate tile grid
    int tile_grid[3];
    int block[3];
    tile_grid[0] = (img_width + BLOCK_X - 1)/BLOCK_X;
    tile_grid[1] = (img_height + BLOCK_Y - 1)/BLOCK_Y;
    tile_grid[2] = 1;
    block[0] = BLOCK_X;
    block[1] = BLOCK_Y;
    block[2] = 1;
    printf("tile_grid[0]: %d\n", tile_grid[0]);
    printf("tile_grid[1]: %d\n", tile_grid[1]);
    printf("tile_grid[2]: %d\n", tile_grid[2]);
    printf("--------\n");
    
    

    // Calculate 2D covariances
    float *cov2D = NULL;
    cov2D = (float *)malloc(num_gaussians * 3 * sizeof(float));
    float *rgb = NULL;
    bool *clamped = NULL;
    rgb = (float *)malloc(num_gaussians * 3 * sizeof(float));
    clamped = (bool *)malloc(num_gaussians * 3 * sizeof(bool));
    for (i = 0; i < num_gaussians*3; i++) {
        clamped[i] = false;
    }
    float *depths = NULL;
    depths = (float *)malloc(num_gaussians * sizeof(float));
    unsigned int *radii = NULL;
    radii = (unsigned int *)malloc(num_gaussians * sizeof(unsigned int));
    float *points_xy_image = NULL;
    points_xy_image = (float *)malloc(num_gaussians * 2 * sizeof(float));
    float *conic_opacity = NULL;
    conic_opacity = (float *)malloc(num_gaussians * 4 * sizeof(float));
    unsigned int *tiles_touched = NULL;
    tiles_touched = (unsigned int *)malloc(num_gaussians * 1 * sizeof(unsigned int));
    unsigned int *g_rect_min = NULL;
    unsigned int *g_rect_max = NULL;
    g_rect_min = (unsigned int *)malloc(num_gaussians * 2 * sizeof(unsigned int));
    g_rect_max = (unsigned int *)malloc(num_gaussians * 2 * sizeof(unsigned int));

    if((cov2D == NULL) || (rgb == NULL) || (clamped == NULL) || (depths == NULL) || (radii == NULL) || (points_xy_image == NULL) || (conic_opacity == NULL) || (tiles_touched == NULL) || (g_rect_min == NULL) || (g_rect_max == NULL) )  {
        printf("mem. Alloc. failed\n");
        return 1;
    }
    

    for (i = 0; i < num_gaussians; i++) {
        // if(i == DEBUG_IDX) { // ignore this ! just debug info ! wont make sense to you ! 
        //     printf("tile_grid[0]: %d\n", tile_grid[0]);
        //     printf("tile_grid[1]: %d\n", tile_grid[1]);
        //     printf("tile_grid[2]: %d\n", tile_grid[2]);
        // }
        if (gaussian_validity[i] == true) {
            float p_view[3];
            
            float p_orig[3];
            p_orig[0] = gaussians[i].x;
            p_orig[1] = gaussians[i].y;
            p_orig[2] = gaussians[i].z;
            transformPoint4x3(p_orig, world_view_transform, p_view);

            float p_hom[4];
            transformPoint4x4(p_orig, full_proj_transform, p_hom);
            if (i == DEBUG_IDX) printf(CYN "p_hom: %f, %f, %f, %f\n" RESET, p_hom[0], p_hom[1], p_hom[2], p_hom[3]);

            float p_w = 1.0 / (p_hom[3] + 0.0000001);
            if (i == DEBUG_IDX) printf(CYN "p_w: %f\n" RESET, p_w);

            float p_proj[3];
            p_proj[0] = p_hom[0] * p_w;
            p_proj[1] = p_hom[1] * p_w;
            p_proj[2] = p_hom[2] * p_w;
            if (i == DEBUG_IDX) printf(CYN "p_proj: %f, %f, %f\n" RESET, p_proj[0], p_proj[1], p_proj[2]);



            float t[3] = {0.0, 0.0, 0.0};
            t[0] =    gaussians[i].x * world_view_transform[0]
                    + gaussians[i].y * world_view_transform[4]
                    + gaussians[i].z * world_view_transform[8]
                    + world_view_transform[12];
            t[1] =    gaussians[i].x * world_view_transform[1]
                    + gaussians[i].y * world_view_transform[5]
                    + gaussians[i].z * world_view_transform[9]
                    + world_view_transform[13];
            t[2] =    gaussians[i].x * world_view_transform[2]
                    + gaussians[i].y * world_view_transform[6]
                    + gaussians[i].z * world_view_transform[10]
                    + world_view_transform[14];

            float limx = 1.3 * tan_fov_x;
            float limy = 1.3 * tan_fov_y;
            float txtz = t[0] / t[2];
            float tytz = t[1] / t[2];
            t[0] = min(limx, max(-1 * limx, txtz)) * t[2];
            t[1] = min(limy, max(-1 * limy, tytz)) * t[2];

            float J[9];
            J[0] = cam_focal_x / t[2];
            J[1] = 0.0;
            J[2] = (-1 * (cam_focal_x * t[0])) / (t[2] * t[2]);
            J[3] = 0.0;
            J[4] = cam_focal_y / t[2];
            J[5] = (-1 * (cam_focal_y * t[1])) / (t[2] * t[2]);
            J[6] = 0.0;
            J[7] = 0.0;
            J[8] = 0.0;
            // if (i == 0) {
            //     printf("J[0]: %f, %f, %f\n", J[0], J[1], J[2]);
            //     printf("J[1]: %f, %f, %f\n", J[3], J[4], J[5]);
            //     printf("J[2]: %f, %f, %f\n", J[6], J[7], J[8]);
            //     printf("--------\n");
            // }

            float W[9];
            W[0] = world_view_transform[0];
            W[1] = world_view_transform[4];
            W[2] = world_view_transform[8];
            W[3] = world_view_transform[1];
            W[4] = world_view_transform[5];
            W[5] = world_view_transform[9];
            W[6] = world_view_transform[2];
            W[7] = world_view_transform[6];
            W[8] = world_view_transform[10];
            // if (i == 0) {
            //     printf("W[0]: %f, %f, %f\n", W[0], W[1], W[2]);
            //     printf("W[1]: %f, %f, %f\n", W[3], W[4], W[5]);
            //     printf("W[2]: %f, %f, %f\n", W[6], W[7], W[8]);
            //     printf("--------\n");
            // }

            float T[9];
            matmul3x3(J, W, T);
            // if (i == 0) {
            //     printf("T[0]: %f, %f, %f\n", T[0], T[1], T[2]);
            //     printf("T[1]: %f, %f, %f\n", T[3], T[4], T[5]);
            //     printf("T[2]: %f, %f, %f\n", T[6], T[7], T[8]);
            //     printf("--------\n");
            // }

            float Vrk[9];
            Vrk[0] = cov3D[i*6 + 0];
            Vrk[1] = cov3D[i*6 + 1];
            Vrk[2] = cov3D[i*6 + 2];
            Vrk[3] = cov3D[i*6 + 1];
            Vrk[4] = cov3D[i*6 + 3];
            Vrk[5] = cov3D[i*6 + 4];
            Vrk[6] = cov3D[i*6 + 2];
            Vrk[7] = cov3D[i*6 + 4];
            Vrk[8] = cov3D[i*6 + 5];
            // if (i == 0) {
            //     printf("Vrk[0]: %f, %f, %f\n", Vrk[0], Vrk[1], Vrk[2]);
            //     printf("Vrk[1]: %f, %f, %f\n", Vrk[3], Vrk[4], Vrk[5]);
            //     printf("Vrk[2]: %f, %f, %f\n", Vrk[6], Vrk[7], Vrk[8]);
            //     printf("--------\n");
            // }

            float T_transpose[9];
            transpose3x3(T, T_transpose);

            // float Vrk_transpose[9];
            // transpose3x3(Vrk, Vrk_transpose);

            float cov_temp[9];
            matmul3x3(T, Vrk, cov_temp);
            
            float cov[9];
            matmul3x3(cov_temp, T_transpose, cov);
            // if (i == 0) {
            //     printf("cov[0]: %f, %f, %f\n", cov[0], cov[1], cov[2]);
            //     printf("cov[1]: %f, %f, %f\n", cov[3], cov[4], cov[5]);
            //     printf("cov[2]: %f, %f, %f\n", cov[6], cov[7], cov[8]);
            //     printf("--------\n");
            // }

            cov[0] = cov[0] + 0.3;
            cov[4] = cov[4] + 0.3;
        
            cov2D[i*3 + 0] = cov[0];
            cov2D[i*3 + 1] = cov[1];
            cov2D[i*3 + 2] = cov[4];


            // Invert covariance
            float det = (cov2D[i*3 + 0] * cov2D[i*3 + 2]) - (cov2D[i*3 + 1]*cov2D[i*3 + 1]);
            if (det == 0.0){
                gaussian_validity[i] = false;
                printf("det zero\n");
                continue;
            }

            float det_inv = 1.0 / det;
            float conic[3];
            conic[0] = cov2D[i*3 + 2] * det_inv;
            conic[1] = -1 * cov2D[i*3 + 1] * det_inv;
            conic[2] = cov2D[i*3 + 0] * det_inv;
            if (i == DEBUG_IDX) printf(CYN "det: %f\n" RESET, det);
            if (i == DEBUG_IDX) printf(CYN "det_inv: %f\n" RESET, det_inv);
            if (i == DEBUG_IDX) printf(CYN "conic: %f, %f, %f\n" RESET, conic[0], conic[1], conic[2]);

            float mid = 0.5 * (cov2D[i*3 + 0] + cov2D[i*3 + 2]);
            float lambda1 = mid + sqrt(max(0.1, mid * mid - det));
            float lambda2 = mid - sqrt(max(0.1, mid * mid - det));
            float my_radius = ceil(3.0 * sqrt(max(lambda1, lambda2)));
            float point_image[2];
            point_image[0] = ndc2Pix(p_proj[0], img_width); // convert ndc point to pixel coordinate ndc2Pix(v, S)- ((v + 1.0) * S - 1.0) * 0.5;
            point_image[1] = ndc2Pix(p_proj[1], img_height);
            unsigned int rect_min[2], rect_max[2];           // THIS RECT REPRESENTS THE GAUSSIAN SIZE IN TILE RANGES, NOT PIXELS
            rect_min[0] = (int)min(tile_grid[0], max((int)0, (int)((point_image[0] - my_radius)/ BLOCK_X)));
            rect_min[1] = (int)min(tile_grid[1], max((int)0, (int)((point_image[1] - my_radius)/ BLOCK_Y)));
            rect_max[0] = (int)min(tile_grid[0], max((int)0.0, (int)((point_image[0] + my_radius + BLOCK_X -1) / BLOCK_X)));
            //              min(grid.x, max((int)0, (int)((p.x + max_radius + BLOCK_X - 1) / BLOCK_X))),
            rect_max[1] = (int)min(tile_grid[1], max((int)0, (int)((point_image[1] + my_radius + BLOCK_Y -1) / BLOCK_Y)));

            if (i == DEBUG_IDX) printf(CYN "mid: %f\n" RESET, mid);
            if (i == DEBUG_IDX) printf(CYN "lamba1: %f\n" RESET, lambda1);
            if (i == DEBUG_IDX) printf(CYN "lamba2: %f\n" RESET, lambda2);
            if (i == DEBUG_IDX) printf(CYN "my_radius: %f\n" RESET, my_radius);
            if (i == DEBUG_IDX) printf(YEL "W: %d\n" RESET, img_width);
		    if (i == DEBUG_IDX) printf(YEL "H: %d\n" RESET, img_height);
		    if (i == DEBUG_IDX) printf(CYN "point_image.x: %f\n" RESET, point_image[0]);
		    if (i == DEBUG_IDX) printf(CYN "point_image.y: %f\n" RESET, point_image[1]);
		    if (i == DEBUG_IDX) printf(CYN "rect_min.x: %d\n" RESET, rect_min[0]);
		    if (i == DEBUG_IDX) printf(CYN "rect_min.y: %d\n" RESET, rect_min[1]);
		    if (i == DEBUG_IDX) printf(CYN "rect_max.x: %d\n" RESET, rect_max[0]);
		    if (i == DEBUG_IDX) printf(CYN "rect_max.y: %d\n" RESET, rect_max[1]);

            if ((rect_max[0] - rect_min[0]) * (rect_max[1] - rect_min[1]) == 0) {
                gaussian_validity[i] = false;
                continue;
            }

            // computeColorFromSH(int idx, int Degree, int Max_coef, (glm::vec3*)orig_points, vec3 *camera_position, float* shs, bool* clamped);
            // degree 3
            // p_orig
		    if (i == DEBUG_IDX) printf(GRN "p_orig: %f, %f, %f\n" RESET, p_orig[0], p_orig[1], p_orig[2]);
		    if (i == DEBUG_IDX) printf(GRN "cam_pos: %f, %f, %f\n" RESET, cam_pos[0], cam_pos[1], cam_pos[2]);
            float dir[3];
            dir[0] = p_orig[0] - cam_pos[0];
            dir[1] = p_orig[1] - cam_pos[1];
            dir[2] = p_orig[2] - cam_pos[2];
		    if (i == DEBUG_IDX) printf(GRN "dir: %f, %f, %f\n" RESET, dir[0], dir[1], dir[2]);
            float dir_lenght = sqrt((dir[0] * dir[0]) + (dir[1] * dir[1]) + (dir[2] * dir[2]));
            dir[0] = dir[0] / dir_lenght;
            dir[1] = dir[1] / dir_lenght;
            dir[2] = dir[2] / dir_lenght;
		    if (i == DEBUG_IDX) printf(GRN "dir: %f, %f, %f\n" RESET, dir[0], dir[1], dir[2]);

            int coeff = (MAX_SH_DEGREE + 1) * (MAX_SH_DEGREE + 1);
		    if (i == DEBUG_IDX) printf(CYN "coeff: %d\n" RESET, coeff);

            float sh_temp[48];
            for (j = 0; j < 3; j++) {
                sh_temp[j] = gaussians[i].sh_coeff[j];
                // printf(MAG"%d\n"RESET, j);
            }
            // printf("--------\n");
            for (j = 0; j < 3; j++) {
                for (k = 1; k <= 15; k++){
                    sh_temp[k*3 + j] = gaussians[i].sh_coeff[((j*15)+k)+2];
                    // printf(MAG "sh_temp[%d][%d](%d) = g[%d] %f\n" RESET, k, j, k*3 + j,((j*15)+k)+2, gaussians[i].sh_coeff[((j*15)+k)+2]);
                }
                // printf("--------\n");
            }
            // for (j = 0; j < 16; j++) {
            //     printf(BLU"sh[%d, %d, %d]\t: %.10f, %.10f, %.10f\n", (j*3)+0,(j*3)+1,(j*3)+2, sh_temp[(j*3)+0], sh_temp[(j*3)+1], sh_temp[(j*3)+2]);    
            // }

            float rgb_result[3]; // to store RGB
            for (j = 0; j < 3; j++) {
                rgb_result[j] = SH_C0 * sh_temp[j];
            }
            if (i == DEBUG_IDX) printf(MAG "result: %f, %f, %f\n" RESET, rgb_result[0], rgb_result[1], rgb_result[2]);

            if(MAX_SH_DEGREE > 0) {
                float x = dir[0];
                float y = dir[1];
                float z = dir[2];
                for (j = 0; j < 3; j++) {
                    rgb_result[j] = rgb_result[j] - SH_C1 * y * sh_temp[3 + j] + SH_C1 * z * sh_temp[6 + j] - SH_C1 * x * sh_temp[9 + j];  
                }
                if (i == DEBUG_IDX) printf(MAG "result: %f, %f, %f\n" RESET, rgb_result[0], rgb_result[1], rgb_result[2]);

                if(MAX_SH_DEGREE > 1) {
                    float xx = x * x, yy = y * y, zz = z * z;
			        float xy = x * y, yz = y * z, xz = x * z;
                    for (j = 0; j < 3; j++) {
                        rgb_result[j] = rgb_result[j] + 
                                        SH_C2[0] * xy * sh_temp[12 + j] +
                                        SH_C2[1] * yz * sh_temp[15 + j] +
                                        SH_C2[2] * (2.0 * zz - xx - yy) * sh_temp[18 + j] +
                                        SH_C2[3] * xz * sh_temp[21 + j] +
                                        SH_C2[4] * (xx - yy) * sh_temp[24 + j];
                    }
                    if (i == DEBUG_IDX) printf(MAG "result: %f, %f, %f\n" RESET, rgb_result[0], rgb_result[1], rgb_result[2]);

                    if(MAX_SH_DEGREE > 2) {
                        for (j = 0; j < 3; j++) {
                            rgb_result[j] = rgb_result[j] +
                                            SH_C3[0] * y * (3.0 * xx - yy) * sh_temp[27 + j] +
                                            SH_C3[1] * xy * z * sh_temp[30 + j] +
                                            SH_C3[2] * y * (4.0 * zz - xx - yy) * sh_temp[33 + j] +
                                            SH_C3[3] * z * (2.0 * zz - 3.0 * xx - 3.0 * yy) * sh_temp[36 + j] +
                                            SH_C3[4] * x * (4.0 * zz - xx - yy) * sh_temp[39 + j] +
                                            SH_C3[5] * z * (xx - yy) * sh_temp[42 + j] +
                                            SH_C3[6] * x * (xx - 3.0 * yy) * sh_temp[45 + j];
                        }
                        if (i == DEBUG_IDX) printf(MAG "result: %f, %f, %f\n" RESET, rgb_result[0], rgb_result[1], rgb_result[2]);            
                    }
                }
            }    
            for (j = 0; j < 3; j++) {
                rgb_result[j] = rgb_result[j] + 0.5;
            }
            if (i == DEBUG_IDX) printf(MAG "result: %f, %f, %f\n" RESET, rgb_result[0], rgb_result[1], rgb_result[2]);            
            
            // do clamping here
            if (rgb_result[0] < 0) clamped[i*3 + 0] = true;
            if (rgb_result[1] < 0) clamped[i*3 + 1] = true;
            if (rgb_result[2] < 0) clamped[i*3 + 2] = true;

            for (j = 0; j < 3; j++) {
                rgb_result[j] = max(rgb_result[j], 0.0);
                rgb[i*3 + j] = rgb_result[j];
            }
            if (i == DEBUG_IDX) printf(MAG "rgb: %f, %f, %f\n" RESET, rgb[i*3 + 0], rgb[i*3 + 1], rgb[i*3 + 2]);            

            depths[i] = p_view[2];
            radii[i] = (int)my_radius;
            points_xy_image[i*2 + 0] = point_image[0];
            points_xy_image[i*2 + 1] = point_image[1];
            conic_opacity[i*4 + 0] = conic[0];
            conic_opacity[i*4 + 1] = conic[1];
            conic_opacity[i*4 + 2] = conic[2];
            conic_opacity[i*4 + 3] = gaussians[i].opacity;
            tiles_touched[i] = (rect_max[1] - rect_min[1]) * (rect_max[0] - rect_min[0]);
            g_rect_min[i*2 + 0] = rect_min[0];
            g_rect_min[i*2 + 1] = rect_min[1];
            g_rect_max[i*2 + 0] = rect_max[0];
            g_rect_max[i*2 + 1] = rect_max[1];

            if (i == DEBUG_IDX) printf(MAG "tiles_touches[%d]: %d\n" RESET, i, tiles_touched[i]);            
        } else {
            // printf("HEREEEEEEE\n");
            continue;
        }

        // break;
    }


    // printf(CYN "DEBUG_IDX - %d\n" RESET, DEBUG_IDX);
    // printf(CYN "Depths[%d]: %f\n" RESET, DEBUG_IDX, depths[DEBUG_IDX]);
    // printf(CYN "Radii[%d]: %d\n" RESET, DEBUG_IDX, radii[DEBUG_IDX]);
    // printf(CYN "tiles_touched[%d]: %d\n" RESET, DEBUG_IDX, tiles_touched[DEBUG_IDX]);
    // printf(CYN "conic_opacity[%d].x: %f\n" RESET, DEBUG_IDX, conic_opacity[DEBUG_IDX*4 + 0]);
    // printf(CYN "conic_opacity[%d].y: %f\n" RESET, DEBUG_IDX, conic_opacity[DEBUG_IDX*4 + 1]);
    // printf(CYN "conic_opacity[%d].z: %f\n" RESET, DEBUG_IDX, conic_opacity[DEBUG_IDX*4 + 2]);
    // printf(CYN "conic_opacity[%d].w: %f\n" RESET, DEBUG_IDX, conic_opacity[DEBUG_IDX*4 + 3]);
    // printf(CYN "Points_xy_image[%d].x: %f\n" RESET, DEBUG_IDX, points_xy_image[DEBUG_IDX*2 + 0]);
    // printf(CYN "Points_xy_image[%d].y: %f\n" RESET, DEBUG_IDX, points_xy_image[DEBUG_IDX*2 + 1]);



    // printf("cov2d[%d] - %f\n", (i*3 + 0), cov2D[i*3 + 0]);
    // printf("cov2d[%d] - %f\n", (i*3 + 1), cov2D[i*3 + 1]);
    // printf("cov2d[%d] - %f\n", (i*3 + 2), cov2D[i*3 + 2]);
    // printf("--------\n");

    // printf("cov3d[%d] - %f\n", (i*6 + 0), cov3D[i*6 + 0]);
    // printf("cov3d[%d] - %f\n", (i*6 + 1), cov3D[i*6 + 1]);
    // printf("cov3d[%d] - %f\n", (i*6 + 2), cov3D[i*6 + 2]);
    // printf("cov3d[%d] - %f\n", (i*6 + 3), cov3D[i*6 + 3]);
    // printf("cov3d[%d] - %f\n", (i*6 + 4), cov3D[i*6 + 4]);
    // printf("cov3d[%d] - %f\n", (i*6 + 5), cov3D[i*6 + 5]);
    // printf("--------\n");

    // Now the remaining forward rasterization
    // inclusive sum
    // in = [2, 3, 0, 1, 4]
    // out = [2, 5, 5, 6, 10]
    int *inclusive_sum_results = NULL;
    inclusive_sum_results = (int *)malloc(num_gaussians * sizeof(int));
    if((inclusive_sum_results == NULL)) {
        printf("mem. Alloc. failed\n");
        return 1;
    }
    InclusiveSum(tiles_touched, inclusive_sum_results, num_gaussians);


    // int count = 0;
    // count the number of valid gaussians
    // for (i = 0; i < num_gaussians; i++){
    //     if(gaussian_validity[i] == true) {
    //         count += 1;
    //     }
    // }
    // printf("Valid Count: %d/%d\n", count, num_gaussians);

    //count valid gaussians with >0 radius
    // count = 0;
    // for (i = 0; i < num_gaussians; i++){
    //     if(gaussian_validity[i] == true) {
    //         if(radii[i] > 0) {
    //             count += 1;
    //         }
    //     }
    // }
    // printf("Valid Count: %d/%d\n", count, num_gaussians);

    // int valid_indices[count];
    // count = 0;
    // store all the valid gaussian indexes of master array
    // for (i = 0; i < num_gaussians; i++){
    //     if(gaussian_validity[i] == true) {
    //         if(radii[i] > 0) {
    //             valid_indices[count] = i;
    //             count += 1;
    //         }
    //     }
    // }
    // printf("Valid Count: %d/%d\n", count, num_gaussians);


    // count total tiles touched
    // int num_rendered = 0;
    // for (i = 0; i < num_gaussians; i++){
    //     if(gaussian_validity[i] == true) {
    //         if(radii[i] > 0) {
    //             // valid_indices[count] = i;
    //             // count += 1;
    //             num_rendered += tiles_touched[i];
    //         }
    //     }
    // }
    // printf("num_rendered: %d\n", num_rendered);

    // int idx = 285983;
    // printf(RED "tmp_tile_t[%d]: %d\n"RESET, idx, tiles_touched[idx]);
    // printf(RED "tmp_pt_off[%d]: %d\n"RESET, idx, inclusive_sum_results[idx]);
    // printf("--------\n");
    // printf(RED "tmp_tile_t[%d]: %d\n"RESET, idx, tiles_touched[idx-1]);
    // printf(RED "tmp_pt_off[%d]: %d\n"RESET, idx, inclusive_sum_results[idx-1]);
    int num_rendered = inclusive_sum_results[num_gaussians-1];
    printf(RED "num_rendered: %d\n"RESET, num_rendered);

    printf(GRN "Preprocessing done!\n" RESET);
    tik2 = clock();
    printf(CYN "Preprocessing Time: %f secs\n" RESET, (double)(tik2 - tik1) / CLOCKS_PER_SEC);

    // create duplicate with keys for sorting ! 
    // Create key value for all gaussian/tile overlap
    // g_kv *gaussian_keys = NULL;
    // gaussian_keys = (g_kv*)malloc(num_rendered * sizeof(g_kv));
    uint64_t *keys = NULL;
    keys = (uint64_t *)malloc(num_rendered * sizeof(uint64_t));
    unsigned int *gaussian_values = NULL;
    gaussian_values = (unsigned int *)malloc(num_rendered * sizeof(unsigned int));
    
    if((keys == NULL) || (gaussian_values == NULL)) {
        printf("mem. Alloc. failed\n");
        return 1;
    }

    int offset = 0;
    for(i = 0; i < num_gaussians; i++) {
        if (radii[i] > 0) {
		    // if (i == DEBUG_IDX) printf(BLU "===== min(%d, %d), max(%d,%d)\n" RESET, g_rect_min[i*2 + 0], g_rect_min[i*2 + 1], g_rect_max[i*2 + 0], g_rect_max[i*2 + 1]);
            for(int y = g_rect_min[i*2 + 1]; y < g_rect_max[i*2 + 1]; y++) {
                for(int x = g_rect_min[i*2 + 0]; x < g_rect_max[i*2 + 0]; x++) {
                    uint64_t key = y * tile_grid[0] + x;
				    if (i == DEBUG_IDX) printf(BLU "===== key1: %lu\n" RESET, key);
                    key = key << 32;
				    if (i == DEBUG_IDX) printf(BLU "===== key2: %lu\n" RESET, key);
				    if (i == DEBUG_IDX) printf(BLU "===== dep: %f\n" RESET, depths[i]);
				    if (i == DEBUG_IDX) printf(BLU "===== dep: %p\n" RESET, &depths[i]);
				    if (i == DEBUG_IDX) printf(BLU "===== dep: %p\n" RESET, (uint32_t*)&depths[i]);
				    if (i == DEBUG_IDX) printf(BLU "===== dep: %d\n" RESET, *((uint32_t*)&depths[i]));
                    key = key | *((uint32_t*)&depths[i]);
				    if (i == DEBUG_IDX) printf(BLU "===== key3: %lu\n" RESET, key);
                    if (i == DEBUG_IDX) printf("--------\n");
                    keys[offset] = key;
                    gaussian_values[offset] = i;
                    offset++;
                }
            }
        }
    }

    int highestBit = getHigherMsb(tile_grid[0] * tile_grid[1]);
    radixSortPairs(keys, gaussian_values, num_rendered, (32+highestBit));

    printf(GRN "Sorting done!\n" RESET);
    tik3 = clock();
    printf(CYN "Sorting Time: %f secs\n" RESET, (double)(tik3 - tik2) / CLOCKS_PER_SEC);

    // for (i = 9000; i < 10000; i++) {
    //     uint32_t current_tile = keys[i] >> 32;
    //     printf("%d: %d\n", i, current_tile);
    // }

    // for (i = 0; i < 10; i++) {
    //     printf("%d: %d\n", i, gaussian_values[i]);
    // }



    // Identify tile ranges ! 
    int *tile_ranges = NULL;
    tile_ranges = (int *)malloc(tile_grid[0] * tile_grid[1] * 2 * sizeof(int));
    if((tile_ranges == NULL)) {
        printf("mem. Alloc. failed\n");
        return 1;
    }

    for (i = 0; i < num_rendered; i++) {
        uint32_t current_tile = keys[i] >> 32;

        if(i == 0) {
            tile_ranges[current_tile*2 + 0] = 0;
        } else {
            uint32_t previous_tile = keys[i - 1] >> 32;
            if (current_tile != previous_tile) {
                tile_ranges[previous_tile*2 + 1] = i;
                tile_ranges[current_tile*2 + 0] = i;
            }
        }

        if (i == (num_rendered - 1)) {
            tile_ranges[current_tile*2 + 1] = num_rendered;
        }
    }

    // int idx = 1;
    // printf("tile_ranges[%d].x: %d\n", idx, tile_ranges[idx*2 + 0]);
    // printf("tile_ranges[%d].y: %d\n", idx, tile_ranges[idx*2 + 1]);

    // finally fucking render!!!
    // things needed
    // tile_grid, block size, tile_ranges, keys, img_width, img_height,
    // points_xy_image(means2d in off impl.)
    // rgb(feature_ptr)
    // conic_opacity
    // return:
    //      accum_alpha(float)?
    //      n_contrib(uint32_t)?
    //      out_color
    //
    // one tile = one block in CUDA
    float* blocks_rgb = NULL;
    blocks_rgb = (float*)malloc(img_width * img_height * 3 * sizeof(float));
    if((blocks_rgb == NULL)) {
        printf("mem. Alloc. failed\n");
        return 1;
    }

    // setting all RGB vals to zero
    for (i = 0; i < (img_width * img_height * 3); i++) {
        blocks_rgb[i] = 0;
    }
    
    

    // int idx_i = 1;
    // int idx_j = 0;
    int idx_i = 9999;
    int idx_j = 9999;
    int pixel_idx;
    for (i = 0; i < tile_grid[0]; i++) {
        for (j = 0; j < tile_grid[1]; j++) {
            tik4 = clock();
            //get the range for this block
            uint32_t horizontal_blocks = (img_width + BLOCK_X - 1) / BLOCK_X;
            int pix_min[2], pix_max[2], pix[2];
            pix_min[0] = i * BLOCK_X;
            pix_min[1] = j * BLOCK_Y;
            pix_max[0] = (int)min(pix_min[0] + BLOCK_X, img_width);
            pix_max[1] = (int)min(pix_min[1] + BLOCK_Y, img_height);

            int from = tile_ranges[(j * horizontal_blocks + i)*2 + 0];
            int to = tile_ranges[(j * horizontal_blocks + i)*2 + 1];

            // int from = tile_ranges[(*tile_grid[1] + j)*2 + 0];
            // int to = tile_ranges[(i*tile_grid[1] + j)*2 + 1];
            
            if ((idx_i == i) && (idx_j == j)) printf("B-(%d,%d) horizontal_blocks: %d\n", i, j, horizontal_blocks);
            if ((idx_i == i) && (idx_j == j)) printf("B-(%d,%d) pix_min: %d, %d\n", i, j, pix_min[0], pix_min[1]);
            if ((idx_i == i) && (idx_j == j)) printf("B-(%d,%d) pix_max: %d, %d\n", i, j, pix_max[0], pix_max[1]);
            if ((idx_i == i) && (idx_j == j)) printf("B-(%d,%d) range: %d - %d\n", i, j, from, to);

            // The CUDA implementation does this differently. But its all the same algorithm in the official paper. Pg.14
            // for (k = pix_min[0]; k < pix_max[0]; k++) {
                // for (l = pix_min[1]; l < pix_max[1]; l++) {
                for (l = 0; l < BLOCK_Y; l++) {
            for (k = 0; k < BLOCK_X; k++) {
                    unsigned int pix[2];
                    pix[0] = pix_min[0] + k;
                    pix[1] = pix_min[1] + l;
                    pixel_idx = img_width * pix[1] + pix[0];

                    if ((idx_i == i) && (idx_j == j) && (k == 0) && (l == 0)) printf("B-(%d,%d),k-(%d,%d) pix: %d, %d\n", i, j, k, l, pix[0], pix[1]);
                    if ((idx_i == i) && (idx_j == j) && (k == 0) && (l == 0)) printf("B-(%d,%d),k-(%d,%d) pix_id: %d\n", i, j, k, l, pixel_idx);
                    
                    if ((idx_i == i) && (idx_j == j) && (k == 15) && (l == 15)) printf("B-(%d,%d),k-(%d,%d) pix: %d, %d\n", i, j, k, l, pix[0], pix[1]);
                    if ((idx_i == i) && (idx_j == j) && (k == 15) && (l == 15)) printf("B-(%d,%d),k-(%d,%d) pix_id: %d\n", i, j, k, l, pixel_idx);


                    bool inside = (pix[0] < img_width) && (pix[1] < img_height);
                    bool done = !inside;


                    // pixel_idx = (k*pix_max[1] + l);
                    // printf("pixel_idx: %d\n", pixel_idx);
                    // printf("(k-%d, l-%d, pix_min-(%d, %d), pix_max-(%d, %d))pix_idx = %d/%d\n", k, l, pix_min[0], pix_min[1], pix_max[0], pix_max[1], pixel_idx, (img_height*img_width));
                    // return 0;

                    // iterate over every tile/gaussian overlap
                    float T = 1.0;
                    float pixel_color[3] = {0, 0, 0};
                    for (int kk = from; kk < to && !done; kk++) {
                        int gaussian_id = gaussian_values[kk];

                        float xy[2];
                        xy[0] = points_xy_image[gaussian_id*2 + 0];
                        xy[1] = points_xy_image[gaussian_id*2 + 1];

                        float d[2];
                        d[0] = xy[0] - (float)pix[0];
                        d[1] = xy[1] - (float)pix[1];

                        float con_o[4];
                        con_o[0] = conic_opacity[gaussian_id*4 + 0];
                        con_o[1] = conic_opacity[gaussian_id*4 + 1];
                        con_o[2] = conic_opacity[gaussian_id*4 + 2];
                        con_o[3] = conic_opacity[gaussian_id*4 + 3];

                        if ((idx_i == i) && (idx_j == j)) {
                            if ((k == 0) && (l == 0)) {
                                if(kk == from) {
                                    printf("B-(%d,%d),k-(%d,%d)[0] gaussian_id: %d\n", i, j, k, l, gaussian_id);
                                    printf("B-(%d,%d),k-(%d,%d)[0] xy: %f, %f\n", i, j, k, l, xy[0], xy[1]);
                                    printf("B-(%d,%d),k-(%d,%d)[0] con_o: %f, %f, %f, %f\n", i, j, k, l, con_o[0], con_o[1], con_o[2], con_o[3]);
                                    printf("B-(%d,%d),k-(%d,%d)[0] d: %f, %f\n", i, j, k, l, d[0], d[1]);
                                    // printf("--------\n");
                                }
                                if(kk == from+12) {
                                    printf("B-(%d,%d),k-(%d,%d)[12] gaussian_id: %d\n", i, j, k, l, gaussian_id);
                                    printf("B-(%d,%d),k-(%d,%d)[12] xy: %f, %f\n", i, j, k, l, xy[0], xy[1]);
                                    printf("B-(%d,%d),k-(%d,%d)[12] con_o: %f, %f, %f, %f\n", i, j, k, l, con_o[0], con_o[1], con_o[2], con_o[3]);
                                    printf("B-(%d,%d),k-(%d,%d)[12] d: %f, %f\n", i, j, k, l, d[0], d[1]);
                                    // printf("--------\n");
                                }
                                if(kk == from+100) {
                                    printf("B-(%d,%d),k-(%d,%d)[100] gaussian_id: %d\n", i, j, k, l, gaussian_id);
                                    printf("B-(%d,%d),k-(%d,%d)[100] xy: %f, %f\n", i, j, k, l, xy[0], xy[1]);
                                    printf("B-(%d,%d),k-(%d,%d)[100] con_o: %f, %f, %f, %f\n", i, j, k, l, con_o[0], con_o[1], con_o[2], con_o[3]);
                                    printf("B-(%d,%d),k-(%d,%d)[100] d: %f, %f\n", i, j, k, l, d[0], d[1]);
                                    // printf("--------\n");
                                    // return 1;
                                }
                            }
                         }




                        float power;
                        power = -0.5 * (con_o[0] * d[0] * d[0] + con_o[2] * d[1] * d[1]) - con_o[1] * d[0] * d[1]; 
                        if ((idx_i == i) && (idx_j == j)) {
                            if ((k == 0) && (l == 0)) {
                                if(kk == from) {
                                    printf("B-(%d,%d),k-(%d,%d)[0] power: %f\n", i, j, k, l, power);
                                    printf("--------\n");
                                }
                                if(kk == from+12) {
                                    printf("B-(%d,%d),k-(%d,%d)[12] power: %f\n", i, j, k, l, power);
                                    printf("--------\n");
                                }
                                if(kk == from+100) {
                                    printf("B-(%d,%d),k-(%d,%d)[100] power: %f\n", i, j, k, l, power);
                                    printf("--------\n");
                                    // return 1;
                                }
                            }
                         }

                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) g_id[%d] = %d\n", i, j, k, l, kk, gaussian_id);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) x[%d] = %f\n", i, j, k, l, kk, xy[0]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) y[%d] = %f\n", i, j, k, l, kk, xy[1]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) d[%d].x = %f\n", i, j, k, l, kk, d[0]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) d[%d].y = %f\n", i, j, k, l, kk, d[1]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) con_o[%d].x = %f\n", i, j, k, l, kk, con_o[0]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) con_o[%d].y = %f\n", i, j, k, l, kk, con_o[1]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) con_o[%d].z = %f\n", i, j, k, l, kk, con_o[2]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) con_o[%d].w = %f\n", i, j, k, l, kk, con_o[3]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) power[%d] = %f\n", i, j, k, l, kk, power);
                    
                        // printf("(%d, %d | %d, %d) xy[0] = %f\n", i, j, k, l, xy[0]);
                        // printf("(%d, %d | %d, %d) xy[1] = %f\n", i, j, k, l, xy[1]);

                        if (power > 0.0) continue;

                        float alpha;
                        alpha = min(0.99, con_o[3] * exp(power));
                        // printf("power: %f,\texp: %f\n", power, exp(power));
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) alpha[%d] = %f\n", i, j, k, l, kk, alpha);
                        // if (kk < 256) printf("--------\n");


                        if (alpha < 1.0/255.0) {
                            // printf("ALPHA LOW++++\n");
                            continue;
                        }

			            float test_T = T * (1 - alpha);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) test_T[%d] = %f\n", i, j, k, l, kk, test_T);
                        // if (kk < 256) printf("--------\n");

                        if (test_T < 0.0001) {
                            // printf("test_T LOW++++\n");
                            done = true;
                            continue;
                        }


                        blocks_rgb[pixel_idx*3 + 0] += rgb[gaussian_id*3 + 0] * alpha * T;
                        blocks_rgb[pixel_idx*3 + 1] += rgb[gaussian_id*3 + 1] * alpha * T;
                        blocks_rgb[pixel_idx*3 + 2] += rgb[gaussian_id*3 + 2] * alpha * T;
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) C_0[%d] = %f\n", i, j, k, l, kk, blocks_rgb[pixel_idx*3 + 0]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) C_1[%d] = %f\n", i, j, k, l, kk, blocks_rgb[pixel_idx*3 + 1]);
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("(%d, %d | %d, %d) C_2[%d] = %f\n", i, j, k, l, kk, blocks_rgb[pixel_idx*3 + 2]);




                        T = test_T;
                        // if ((idx_i == i) && (idx_j == j)) if (kk < (from + 256)) printf("--------\n");
                    }

                    // break;
                }
                // break;
            }

            tik5 = clock();
            double time_taken_for_block = (tik5 - tik4) / CLOCKS_PER_SEC;
            printProgress((i*tile_grid[1] + j + 1), (tile_grid[0] * tile_grid[1]), time_taken_for_block);
            // break;
        }
        // break;
    }
    printf("\n");
    printf(GRN "Rasterization done!\n" RESET);
    tik4 = clock();
    printf(CYN "Rasterization Time: %f secs\n" RESET, (double)(tik4 - tik3) / CLOCKS_PER_SEC);

    // printf("last pix: %d\n", pixel_idx);



    // scale and clip
    unsigned char *rgb_data = NULL;
    rgb_data = (unsigned char *)malloc(3 * img_height * img_width * sizeof(unsigned char));
    if((rgb_data == NULL)) {
        printf("rgb_data: mem. Alloc. failed\n");
        return 1;
    }

    for (i = 0; i < (img_height * img_width * 3); i++) {
        int scaled_val = (int)ceil(blocks_rgb[i] * 255);

        // if (scaled_val > 255) blocks_rgb[i] = 255;
        // else if (scaled_val < 0) blocks_rgb[i] = 0;
        // else blocks_rgb[i] = scaled_val;

        if (scaled_val > 255) rgb_data[i] = (unsigned char)255;
        else if (scaled_val < 0) rgb_data[i] = (unsigned char)0;
        else rgb_data[i] = (unsigned char)scaled_val;
    }

    // for (i = 0; i < 16; i++) {
    //     printf("[%d]color: %f, %f, %f\n", i, blocks_rgb[i*3 + 0], blocks_rgb[i*3 + 1], blocks_rgb[i*3 + 2]);
    //     // if (i % 16 == 0) printf("--------\n");
    // }
    
    // for (i = 1557; i < 1573; i++) {
    //     printf("[%d]color: %f, %f, %f\n", i, blocks_rgb[i*3 + 0], blocks_rgb[i*3 + 1], blocks_rgb[i*3 + 2]);
    //     // if (i % 16 == 0) printf("--------\n");
    // }


    // write to a PPM file
    FILE *ppm_image = fopen("render_gs_c.ppm", "wb");
    if (!ppm_image) {
        printf("Unable to open file for writing..\n");
        // return 1;
    } else {
        fprintf(ppm_image, "P6\n%d %d\n255\n", img_width, img_height);
        fwrite(rgb_data, 3 * img_height * img_width, 1, ppm_image);
        fclose(ppm_image);
    }

    printf(GRN "Image written!\n" RESET);

//========================== Cleaning =========================================
    if (feof(file) || feof(csv_file))
    {
        printf("End of file reached.\n");
    }
    else if (ferror(file) || ferror(csv_file))
    {
        printf("An error occurred.\n");
    }

    fclose(file);
    fclose(csv_file);
    free(gaussians);
    free(rotations);
    free(scales);
    free(cov3D);
    free(cam_translation);
    free(world_view_transform);
    free(projection_matrix);
    free(cov2D);
    free(full_proj_transform);
    free(inclusive_sum_results);
    free(g_rect_min);
    free(g_rect_max);
    free(gaussian_values);
    free(tile_ranges);
    free(blocks_rgb);
    free(rgb_data);

    printf(GRN "ALL DONE !\n" RESET);
    tok = clock();
    printf(CYN "Total Time: %f secs\n" RESET, (double)(tok - tik1) / CLOCKS_PER_SEC);
    return 0;
}
