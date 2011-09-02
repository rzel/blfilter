#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <sys/time.h>
#include <ostream>
#include "CImg.h"
#include <omp.h>
#include <cuda.h>
using namespace cimg_library;
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define INDEX(i, j, width) i*width + j

/*Function Declarations */
void blfilter(float *paddedImg, float *gaussian_sp, int width, int height, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img);
/* Explicit tiling uses arrays to store tiles */
__global__ void blfilter_cuda(float *paddedImg, float *gaussian_sp, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img);

double my_difftime();
void *MyCalloc(int MemSize);
void  *MyMalloc(int MemSize);

/* Main Function*/
int main(int argc, char* argv[])
{
    int i = 0, j = 0;
    float *gaussian_sp;
    float *Ip_Img, *Op_Img, *paddedImg;

#if defined USE_TILES || EXP_TILES
    if (argc != 7)
    {
        fprintf(stderr, "Usage: %s <Image File> <Filter Radius> <Sigma Spatial> <Sigma Photometric> <Display Output(0/1)> <TileSize>\n", argv[0]);
        exit(0);
    }
    int TileSize = atoi(argv[6]);
#else
    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s <Image File> <Filter Radius> <Sigma Spatial> <Sigma Photometric> <Display Output(0/1)>\n", argv[0]);
        exit(0);
    }
#endif

    char* fileName = argv[1];//1. input image
    int filter_hw = atoi(argv[2]); // filter half width
    float sigma_sp = (float) atof(argv[3]); // sigma_sp (standard deviation for the spatial distance)
    float sigma_ph = (float) atof(argv[4]); // sigma_ph standard deviation for the photometric distance
    int ShowOutput = atoi(argv[5]);
    double start, end;/* To measure Time */

    /* Load image using CImg library */
   CImg<float> image(fileName);

    int cols = image.width();
    int rows = image.height();
    int ds = image.spectrum(); //Gray or RGB (or could be more channels);


/*    int rows = 5000;
    int cols = 5000;
*/
#ifdef USE_TILES
    if ((cols % TileSize != 0) || (rows % TileSize !=0))
    {
        printf("Image width = %d or Height = %d is not divisible by TileSize\nexiting..",cols, rows);
        exit(1);
    }
#endif

    /* Check for erros */
    Ip_Img = (float *) MyMalloc(rows * cols * sizeof(float));
    Op_Img = (float *) MyMalloc(rows * cols * sizeof(float));
    CImg<float> grayImg(rows, cols, 1, 1, 0); //create a single image dz = 1, rgb channel  = 1 (only gray); fill it with 0s.
    CImg<float> filteredImg(rows, cols, 1, 1, 0); //create a single image dz = 1, rgb channel  = 1 (only gray); fill it with 0s.

    if(ds == 3)
    {   //convert the input RGB to gray image. (work with gray image for now)
        for(i = 0; i < rows; i++)
            for(j  = 0; j < cols; j++)
            {
                float gray_pixel = 0.299 * image(i,j,0,0) + 0.587 * image(i,j,0,1) + 0.114 * image(i,j,0,2);
                grayImg(i,j,0,0) = gray_pixel;
                Ip_Img[INDEX(i,j,cols)] = gray_pixel;
            }
    }
    else
    {
        for(i = 0; i < rows; i++)
            for(j  = 0; j < cols; j++) {
                Ip_Img[INDEX(i,j,cols)] = image(i,j,0,0);
            }
        grayImg = image;
    }

    //gaussian kernel for the spatial distances
    gaussian_sp = (float *) MyMalloc( (2 * filter_hw + 1) *  (2 * filter_hw + 1) * sizeof(float));

    for( i = -filter_hw; i<= filter_hw ; i++)
    {
        for(j = -filter_hw; j<= filter_hw; j++)
        {
            gaussian_sp[INDEX((i+filter_hw), (j+filter_hw), (2 * filter_hw + 1))] = exp( -(j*j + i*i) / (2*sigma_sp*sigma_sp) );
        }
    }

    //zero pad the input image to accomodate filtering at the boundary
    paddedImg = (float *) MyMalloc( (rows + 2 * filter_hw) * (cols + 2 * filter_hw) * sizeof(float));
    memset(paddedImg, 0,((rows + 2 * filter_hw) * (cols + 2 * filter_hw) * sizeof(float)));

    for( i = 0; i < rows; i++)
        memcpy(&(paddedImg[INDEX((i+filter_hw),filter_hw,(cols + filter_hw + filter_hw))]), &(Ip_Img[INDEX(i,0,cols)]), cols * sizeof(float));

    //free(Ip_Img);

    printf("start");
    float *d_paddedImg, *d_gaussian_sp, *d_Op_Img; 
    size_t mem_size = (2 * filter_hw + 1) *  (2 * filter_hw + 1) * sizeof(float);
    (cudaMalloc( (void **) &d_gaussian_sp, mem_size)) ; 
    mem_size = (rows + 2 * filter_hw) * (cols + 2 * filter_hw) * sizeof(float);
    (cudaMalloc( (void **) &d_paddedImg,mem_size)); 
    mem_size = (rows * cols * sizeof(float));
    (cudaMalloc((void **) &d_Op_Img, mem_size)); 
    
    (cudaMemcpy( d_gaussian_sp, gaussian_sp, (2 * filter_hw + 1) *  (2 * filter_hw + 1) * sizeof(float) , cudaMemcpyHostToDevice));
    (cudaMemcpy( d_paddedImg, paddedImg, (rows + 2 * filter_hw) * (cols + 2 * filter_hw) * sizeof(float), cudaMemcpyHostToDevice));
    dim3 blockSize = dim3(32,32);
    dim3 numBlocks = dim3(rows/blockSize.x, cols/blockSize.y);
    start = my_difftime();

    blfilter_cuda<<<numBlocks, blockSize>>>(d_paddedImg, d_gaussian_sp, rows, cols, filter_hw, sigma_sp, sigma_ph, d_Op_Img);
    (cudaMemcpy( paddedImg, d_paddedImg, (rows + 2 * filter_hw) * (cols + 2 * filter_hw) * sizeof(float), cudaMemcpyDeviceToHost));

    end = my_difftime() - start;
    printf("ENd");
    printf("Time taken for filtering = %6.2f \n", end);

    for(i = 0; i < rows; i++)
        for(j  = 0; j < cols; j++)
        {
            filteredImg(i,j,0,0) = paddedImg[INDEX(i,j,(cols + 2 * filter_hw))];
        }

    if (ShowOutput == 1)
    {
        CImgDisplay gray_disp(grayImg, "Gray Image");
        CImgDisplay flt_disp(filteredImg, "Filtered Image");

        //create event loop to display the image
        while (!gray_disp.is_closed() ) {
            gray_disp.wait();
        }
    }

    return 0;
}




void blfilter(float *paddedImg, float *gaussian_sp, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img)
{
    float filter_current_pixel = 0; //local window
    float pd = 0;          //photometric distance
    float gaussian_ph = 0; //photometric filter coefficient
    float gaussian_bl = 0; //bilateral filter coefficient
    float filtered_pixel = 0; //filtered value of pixel
    float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient
    int i = 0,j = 0, k = 0, l = 0;

    //image height and width
    int Rows = rows;
    int Cols = cols;
    float sigmasq  = 1/(-2*sigma_ph*sigma_ph);

#ifdef USE_OMP
#pragma omp parallel shared(paddedImg, Op_Img, gaussian_sp, Rows, Cols, filter_hw, sigma_ph)
    {
#pragma omp for private( filter_current_pixel,i,j,k,l,normal_factor,gaussian_ph, gaussian_bl, filtered_pixel,pd) schedule(dynamic) 
#endif
        for( i = filter_hw; i < Rows + filter_hw ; i++)
        {
            for( j = filter_hw; j < Cols + filter_hw; j++)
            {
                filtered_pixel = 0;
                normal_factor  = 0;

                for( k = -filter_hw ; k <= filter_hw ; k++)
                {
                    for( l = -filter_hw ; l <= filter_hw ; l++)
                    {
                        filter_current_pixel = paddedImg[INDEX((i + k), (j + l), (Cols + filter_hw + filter_hw))];
                        //photometric distance = difference in pixel values (squared)
                        pd = paddedImg[INDEX( i, j, (Cols + filter_hw +  filter_hw))] - filter_current_pixel;
                        pd = pd * pd;

                        //gaussian_ph = exp( pd * sigmasq );
                        gaussian_bl = gaussian_ph * gaussian_sp[INDEX((k+filter_hw), (l+filter_hw), (2*filter_hw + 1))] ;

                        filtered_pixel += gaussian_bl * filter_current_pixel;
                        normal_factor += gaussian_bl;
                    }
                }
                filtered_pixel = filtered_pixel / normal_factor;
                Op_Img[INDEX((i-filter_hw), (j-filter_hw), Cols)] = filtered_pixel;
            }
        }
#ifdef USE_OMP
    }
#endif
}






/*----------------------------------------------------------------------------
 *  * Function ........... MyMalloc
 *   *
 *    * Purpose ............ Allocates memory using malloc and does error checking
 *     *
 *      * Arguments Passed ... MemSize - size of memory to be allocated
 *       *
 *        * Return Value ....... Pointer to the allocated memory
 *         *
 *          * Globals Modified ... None
 *           *---------------------------------------------------------------------------*/
void *MyMalloc(int MemSize)
{
    void  *ptr;

    ptr = (void *) malloc(MemSize);
    if (ptr == (void *) NULL)
    {
        printf("Error (MyMalloc) : Unable to allocate memory\n");
        exit(0);
    }
    return ptr;
}


/*----------------------------------------------------------------------------
 *  * Function ........... MyCalloc
 *   *
 *    * Purpose ............ Allocates memory using calloc (for 1 element only)
 *     *                                                       and does error checking
 *      *
 *       * Arguments Passed ... MemSize - size of memory to be allocated
 *        *
 *         * Return Value ....... Pointer to the allocated memory
 *          *
 *           * Globals Modified ... None
 *            *---------------------------------------------------------------------------*/
void *MyCalloc(int MemSize)
{
    void *ptr;

    if ( (ptr = (void *) calloc(1, MemSize)) == (void *) NULL)
    {
        printf("Error (MyCalloc) : Unable to allocate memory\n");
        exit(0);
    }
    return ptr;
}

double my_difftime ()
{
    struct timeval tp;
    struct timezone tzp;
    int i;

    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
};

__global__ void blfilter_cuda(float *paddedImg, float *gaussian_sp, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img)
{
    float filter_current_pixel = 0; //local window
    float pd = 0;          //photometric distance
    float gaussian_ph = 0; //photometric filter coefficient
    float gaussian_bl = 0; //bilateral filter coefficient
    float filtered_pixel = 0; //filtered value of pixel
    float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient
    int i = 0,j = 0, k = 0, l = 0;
    uint tid_i = (blockIdx.x * blockDim.x) + threadIdx.x;
    uint tid_j = (blockIdx.y * blockDim.y) + threadIdx.y;

    //image height and width
    int Rows = rows;
    int Cols = cols;
    float sigmasq  = 1/(-2*sigma_ph*sigma_ph);

    filtered_pixel = 0;
    normal_factor  = 0;
    for( k = -filter_hw ; k <= filter_hw ; k++)
    {
        for( l = -filter_hw ; l <= filter_hw ; l++)
        {
            filter_current_pixel = paddedImg[INDEX((tid_i+k), (tid_j+l), (Cols + filter_hw + filter_hw))];
            //photometric distance = difference in pixel values (squared)
            pd = paddedImg[INDEX( tid_i, tid_j, (Cols + filter_hw +  filter_hw))] - filter_current_pixel;
            pd = pd * pd;

            gaussian_ph = exp( pd * sigmasq );
            gaussian_bl = gaussian_ph * gaussian_sp[INDEX((k+filter_hw), (l+filter_hw), (2*filter_hw + 1))] ;

            filtered_pixel += gaussian_bl * filter_current_pixel;
            normal_factor += gaussian_bl;
        }
    }
    filtered_pixel = filtered_pixel / normal_factor;
    Op_Img[INDEX((tid_i-filter_hw), (tid_j-filter_hw), Cols)] = filtered_pixel;
}
