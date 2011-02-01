#include "def.h"
#include <sys/time.h>
#include "CImg.h"
#include <omp.h>
#include <ia32intrin.h>
using namespace cimg_library;

/*Function Declarations */
void blfilter(float *paddedImg, float *gaussian_sp, int width, int height, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img);
/* Explicit tiling uses arrays to store tiles */
#ifdef EXP_TILES
void blfilter_Read_Tiles(float *paddedImg, float *gaussian_sp, int TileSize, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img);
#endif
#ifdef USE_TILES
void blfilter_Tiles(float *paddedImg, float *gaussian_sp, int TileSize, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img);
#endif

_MM_ALIGN16 float getFilteredPixel(float* Ip_Tile, float* gaussian_sp, float sigmasq, int i, int j, int filter_hw, int rows, int cols, int TileSize);
CImg<float> writeCIMGFromArr(float* imgArr, int rows, int cols);

double my_difftime();
void *MyCalloc(int MemSize);
void *MyMalloc(int MemSize);
#if defined USE_TILES || EXP_TILES
void ReadTile(float *Tile, int Imgsize, int TileSize, int filter_hw , float *paddedImg, int tileId);
void WriteTile(float *Tile, int Cols, int filter_hw, int TileSize, float *Img, int tileId);
#endif

/* Main Function*/
int main(int argc, char* argv[])
{
    int i = 0, j = 0;
    float *gaussian_sp;
    float *Ip_Img, *Op_Img, *paddedImg;

#if defined USE_TILES || EXP_TILES
    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s <Image File> <Filter Radius> <Sigma Spatial> <Sigma Photometric> <TileSize>\n", argv[0]);
        exit(0);
    }
    int TileSize = atoi(argv[5]);
#else
    if (argc != 5)
    {
        fprintf(stderr, "Usage: %s <Image File> <Filter Radius> <Sigma Spatial> <Sigma Photometric>\n", argv[0]);
        exit(0);
    }
#endif

    char* fileName = argv[1];//1. input image
    int filter_hw = atoi(argv[2]); // filter half width
    float sigma_sp = (float) atof(argv[3]); // sigma_sp (standard deviation for the spatial distance)
    float sigma_ph = (float) atof(argv[4]); // sigma_ph standard deviation for the photometric distance
    double start, end;/* To measure Time */

    /* Load image using CImg library */
    CImg<float> image(fileName);

    int cols = image.width();
    int rows = image.height();
    int ds = image.spectrum(); //Gray or RGB (or could be more channels);

#ifdef USE_TILES
    if ((cols % TileSize != 0) || (rows % TileSize !=0))
    {
        printf("Image width = %d or Height = %d is not divisible by TileSize\nexiting..",cols, rows);
        exit(1);
    }
#endif

    /* Check for erros */
    Ip_Img = (float *) _mm_malloc(rows * cols * sizeof(float),16);

    Op_Img = (float *) _mm_malloc(rows * cols * sizeof(float),16);

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
    gaussian_sp = (float *) _mm_malloc( (2 * filter_hw + 1) *  (2 * filter_hw + 1) * sizeof(float), 16);

    for( i = -filter_hw; i<= filter_hw ; i++)
    {
        for(j = -filter_hw; j<= filter_hw; j++)
        {
            gaussian_sp[INDEX((i+filter_hw), (j+filter_hw), (2 * filter_hw + 1))] = exp( -(j*j + i*i) / (2*sigma_sp*sigma_sp) );
        }
    }

    //zero pad the input image to accomodate filtering at the boundary
    paddedImg = (float *) _mm_malloc( (rows + 2 * filter_hw) * (cols + 2 * filter_hw) * sizeof(float), 16);
    memset(paddedImg, 0,((rows + 2 * filter_hw) * (cols + 2 * filter_hw) * sizeof(float))); 

    for( i = 0; i < rows; i++)
        memcpy(&(paddedImg[INDEX((i+filter_hw),filter_hw,(cols + filter_hw + filter_hw))]), &(Ip_Img[INDEX(i,0,cols)]), cols * sizeof(float));

    //free(Ip_Img);

    start = my_difftime();
#ifdef USE_TILES
    blfilter_Tiles(paddedImg, gaussian_sp, TileSize, rows, cols, filter_hw, sigma_sp, sigma_ph, Op_Img);
#elif EXP_TILES
    blfilter_Read_Tiles(paddedImg, gaussian_sp, TileSize, rows, cols, filter_hw, sigma_sp, sigma_ph, Op_Img);
#else
    blfilter(paddedImg, gaussian_sp, rows, cols, filter_hw, sigma_sp, sigma_ph, Op_Img);
#endif
    end = my_difftime() - start;
    printf("Time taken for filtering = %6.2f \n", end);
    for(i = 0; i < rows; i++)
        for(j  = 0; j < cols; j++)
        {
#if defined USE_TILES || EXP_TILES
            filteredImg(i,j,0,0) = paddedImg[INDEX(i,j,(cols+2*filter_hw))];
#else
            filteredImg(i,j,0,0) = Op_Img[INDEX(i,j,cols)];
#endif
        }

/*	free(paddedImg);
	free(Op_Img);
	free(gaussian_sp);
*/
    CImgDisplay gray_disp(grayImg, "Gray Image");
    CImgDisplay flt_disp(filteredImg, "Filtered Image");

    //create event loop to display the image
      while (!gray_disp.is_closed() ) {
        gray_disp.wait();
        //flt_disp.wait();
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
#pragma omp parallel shared(paddedImg, gaussian_sp, Rows, Cols, filter_hw, sigma_ph)
    {
#pragma omp for private( filter_current_pixel,i,j,k,l,normal_factor,gaussian_ph, gaussian_bl, filtered_pixel,pd) 
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

                    gaussian_ph = exp( pd * sigmasq );
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






/* Funtion to Apply Bilateral filtering using Tiles without reading tiles 
* ----------------------------------------------------------------------*
* ----------------------------------------------------------------------*
* ----------------------------------------------------------------------*
* ----------------------------------------------------------------------*/

#ifdef USE_TILES
void blfilter_Tiles(float *paddedImg, float *gaussian_sp, int TileSize, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img)
{
    float filter_current_pixel = 0; //local window
    float pd = 0, pixel_value=0;          //photometric distance
    float gaussian_ph = 0; //photometric filter coefficient
    float gaussian_bl = 0; //bilateral filter coefficient
    float filtered_pixel = 0; //filtered value of pixel
    float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient
    int i = 0,j = 0, k = 0, l = 0;

    //image height and width
    int Rows = rows;
    int Cols = cols;

    int tiles;
    /* A tile has tile size full of pixels from input image and every other pixel needed to calculate filtered tile */
    int numTiles = (Rows * Cols) / (TileSize * TileSize);
    float sigmasq  = 1/(-2*sigma_ph*sigma_ph);	
#ifdef USE_OMP
#pragma omp parallel shared(paddedImg, gaussian_sp, Rows, Cols, TileSize, filter_hw, sigma_ph)
        {
	#pragma omp for private( filter_current_pixel,i,j,k,l,normal_factor,gaussian_ph, gaussian_bl, filtered_pixel,pd) 
#	endif
    for( tiles = 0; tiles < numTiles; tiles++)
    {

        int TilesPerRow = Cols/TileSize;
        int Tile_row = floor(tiles / TilesPerRow);
        int Tile_col = tiles % TilesPerRow;

        for( i = filter_hw; i < TileSize + filter_hw ; i++)
        {
            /* Get the first pixel in current row  */
            int ImageIndex = (Tile_row * TileSize + i) * (Cols + 2 * filter_hw) + Tile_col * TileSize ;
            int OpImgIndex = (Tile_row * TileSize + (i-filter_hw)) * (Cols + 2*filter_hw) + Tile_col * TileSize;
            for( j = filter_hw; j < TileSize + filter_hw; j++)
            {
                filtered_pixel = 0;
                normal_factor  = 0;
                pixel_value = paddedImg[(ImageIndex + j) ];
                for( k = -filter_hw ; k <= filter_hw ; k++)
                {
                    for( l = -filter_hw ; l <= filter_hw ; l++)
                    {
                        filter_current_pixel = paddedImg[ImageIndex +  INDEX(k, (j + l), (Cols+2*filter_hw))];
                        //photometric distance = difference in pixel values (squared)
                        pd = pixel_value - filter_current_pixel;
                        pd = pd * pd;
                        gaussian_ph = exp(pd * sigmasq);
                        gaussian_bl = gaussian_ph * gaussian_sp[INDEX((k+filter_hw), (l+filter_hw), (2*filter_hw + 1))] ;
                        filtered_pixel += gaussian_bl * filter_current_pixel;
                        normal_factor += gaussian_bl;
                    }
                }
                filtered_pixel = filtered_pixel / normal_factor;
                /* Write out the new value to output image */
                paddedImg[(OpImgIndex + j-filter_hw)] = filtered_pixel;
               
            }
        }
    }
#ifdef USE_OMP
    }
#endif
}
#endif







/* Funtion to Apply Bilateral filtering using Tiles without reading tiles 
* ----------------------------------------------------------------------*
* ----------------------------------------------------------------------*
* ----------------------------------------------------------------------*
* ----------------------------------------------------------------------*/

#ifdef EXP_TILES
void blfilter_Read_Tiles(float *paddedImg, float *gaussian_sp, int TileSize, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img)
{
    float filter_current_pixel = 0; //local window
    float pd = 0, pixel_value=0;          //photometric distance
    float gaussian_ph = 0; //photometric filter coefficient
    float gaussian_bl = 0; //bilateral filter coefficient
    float filtered_pixel = 0; //filtered value of pixel
    float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient
    int i = 0,j = 0, k = 0, l = 0;

    //image height and width
    int Rows = rows;
    int Cols = cols;

    float *Ip_Tile, *Op_Tile;
    int tiles;
    /* A tile has tile size full of pixels from input image and every other pixel needed to calculate filtered tile */
    int numTiles = (Rows * Cols) / (TileSize * TileSize);
    float sigmasq  = 1/(-2*sigma_ph*sigma_ph);		
#ifdef USE_OMP
#pragma omp parallel  private(Ip_Tile, Op_Tile) shared(paddedImg, gaussian_sp, Rows, Cols, TileSize, filter_hw, sigma_ph)
        {
#endif
    Ip_Tile = (float *) _mm_malloc( (TileSize + 2 * filter_hw ) * ( TileSize + 2 * filter_hw) * sizeof(float), 16);
    Op_Tile = (float *) _mm_malloc( TileSize  * TileSize  * sizeof(float), 16);
#ifdef USE_OMP
#pragma omp for private(tiles, filter_current_pixel,i,j,k,l,normal_factor,gaussian_ph, gaussian_bl, filtered_pixel,pd) 
#endif
    for( tiles = 0; tiles < numTiles; tiles++)
    {
        ReadTile(Ip_Tile, Cols, TileSize, filter_hw , paddedImg, tiles);
        for( i = filter_hw; i < TileSize + filter_hw ; i++)
        {
            for( j = filter_hw; j < TileSize + filter_hw; j++)
            {
		
                filtered_pixel = getFilteredPixel(Ip_Tile, gaussian_sp, sigmasq, i, j, filter_hw, rows, cols, TileSize);
                Op_Tile[INDEX((i-filter_hw), (j-filter_hw), TileSize)] = filtered_pixel;
 
            }
        }
        WriteTile(Op_Tile, Cols, filter_hw, TileSize, paddedImg, tiles);
	
	/** To visualize how the tiles are written in image
	CImg<float> partialImg = writeCIMGFromArr(Op_Img, rows, cols);
	CImgDisplay flt_disp(partialImg, "Filtered Image");		

	while (!flt_disp.is_closed() ) {
	flt_disp.wait();  
        }
	**/	

    }
#ifdef USE_OMP
    }
#endif
}
#endif

CImg<float> writeCIMGFromArr(float* imgArr, int rows, int cols) {
	CImg<float> img(cols, rows,1,1,0);

	for(int i = 0; i < rows; i++)
		for(int j = 0; j < cols; j++) {
			img(i, j) = imgArr[INDEX(i,j,cols)];
	}

	return img;

}

_MM_ALIGN16 float getFilteredPixel(float* Ip_Tile, float* gaussian_sp, float sigmasq, int i, int j, int filter_hw, int rows, int cols, int TileSize) {

	__declspec(align(16)) float local_pixel0 = 0; //local window
	__declspec(align(16)) float local_pixel1 = 0;
	__declspec(align(16)) float local_pixel2 = 0;
	__declspec(align(16)) float local_pixel3 = 0;

	__declspec(align(16)) float pd0 = 0; //photometric distance
	__declspec(align(16)) float pd1 = 0;
	__declspec(align(16)) float pd2 = 0;
	__declspec(align(16)) float pd3 = 0;

	__declspec(align(16)) float gaussian_ph0 = 0; //photometric filter coefficient
	__declspec(align(16)) float gaussian_ph1 = 0;
	__declspec(align(16)) float gaussian_ph2 = 0;
	__declspec(align(16)) float gaussian_ph3 = 0;

	__declspec(align(16)) float gaussian_bl0 = 0; //bilateral filter coefficient
	__declspec(align(16)) float gaussian_bl1 = 0;
	__declspec(align(16)) float gaussian_bl2 = 0;
	__declspec(align(16)) float gaussian_bl3 = 0;

	__declspec(align(16)) float filtered_pixel0 = 0; //filtered value of pixel 
	__declspec(align(16)) float filtered_pixel1 = 0;
	__declspec(align(16)) float filtered_pixel2 = 0;
	__declspec(align(16)) float filtered_pixel3 = 0;	

	__declspec(align(16)) float filtered_pixel = 0;

	__declspec(align(16)) float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient

	_MM_ALIGN16 float spk0 = 0, spk1 = 0, spk2 = 0, spk3 = 0; 

	
	__declspec(align(16)) float pixel_value = Ip_Tile[INDEX( i, j, (TileSize + filter_hw +  filter_hw))];

		//1. calculate the photometric gaussian distances and the filter coefficent (no need to create a full kernel)	
		//2. multiply each pixel of the window with the spatial and photometric coefficients

		//#pragma ivdep
		for(int k = -filter_hw ; k <= filter_hw ; k+=2)
			for(int l = -filter_hw ; l <= filter_hw ; l+=2) {
			
			//local_window[k + filter_hw][l + filter_hw] = paddedImg(i + k, j + l);
						
			local_pixel0 = Ip_Tile[INDEX((i+k), (j+l), (TileSize + filter_hw + filter_hw))];
			local_pixel1 = Ip_Tile[INDEX((i+k), (j+l+1), (TileSize + filter_hw + filter_hw))];
			local_pixel2 = Ip_Tile[INDEX((i+k+1), (j+l), (TileSize + filter_hw + filter_hw))];
			local_pixel3 = Ip_Tile[INDEX((i+k+1), (j+l+1), (TileSize + filter_hw + filter_hw))];

			spk0 = gaussian_sp[INDEX((k+filter_hw), (l+filter_hw), (2*filter_hw + 1))];
			spk1 = gaussian_sp[INDEX((k+filter_hw), (l+1+filter_hw), (2*filter_hw + 1))];
			spk2 = gaussian_sp[INDEX((k+1+filter_hw), (l+filter_hw), (2*filter_hw + 1))];
			spk3 = gaussian_sp[INDEX((k+1+filter_hw), (l+1+filter_hw), (2*filter_hw + 1))];
			
			//photometric distance = difference in pixel values (squared)
			pd0 = pixel_value - local_pixel0;
			pd1 = pixel_value - local_pixel1;
			pd2 = pixel_value - local_pixel2;
			pd3 = pixel_value - local_pixel3;

			pd0 = pd0 * pd0;  
			pd1 = pd1 * pd1;  
			pd2 = pd2 * pd2;  
			pd3 = pd3 * pd3;  

			gaussian_ph0 = exp( pd0 * sigmasq );
			gaussian_ph1 = exp( pd1 * sigmasq );
			gaussian_ph2 = exp( pd2 * sigmasq );
			gaussian_ph3 = exp( pd3 * sigmasq );			
	
			gaussian_bl0 = gaussian_ph0 * spk0;
			gaussian_bl1 = gaussian_ph1 * spk1;
			gaussian_bl2 = gaussian_ph2 * spk2;
			gaussian_bl3 = gaussian_ph3 * spk3;

			if(k == filter_hw){gaussian_bl2 = 0; gaussian_bl3 = 0;}
			if(l == filter_hw){gaussian_bl1 = 0; gaussian_bl3 = 0;}

			filtered_pixel0 += gaussian_bl0 * local_pixel0;	
			filtered_pixel1 += gaussian_bl1 * local_pixel1;
			filtered_pixel2 += gaussian_bl2 * local_pixel2;	
			filtered_pixel3 += gaussian_bl3 * local_pixel3;

			normal_factor += (gaussian_bl0 + gaussian_bl1 + gaussian_bl2 + gaussian_bl3);	
							
		} 

		filtered_pixel = (filtered_pixel0 + filtered_pixel1 + filtered_pixel1 + filtered_pixel3) / normal_factor;

	return filtered_pixel;


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
    void *ptr;

    if ( (ptr = (void *) malloc(MemSize)) == (void *) NULL)
    {
        fprintf(stderr, "Error (MyMalloc) : Unable to allocate memory\n");
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
        fprintf(stderr, "Error (MyCalloc) : Unable to allocate memory\n");
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


void ReadTile(float *Tile, int Cols, int TileSize, int filter_hw , float *paddedImg, int tileId)
{
    int i,j;
    int TilesPerRow, row, col;
    int ImageIndex = 0;
    int TileCols = TileSize + 2 * filter_hw;
    /* Find Number of tiles - depends on Ip_Img not padded Image */
    TilesPerRow = Cols/TileSize;
    row = floor(tileId / TilesPerRow);
    col = tileId % TilesPerRow;
    ImageIndex = (row * TileSize) * (Cols + 2 * filter_hw) + col * TileSize ;
    for (i = 0; i < TileCols; i++)
    {
        memcpy(&(Tile[(i* TileCols)]), &(paddedImg[(ImageIndex + (i * (Cols + 2 * filter_hw )) )]), (TileCols * sizeof(float)));
    }
    /*
    CImg<float> tmp(TileCols, TileCols, 1, 1, 0); //create a single image dz = 1, rgb channel  = 1 (only gray); fill it with 0s.
    for (i = 0; i< TileCols; i++)
        for (j = 0; j< TileCols; j++)
        {
            tmp(i,j,0,0) = Tile[INDEX(i,j,TileCols)];
            //printf("%f\n", Tile[ImageIndex-1]);
        }

    CImgDisplay gray_disp(tmp, "I/p Tile");
    while (!gray_disp.is_closed() ) {
        gray_disp.wait();
        //flt_disp.wait();
    }
    */
}

void WriteTile(float *Tile, int Cols, int filter_hw, int TileSize, float *Img, int tileId)
{
    int i,j;
    int TilesPerRow, row, col;
    int ImageIndex;
    TilesPerRow = Cols / TileSize;
    row = floor(tileId / TilesPerRow);
    col = tileId % TilesPerRow;
    ImageIndex = row * TileSize * (Cols+2*filter_hw) + col * TileSize;
    for (i = 0; i < TileSize; i++)
        memcpy(&(Img[(ImageIndex + i * (Cols+2*filter_hw))]), &(Tile[(i*TileSize)]), TileSize * sizeof(float));
    /*
    CImg<float> tmp(TileSize, TileSize, 1, 1, 0); //create a single image dz = 1, rgb channel  = 1 (only gray); fill it with 0s.
    for (i = 0; i< TileSize; i++)
        for (j = 0; j< TileSize; j++)
        {
            tmp(i,j,0,0) = Tile[INDEX(i,j,TileSize)];
            //printf("%f\n", Tile[ImageIndex-1]);
        }

    CImgDisplay gray_disp(tmp, "O/p Tile");
    while (!gray_disp.is_closed() ) {
        gray_disp.wait();
    }
    */
}













