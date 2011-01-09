#include "def.h"
#include <sys/time.h>
#include "CImg.h"
using namespace cimg_library;

/*Function Declarations */
void blfilter(float *Ip_Img, int width, int height, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img);
#ifdef USE_TILES
void blfilter_Tiles(float *Ip_Img, int TileSize, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img);
#endif
double my_difftime();
void *MyCalloc(int MemSize);
void *MyMalloc(int MemSize);
#ifdef USE_TILES
void ReadTile(float *Tile, int Imgsize, int TileSize, int filter_hw , float *paddedImg, int tileId);
void WriteTile(float *Tile, int Cols, int TileSize, float *Img, int tileId);
#endif

/* Main Function*/
int main(int argc, char* argv[])
{
    int i = 0, j = 0;

#ifdef USE_TILES
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

    float *Ip_Img, *Op_Img;

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
                Ip_Img[INDEX(i,j,cols)] = image(i,j);
            }
        grayImg = image;
    }

    start = my_difftime();
#ifdef USE_TILES
    blfilter_Tiles(Ip_Img, TileSize, rows, cols, filter_hw, sigma_sp, sigma_ph, Op_Img);
#else
    blfilter(Ip_Img, rows, cols, filter_hw, sigma_sp, sigma_ph, Op_Img);
#endif
    end = my_difftime() - start;
    printf("Time taken for filtering = %6.2f \n", end);
    for(i = 0; i < rows; i++)
        for(j  = 0; j < cols; j++)
        {
            filteredImg(i,j) = Op_Img[INDEX(i,j,cols)];
        }

    CImgDisplay gray_disp(grayImg, "Gray Image");
    CImgDisplay flt_disp(filteredImg, "Filtered Image");

    //create event loop to display the image
    while (!gray_disp.is_closed() ) {
        gray_disp.wait();
        //flt_disp.wait();
    }
    return 0;
}

void blfilter(float *Ip_Img, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img)
{
    float filter_current_pixel = 0; //local window
    float pd = 0;          //photometric distance
    float gaussian_ph = 0; //photometric filter coefficient
    float gaussian_bl = 0; //bilateral filter coefficient
    float filtered_pixel = 0; //filtered value of pixel
    float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient
    float *paddedImg;
    int i = 0,j = 0, k = 0, l = 0;

    //image height and width
    int Rows = rows;
    int Cols = cols;

    //gaussian kernel for the spatial distances
    float *gaussian_sp;
    gaussian_sp = (float *) MyMalloc( (2 * filter_hw + 1) *  (2 * filter_hw + 1) * sizeof(float));

    for( i = -filter_hw; i<= filter_hw ; i++)
    {
        for(j = -filter_hw; j<= filter_hw; j++)
        {
            gaussian_sp[INDEX((i+filter_hw), (j+filter_hw), (2 * filter_hw + 1))] = exp( -(j*j + i*i) / (2*sigma_sp*sigma_sp) );
        }
    }

    //zero pad the input image to accomodate filtering at the boundary
    paddedImg = (float *) MyCalloc( (Rows + 2 * filter_hw) * (Cols + 2 * filter_hw) * sizeof(float));

    for( i = 0; i < Rows; i++)
        memcpy(&(paddedImg[INDEX((i+filter_hw),filter_hw,(Cols + filter_hw + filter_hw))]), &(Ip_Img[INDEX(i,0,Cols)]), Cols * sizeof(float));

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

                        gaussian_ph = exp( -(pd) / (2*sigma_ph*sigma_ph) );
                        gaussian_bl = gaussian_ph * gaussian_sp[INDEX((k+filter_hw), (l+filter_hw), (2*filter_hw + 1))] ;

                        filtered_pixel += gaussian_bl * filter_current_pixel;
                        normal_factor += gaussian_bl;
                    }
                }
                filtered_pixel = filtered_pixel / normal_factor;
                Op_Img[INDEX((i-filter_hw), (j-filter_hw), Cols)] = filtered_pixel;
                //filteredImg(i-filter_hw, j-filter_hw) =  filtered_pixel;
            }
        }
}

#ifdef USE_TILES
void blfilter_Tiles(float *Ip_Img, int TileSize, int rows, int cols, int filter_hw, float sigma_sp, float sigma_ph, float *Op_Img)
{
    float filter_current_pixel = 0; //local window
    float pd = 0;          //photometric distance
    float gaussian_ph = 0; //photometric filter coefficient
    float gaussian_bl = 0; //bilateral filter coefficient
    float filtered_pixel = 0; //filtered value of pixel
    float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient
    float *paddedImg;
    int i = 0,j = 0, k = 0, l = 0;

    //image height and width
    int Rows = rows;
    int Cols = cols;

    float *Ip_Tile, *Op_Tile;
    int tiles;
    /* A tile has tile size full of pixels from input image and every other pixel needed to calculate filtered tile */
    int numTiles = (Rows * Cols) / (TileSize * TileSize);
    Ip_Tile = (float *) MyMalloc( (TileSize + 2 * filter_hw ) * ( TileSize + 2 * filter_hw) * sizeof(float));
    Op_Tile = (float *) MyMalloc( TileSize  * TileSize  * sizeof(float));
    //gaussian kernel for the spatial distances
    float *gaussian_sp;
    gaussian_sp = (float *) MyMalloc( (2 * filter_hw + 1) *  (2 * filter_hw + 1) * sizeof(float));

    for( i = -filter_hw; i<= filter_hw ; i++)
    {
        for(j = -filter_hw; j<= filter_hw; j++)
        {
            gaussian_sp[INDEX((i+filter_hw), (j+filter_hw), (2 * filter_hw + 1))] = exp( -(j*j + i*i) / (2*sigma_sp*sigma_sp) );
        }
    }

    //zero pad the input image to accomodate filtering at the boundary
    paddedImg = (float *) MyCalloc( (Rows + 2 * filter_hw) * (Cols + 2 * filter_hw) * sizeof(float));

    for( i = 0; i < Rows; i++)
        memcpy(&(paddedImg[INDEX((i+filter_hw),filter_hw,(Cols + filter_hw + filter_hw))]), &(Ip_Img[INDEX(i,0,Cols)]), Cols * sizeof(float));

    for( tiles = 0; tiles < numTiles; tiles++)
    {
        ReadTile(Ip_Tile, Cols, TileSize, filter_hw , paddedImg, tiles);

        for( i = filter_hw; i < TileSize + filter_hw ; i++)
        {
            for( j = filter_hw; j < TileSize + filter_hw; j++)
            {
                filtered_pixel = 0;
                normal_factor  = 0;
                pd = Ip_Tile[INDEX( i, j, (TileSize + filter_hw +  filter_hw))];

                for( k = -filter_hw ; k <= filter_hw ; k++)
                {
                    for( l = -filter_hw ; l <= filter_hw ; l++)
                    {
                        filter_current_pixel = Ip_Tile[INDEX((i + k), (j + l), (TileSize + filter_hw + filter_hw))];

                        //photometric distance = difference in pixel values (squared)
                        pd = pd - filter_current_pixel;
                        pd = pd * pd;

                        gaussian_ph = exp( -(pd) / (2*sigma_ph*sigma_ph) );
                        gaussian_bl = gaussian_ph * gaussian_sp[INDEX((k+filter_hw), (l+filter_hw), (2*filter_hw + 1))] ;

                        filtered_pixel += gaussian_bl * filter_current_pixel;
                        normal_factor += gaussian_bl;
                    }
                }
                filtered_pixel = filtered_pixel / normal_factor;
                Op_Tile[INDEX((i-filter_hw), (j-filter_hw), TileSize)] = filtered_pixel;
                //filteredImg(i-filter_hw, j-filter_hw) =  filtered_pixel;
            }
        }
        WriteTile(Op_Tile, Cols, TileSize,Op_Img, tiles);
    }
}
#endif


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
    int i;
    int TilesPerRow, row, col;
    int TileCols = TileSize + 2 * filter_hw;
    TilesPerRow = Cols / TileSize;
    row = floor(tileId / TilesPerRow);
    col = tileId % TilesPerRow;
    for (i = 0; i < TileCols; i++)
       {
        memcpy(&(Tile[INDEX(i, 0, TileCols)]), &(paddedImg[INDEX(((row * filter_hw) + i), (col * filter_hw), Cols) ]), TileCols * sizeof(float));
       }
}

void WriteTile(float *Tile, int Cols, int TileSize, float *Img, int tileId)
{
    int i;
    int TilesPerRow, row, col;
    TilesPerRow = Cols / TileSize;
    row = floor(tileId / TilesPerRow);
    col = tileId % TilesPerRow;
    for (i = 0; i < TileSize; i++)
        memcpy(&(Img[INDEX((i + row * TileSize), (col * TileSize),Cols)]), &(Tile[INDEX(i,0, TileSize)]), TileSize * sizeof(float));
}













