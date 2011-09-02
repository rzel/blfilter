
#include "CImg.h"

using namespace cimg_library;
CImg<double>  blfilter(CImg<double>, int, double, double);


int main(int argc, char* args[]) {


//1. input image
    char* fileName = args[1];
//2. filter half width
    int filter_hw = (int)atoi(args[2]);
//3. sigma_sp (standard deviation for the spatial distance)
    double sigma_sp = (double)atoi(args[3]);
//4. sigma_ph standard deviation for the photometric distance
    double sigma_ph = (double)atoi(args[4]);

//load image
    CImg<double> image(fileName);

//load image
//CImg<double> image("lena.png");

    int dx = image.width();
    int dy = image.height();
    int ds = image.spectrum(); //Gray or RGB (or could be more channels);

//create a single image dz = 1, rgb channel  = 1 (only gray); fill it with 0s.
    CImg<double> grayImg(dx,dy,1,1,0);

    if(ds == 3) { //convert the input RGB to gray image. (work with gray image for now)
        for(int i = 0; i < dx; i++)
            for(int j  = 0; j < dy; j++) {
                double gray_pixel = 0.299 * image(i,j,0,0) +
                                    0.587 * image(i,j,0,1) +
                                    0.114 * image(i,j,0,2);
                //printf("(%d,%d) = (%f) \n", i,j, gray_pixel);
                grayImg(i,j,0,0) = gray_pixel;
            }
    }
    else
        grayImg = image;

    printf("spectra = %d\n\n", ds );
    printf("filtering ... \n\n");

// this function returns the filtered image
    clock_t start, end;
    double elapsed;

    start = clock();
//call the filtering function
    CImg<double> filteredImg = blfilter(grayImg, filter_hw, sigma_sp, sigma_ph);
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Total CPU cycles = %ld\n", (end - start));
    printf("Time taken for filtering = %f \n", elapsed);

//create display windows
//CImgDisplay main_disp(image, "Lena");
    CImgDisplay gray_disp(grayImg, "Gray Image");
    CImgDisplay flt_disp(filteredImg, "Filtered Image");

//create event loop to display the image
    while (!gray_disp.is_closed() ) {
        gray_disp.wait();
        //flt_disp.wait();
    }

}


CImg<double>  blfilter(CImg<double> grayImg, int filter_hw, double sigma_sp, double sigma_ph)
{
    //image height and width
    int dx = grayImg.width();
    int dy = grayImg.height();


    //create a distance grid(x and y distance from the filter core)
    //(dont need separate grids as in matlab)
    //double xgrid[2*filter_hw + 1][2*filter_hw + 1];
    //double ygrid[2*filter_hw + 1][2*filter_hw + 1];

    //gaussian kernel for the spatial distances
    double gaussian_sp[2*filter_hw + 1][2*filter_hw + 1];

    for(int i = -filter_hw; i<= filter_hw ; i++)
        for(int j = -filter_hw; j<= filter_hw; j++) {
            //xgrid[i+filter_hw][j+filter_hw] = j;
            //ygrid[i+filter_hw][j+filter_hw] = i;
            gaussian_sp[i+filter_hw][j+filter_hw] = exp( -(j*j + i*i) / (2*sigma_sp*sigma_sp) );

        }


    /**
    //print the kernel for splot (to visualize the kernel using gnuplot)
    for(int i = -filter_hw; i<= filter_hw ; i++) {
    	for(int j = -filter_hw; j<= filter_hw; j++) {
    		printf(" %d\t%d\t%f \n ", i, j,gaussian_sp[i+filter_hw][j+filter_hw]);
    	}
    }
    **/

    //create empty image for filtered output dz = 1, rgb channel  = 1 (only gray); fill it with 0s.
    CImg<double> filteredImg(dx,dy,1,1,0);

    //zero pad the input image to accomodate filtering at the boundary
    CImg<double> paddedImg(dx + 2 * filter_hw, dy + 2 * filter_hw, 1, 1, 0);

    for(int i = 0; i < dx; i++)
        for(int j = 0; j < dy; j++)
            paddedImg(i+filter_hw, j+filter_hw) = grayImg(i,j);

    /**
    CImgDisplay pad_disp(paddedImg, "padded image");
    //create event loop to display the image
    while (!pad_disp.is_closed()) {
    	pad_disp.wait();
    	//fil_disp.wait();
        }
    **/

    double local_window[2*filter_hw + 1][2*filter_hw + 1];
    double local_pixel = 0; //local window
    double pd = 0; //photometric distance
    double gaussian_ph = 0; //photometric filter coefficient
    double gaussian_bl = 0; //bilateral filter coefficient
    double filtered_pixel = 0; //filtered value of pixel
    double normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient

    for(int i = filter_hw; i < dx + filter_hw ; i++)
        for(int j = filter_hw; j < dy + filter_hw; j++) {

            filtered_pixel = 0;
            normal_factor = 0;

            //1. calculate the photometric gaussian distances and the filter coefficent (no need to create a full kernel)
            //2. multiply each pixel of the window with the spatial and photometric coefficients

            for(int k = -filter_hw ; k <= filter_hw ; k++)
                for(int l = -filter_hw ; l <= filter_hw ; l++) {

                    //local_window[k + filtered_hw][l + filtered_hw] = paddedImg(i + k, j + l);
                    local_pixel = paddedImg(i + k, j + l);

                    //photometric distance = difference in pixel values (squared)
                    pd = paddedImg(i,j) - local_pixel;
                    pd = pd * pd;

                    gaussian_ph = exp( -(pd) / (2*sigma_ph*sigma_ph) );
                    gaussian_bl = gaussian_ph * gaussian_sp[k+filter_hw][l+filter_hw] ;

                    filtered_pixel += gaussian_bl * local_pixel;
                    normal_factor += gaussian_bl;

                }

            filtered_pixel = filtered_pixel / normal_factor;
            //printf("(%d, %d) = %f", i, j, normal_factor);
            filteredImg(i-filter_hw, j - filter_hw) =  filtered_pixel;

        }



    return filteredImg;

}









