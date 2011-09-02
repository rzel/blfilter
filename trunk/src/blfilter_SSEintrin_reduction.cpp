
#include "CImg.h"
#include <omp.h>
#include<pmmintrin.h> //SSE3
#include<xmmintrin.h>
#include<ia32intrin.h> //general intrinsics


using namespace cimg_library;
CImg<float>  blfilter(CImg<float>, int, float, float);
float* intrinadd(float* pa, float* pb);


int main(int argc, char* args[]) {


//1. input image
    char* fileName = args[1];
//2. filter half width
    int filter_hw = (int)atoi(args[2]);
//3. sigma_sp (standard deviation for the spatial distance)
    float sigma_sp = (float)atoi(args[3]);
//4. sigma_ph standard deviation for the photometric distance
    float sigma_ph = (float)atoi(args[4]);

//load image
    CImg<float> image(fileName);

//load image
//CImg<float> image("lena.png");

    int dx = image.width();
    int dy = image.height();
    int ds = image.spectrum(); //Gray or RGB (or could be more channels);

//create a single image dz = 1, rgb channel  = 1 (only gray); fill it with 0s.
    CImg<float> grayImg(dx,dy,1,1,0);

    if(ds == 3) { //convert the input RGB to gray image. (work with gray image for now)
        for(int i = 0; i < dx; i++)
            for(int j  = 0; j < dy; j++) {
                float gray_pixel = 0.299 * image(i,j,0,0) +
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


    clock_t start, end;
    float elapsed;
    start = clock();

//__int64 start, stop, elapsed1, elapsed2;
//start= _rdtsc();

//call the filtering function

    CImg<float> filteredImg = blfilter(grayImg, filter_hw, sigma_sp, sigma_ph);
    end = clock();
    elapsed = ((float) (end - start)) / CLOCKS_PER_SEC;
//stop= _rdtsc();

    printf("Total CPU cycles = %ld\n", (end - start));
    printf("Time taken for filtering = %f \n", elapsed);
//printf("Processor cycles = %I64u \n",  (stop - start));

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


CImg<float>  blfilter(CImg<float> grayImg, int filter_hw, float sigma_sp, float sigma_ph)
{
    //image height and width
    int dx = grayImg.width();
    int dy = grayImg.height();

    //create a distance grid(x and y distance from the filter core)
    //(dont need separate grids as in matlab)
    //float xgrid[2*filter_hw + 1][2*filter_hw + 1];
    //float ygrid[2*filter_hw + 1][2*filter_hw + 1];

    //gaussian kernel for the spatial distances
    float gaussian_sp[2*filter_hw + 1][2*filter_hw + 1];

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
    CImg<float> filteredImg(dx,dy,1,1,0);

    //zero pad the input image to accomodate filtering at the boundary
    CImg<float> paddedImg(dx + 2 * filter_hw, dy + 2 * filter_hw, 1, 1, 0);

    for(int i = 0; i < dx; i++)
        for(int j = 0; j < dy; j++)
            paddedImg(i + filter_hw, j + filter_hw) = grayImg(i,j);

    /**
    CImgDisplay pad_disp(paddedImg, "padded image");
    //create event loop to display the image
    while (!pad_disp.is_closed()) {
    	pad_disp.wait();
    	//fil_disp.wait();
        }
    **/

    float local_window[2*filter_hw + 1][2*filter_hw + 1];

    float filtered_pixel = 0;

    float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient

    //gaussian window in 1-d array
    int sizewin = (2*filter_hw+1)*(2*filter_hw+1);
    int lwgcount = 0;
    float linwindowg[sizewin];
    for(int k = 0 ; k < 2*filter_hw+1 ; k++)
        for(int l = 0 ; l < 2*filter_hw+1 ; l++)
            linwindowg[lwgcount++] = gaussian_sp[k][l];

    //#pragma ivdep
    //#pragma omp parallel shared(filteredImg, paddedImg)
#ifdef USE_OMP
    #pragma omp parallel shared(filteredImg, paddedImg)
    {
    #pragma omp for
#endif
    for(int i = filter_hw; i < dx + filter_hw  ; i++)
        for(int j = filter_hw; j < dy + filter_hw; j++) {

            filtered_pixel = 0;
            normal_factor = 0;

            int lwcount = 0;

            __declspec(align(16)) float sigmasq  = -2*sigma_ph*sigma_ph;
            __declspec(align(16)) float linwindow[sizewin];

            for(int k = -filter_hw ; k <= filter_hw ; k++)
                for(int m = -filter_hw ; m <= filter_hw ; m++)  {
                    linwindow[lwcount++] = paddedImg(i + k, j + m);
                }


            /** intrinsic method 1 :  by aliasing the pointer with registers (gave strange segmentation fault)
            __declspec(align(16)) float* blfiltered = (float*) _mm_malloc(sizewin * sizeof(float), 16);
            __declspec(align(16)) float* normwin = (float*) _mm_malloc(sizewin * sizeof(float), 16);

            //__m128 pmg =  _mm_set1_ps(paddedImg(i,j));
            //__m128 sgmas = _mm_set1_ps(sigmasq);

            __declspec(align(16)) float blfiltered[sizewin]; //pointer to filtered pixels (to be used in intrinsic)
            __declspec(align(16)) float normwin[sizewin];

            __m128 *lp = (__m128*)linwindow;
            __m128 *spk = (__m128*)linwindowg;
            __m128 *pmg = (__m128*)pimg;
            __m128 *sgmas = (__m128*)sigmaphs;
            __m128 *gfp = (__m128*)blfiltered;
            __m128 *norm = (__m128*)normwin;

            for(int l = 0; l < lwcount ; l += 4) {

            	*gfp = _mm_sub_ps(*pmg, *lp);
            	*gfp = _mm_mul_ps(*gfp, *gfp); //photometric distance squared
            	*gfp = _mm_div_ps(*gfp, *sgmas); //exp power
            	*gfp = _mm_exp_ps(*gfp); //photometric coeff (function is from SVML - no corresponding Assembly is available).
                    *gfp = _mm_mul_ps(*gfp, *spk); //blfilter coeff
            	*gfp = _mm_mul_ps(*gfp, *lp); //filtered coefficient

            	lp++; gfp++; spk++;
            }

            for(int q = 0; q < lwcount ; q++) {
            	filtered_pixel +=  blfiltered[q];
            	normal_factor += normwin[q];
            }
            **/


            //Method 2. intrinisics by using the load / store  command

            __declspec(align(16)) float blfiltered[sizewin]; //pointer to filtered pixels (to be used in intrinsic)
            __declspec(align(16)) float normwin[sizewin];

            __m128 lp, pmg, sgmas, gfp, spk;
            pmg =  _mm_set1_ps(paddedImg(i,j));
            sgmas = _mm_set1_ps(sigmasq);

            /**stop undoing **/
            //__declspec(ali
            __m128 o1, o2; //final aggregate filtered and normal factor vectors
            o1 = _mm_setzero_ps();
            o2 = _mm_setzero_ps();

            for (int l = 0; l < lwcount; l+= 4) {
                lp = _mm_load_ps(linwindow + l);
                spk = _mm_load_ps(linwindowg + l);
                gfp = _mm_sub_ps(pmg, lp);
                gfp = _mm_mul_ps(gfp, gfp);
                gfp = _mm_div_ps(gfp, sgmas);
                gfp = _mm_exp_ps(gfp); //photometric filter coeff
                gfp = _mm_mul_ps(gfp, spk); //bilateral filter coeff

                o1 = _mm_hadd_ps(o1, gfp);	//normalization coefficient using reduction
                gfp = _mm_mul_ps(gfp, lp); //multiply it with the local window (filtered coeff)

                o2 = _mm_hadd_ps(o2, gfp);	//filtered coefficient using reduction
            }

            /**
            // : Use Reduction Here
            for(int q = 0; q < lwcount ; q++) {
            	filtered_pixel +=  blfiltered[q];
            	normal_factor += normwin[q];
            }**/

            float np[4], fp[4];
            _mm_store_ps(np, o1);
            _mm_store_ps(fp, o2);

            filtered_pixel = fp[0];
            normal_factor = np[0];

            filtered_pixel = filtered_pixel / normal_factor;
            filteredImg(i - filter_hw, j - filter_hw) =  filtered_pixel;

            //printf("(%d, %d) = %f\t%f\n", (i - filter_hw), (j - filter_hw), normal_factor, filtered_pixel);
        }
#ifdef USE_OMP
    }
#endif

    return filteredImg;

}


float* intrinadd(float* pa, float* pb) {

    float pc[4];

    __m128 *t0 = (__m128*)pa;
    __m128 *t1 = (__m128*)pb;
    __m128 *t3 = (__m128*)pc;
    *t3 = _mm_add_ps(*t0, *t1);
    *t3 = _mm_exp_ps(*t3);

    return pc;
}








