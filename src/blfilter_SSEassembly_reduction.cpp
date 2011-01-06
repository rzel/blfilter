
#include "CImg.h"
#include <omp.h>
#include<pmmintrin.h> //SSE3
#include<xmmintrin.h>
#include<ia32intrin.h> //general intrinsics
//#include<smmintrin.h> //extract


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
    for(int i = filter_hw; i < dx + filter_hw  ; i++)
        for(int j = filter_hw; j < dy + filter_hw; j++) {

            filtered_pixel = 0;
            normal_factor = 0;

            int lwcount = 0;
            __declspec(align(16)) float linwindow[sizewin];

            __declspec(align(16)) float pimg[] = {paddedImg(i,j), paddedImg(i,j), paddedImg(i,j), paddedImg(i,j)};
            __declspec(align(16)) float sigmasq  = -2*sigma_ph*sigma_ph;
            __declspec(align(16)) float sigmaphs[] = {sigmasq, sigmasq, sigmasq, sigmasq};

            for(int k = -filter_hw ; k <= filter_hw ; k++)
                for(int m = -filter_hw ; m <= filter_hw ; m++)  {
                    linwindow[lwcount++] = paddedImg(i + k, j + m);
                }

            __declspec(align(16)) float blfiltered[sizewin]; //pointer to filtered pixels (to be used in intrinsic)
            //float* blfiltered = (float*) _mm_malloc(sizewin * sizeof(float), 16);
            //__assume_aligned(blfiltered, 16);


            struct Vector4 {
                float x, y, z, w;
            };
            struct Vector4 o1, o2; //to store final filtered and normal factor vectors

            __m128 norm128;
            __m128 blfiltered128;
            //__m128 zero128 = _mm_setzero_ps(); //initialize empty vector (to hold the result of filtered and normal coefficients)

            __m128 pmgc128; //to hold the photometric gaussian coefficient
            __m128 lwin128; //for local window
            __m128 lwing128; //for spatial gaussian

            __m128 sigmaphs128 = _mm_set1_ps(sigmasq);
            __m128 pimg128 = _mm_set1_ps(paddedImg(i,j));


            float normwin[sizewin];

            //clear the outputs o1 and o2
            o1.x = o1.y = o1.z = o1.w = 0;
            o2.x = o2.y = o2.z = o2.w = 0;
            /**
            __asm__(
            "xorps %%xmm5, %%xmm5 ;"  corresponding intrinsic --> _mm_setzero_ps()
            "xorps %%xmm6, %%xmm6 ;"
            ::
            );**/

            for (int l = 0; l < lwcount; l+= 4) {

                lwin128 = _mm_load_ps(linwindow + l); //loads 4 floats from the base pointer
                lwing128 = _mm_load_ps(linwindowg + l);

                __asm__(
                    "movaps %1, %%xmm0 ;"
                    "movaps %2, %%xmm1 ;"
                    "movaps %3, %%xmm2 ;"
                    "subps %%xmm0, %%xmm1 ;" //subtract window from the center pixel (format :: dest = dest - src)
                    "mulps %%xmm1, %%xmm1 ;" //square the difference (photometric dist)
                    "divps %%xmm2, %%xmm1 ;" //divide with -2*sigma_sq
                    "movaps %%xmm1, %0 ;"
                    :"=m"(pmgc128)
                    :"m"(lwin128), "m"(pimg128), "m"(sigmaphs128)
                    : "%xmm1" //clobbered register
                );

                //calculate the exponential using SVML intrinsic (no assembly available)
                pmgc128  = _mm_exp_ps(pmgc128); //photometric gaussian coefficient

                __asm__(
                    "movaps %2, %%xmm0 ;"
                    "movaps %3, %%xmm1 ;"
                    "movaps %4, %%xmm2 ;"
                    "mulps %%xmm1, %%xmm2 ;" //bl filter coefficient
                    "movaps %%xmm2, %0 ;"    //store the filter coefficient
                    "mulps %%xmm0, %%xmm2 ;" //filtered coefficient
                    "movaps %%xmm2, %1;"     //store the filtered coefficient
                    :"=m"(norm128), "=m"(blfiltered128)
                    :"m"(lwin128), "m"(lwing128), "m"(pmgc128)
                    :"%xmm2" //clobbered
                );

                __asm__(
                    "movaps %2, %%xmm3 ;"  //TODO: try using different registers in these 3 fragments of assembly
                    "movaps %3, %%xmm4 ;"
                    "movaps %4, %%xmm5 ;"
                    "movaps %5, %%xmm6 ;"
                    "haddps %%xmm3, %%xmm5 ;"
                    "haddps %%xmm4, %%xmm6 ;" //corresponding intrinsic --> _mm_hadd_ps
                    "movaps %%xmm5, %0 ;"
                    "movaps %%xmm6, %1 ;"
                    :"=m"(o1), "=m"(o2)
                    :"m"(norm128), "m"(blfiltered128), "m"(o1), "m"(o2)
                    :"%xmm5", "%xmm6"
                );

                /**
                if (i == 200 && j == 200) {
                printf("(%f, %f, %f, %f)\n", o1.x, o1.y, o1.z, o1.w);
                printf("(%f, %f, %f, %f)\n", o2.x, o2.y, o2.z, o2.w);
                 }**/

            } //inner loop


            //TODO: talk how this loop is transformed into vector intrinsic in presentation (a form of reduction)//
            /**
            for(int q = 0; q < lwcount ; q++) {
            	filtered_pixel +=  blfiltered[q];
            	normal_factor += normwin[q];
            }**/

            filtered_pixel = o2.x;
            normal_factor = o1.x;
            filtered_pixel = filtered_pixel / normal_factor;

            //printf("(%d, %d) = %f\t%f\n", (i - filter_hw), (j - filter_hw), normal_factor, filtered_pixel);
            filteredImg(i - filter_hw, j - filter_hw) =  filtered_pixel;

        }

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








