
#include "CImg.h"
#include <omp.h>
#include<pmmintrin.h> //SSE3
#include<xmmintrin.h>
#include<ia32intrin.h> //general intrinsics
//#include<smmintrin.h> //extract


using namespace cimg_library;
double my_difftime ();
CImg<float>  blfilter(CImg<float>, int, float, float*);
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


    double start, end;
    double elapsed;

//call the filtering function

    float gaussian_sp[2*filter_hw + 1][2*filter_hw + 1];
    for(int i = -filter_hw; i<= filter_hw ; i++)
        for(int j = -filter_hw; j<= filter_hw; j++) {
            gaussian_sp[i+filter_hw][j+filter_hw] = exp( -(j*j + i*i) / (2*sigma_sp*sigma_sp) );
        }
    start = my_difftime();
    CImg<float> filteredImg = blfilter(grayImg, filter_hw, sigma_ph, gaussian_sp);
    end = my_difftime();
    elapsed = end - start;

    printf("Time taken for filtering = %6.2f \n", elapsed);

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


CImg<float>  blfilter(CImg<float> grayImg, int filter_hw, float sigma_ph, float *gaussian_sp)
{
    //image height and width
    int dx = grayImg.width();
    int dy = grayImg.height();



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

#ifdef USE_OMP
    #pragma omp parallel shared(filteredImg, paddedImg)
    {
    #pragma omp for private(i, j, k, l, 
#endif
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


            struct Vector4 {
                float x, y, z, w;
            };
            struct Vector4 o1, o2; //to store filtered and normal factors

            __m128 norm128;
            __m128 blfiltered128;

            __m128 pmgc128; //to hold the photometric gaussian coefficient
            __m128 lwin128; //for local window
            __m128 lwing128; //for spatial gaussian

            __m128 sigmaphs128 = _mm_set1_ps(sigmasq);
            __m128 pimg128 = _mm_set1_ps(paddedImg(i,j));


            float normwin[sizewin];

            //clear the outputs o1 and o2
            o1.x = o1.y = o1.z = o1.w = 0;
            o2.x = o2.y = o2.z = o2.w = 0;


            for (int l = 0; l < lwcount; l+= 4) {

                lwin128 = _mm_load_ps(linwindow + l); //loads 4 floats from the base pointer
                lwing128 = _mm_load_ps(linwindowg + l);

                __asm__(
                    "movups %1, %%xmm0 ;"
                    "movups %2, %%xmm1 ;"
                    "movups %3, %%xmm2 ;"
                    "subps %%xmm0, %%xmm1 ;" //subtract window from the center pixel (format :: dest = dest - src)
                    "mulps %%xmm1, %%xmm1 ;" //square the difference (photometric dist)
                    "divps %%xmm2, %%xmm1 ;" //divide with -2*sigma_sq
                    "movups %%xmm1, %0 ;"
                    :"=m"(pmgc128)
                    :"m"(lwin128), "m"(pimg128), "m"(sigmaphs128)
                    : "%xmm1" //clobbered register
                );

                //calculate the exponential using SVML intrinsic (no assembly available)
                pmgc128  = _mm_exp_ps(pmgc128); //photometric gaussian coefficient

                __asm__(
                    "movups %2, %%xmm0 ;"
                    "movups %3, %%xmm1 ;"
                    "movups %4, %%xmm2 ;"
                    "mulps %%xmm1, %%xmm2 ;" //bl filter coefficient
                    "movups %%xmm2, %0 ;"    //store the filter coefficient
                    "mulps %%xmm0, %%xmm2 ;" //filtered coefficient
                    "movups %%xmm2, %1;"     //store the filtered coefficient
                    :"=m"(o1), "=m"(o2)
                    :"m"(lwin128), "m"(lwing128), "m"(pmgc128)
                    :"%xmm2" //clobbered
                );


                blfiltered[l + 0] = o2.x;
                blfiltered[l + 1] = o2.y;
                blfiltered[l + 2] = o2.z;
                blfiltered[l + 3] = o2.w;

                normwin[l + 0] = o1.x;
                normwin[l + 1] = o1.y;
                normwin[l + 2] = o1.z;
                normwin[l + 3] = o1.w;


            } //inner loop

            for(int q = 0; q < lwcount ; q++) {
                filtered_pixel +=  blfiltered[q];
                normal_factor += normwin[q];
            }

            filtered_pixel = filtered_pixel / normal_factor;
            filteredImg(i - filter_hw, j - filter_hw) =  filtered_pixel;
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


double my_difftime ()
{
    struct timeval tp;
    struct timezone tzp;
    int i;

    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
};






