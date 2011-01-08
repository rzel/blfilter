
#include "CImg.h"
//#include <omp.h>
#include<pmmintrin.h> //SSE3
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

// this function returns the filtered image
clock_t start, end;
float elapsed;

start = clock();
//call the filtering function
CImg<float> filteredImg = blfilter(grayImg, filter_hw, sigma_sp, sigma_ph);
end = clock();
elapsed = ((float) (end - start)) / CLOCKS_PER_SEC;
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
			paddedImg(i+filter_hw, j+filter_hw) = grayImg(i,j);	

	/**
	CImgDisplay pad_disp(paddedImg, "padded image");
	//create event loop to display the image
	while (!pad_disp.is_closed()) {
		pad_disp.wait();  
		//fil_disp.wait();       
        }
	**/	 	

	float local_window[2*filter_hw + 1][2*filter_hw + 1]; 	 
	
	float local_pixel0 = 0; //local window
	float local_pixel1 = 0;
	float local_pixel2 = 0;
	float local_pixel3 = 0;

	float pd0 = 0; //photometric distance
	float pd1 = 0;
	float pd2 = 0;
	float pd3 = 0;

	float gaussian_ph0 = 0; //photometric filter coefficient
	float gaussian_ph1 = 0;
	float gaussian_ph2 = 0;
	float gaussian_ph3 = 0;

	float gaussian_bl0 = 0; //bilateral filter coefficient
	float gaussian_bl1 = 0;
	float gaussian_bl2 = 0;
	float gaussian_bl3 = 0;

	float filtered_pixel0 = 0; //filtered value of pixel 
	float filtered_pixel1 = 0;
	float filtered_pixel2 = 0;
	float filtered_pixel3 = 0;	

	float filtered_pixel = 0;

	float normal_factor = 0; //total weight of the filter window used to normalize the filtered coefficient


			//linear gaussian window
			int lwgcount = 0;
			float linwindowg[(2*filter_hw+1)*(2*filter_hw+1)];
			for(int k = 0 ; k <= 2*filter_hw+1 ; k++) 
				for(int l = -filter_hw ; l <= 2*filter_hw+1 ; l++) 
					linwindowg[lwgcount++] = gaussian_sp[k][l];

	//#pragma ivdep
	//#pragma omp parallel for
	for(int i = filter_hw; i < dx + filter_hw ; i++)
		for(int j = filter_hw; j < dy + filter_hw; j++) {
				
		filtered_pixel = 0;
		filtered_pixel0 = 0;
		filtered_pixel1 = 0;
		filtered_pixel2 = 0;
		filtered_pixel3 = 0;

		normal_factor = 0;

		//1. calculate the photometric gaussian distances and the filter coefficent (no need to create a full kernel)	
		//2. multiply each pixel of the window with the spatial and photometric coefficients
		int lwcount = 0;
		float linwindow[(2*filter_hw+1)*(2*filter_hw+1)];
		for(int k = -filter_hw ; k <= filter_hw ; k++) 
			for(int m = -filter_hw ; m <= filter_hw ; m++) 
				linwindow[lwcount++] = paddedImg(i + k, j + m);

		//#pragma ivdep
		//for(int k = -filter_hw ; k <= filter_hw ; k++) 
		//	for(int l = -filter_hw ; l <= filter_hw ; l+=4) {			

		for(int l = 0; l < lwcount ; l+= 4) {

			
			//local_pixel0 = paddedImg(i + k , j + l);
			//local_pixel1 = paddedImg(i + k, j + (l + 1));
			//local_pixel2 = paddedImg(i + k, j + (l + 2));
			//local_pixel3 = paddedImg(i + k, j + (l + 3));			
			
			local_pixel0 = linwindow[l];
			local_pixel1 = linwindow[l + 1];
			local_pixel2 = linwindow[l + 2];
			local_pixel3 = linwindow[l + 3];			
			
			//photometric distance = difference in pixel values (squared)
			pd0 = paddedImg(i,j) - local_pixel0;
			pd0 = pd0 * pd0;  
			pd1 = paddedImg(i,j) - local_pixel1;
			pd1 = pd1 * pd1;  
			pd2 = paddedImg(i,j) - local_pixel2;
			pd2 = pd2 * pd2;  
			pd3 = paddedImg(i,j) - local_pixel3;
			pd3 = pd3 * pd3;  

			gaussian_ph0 = exp( -(pd0) / (2*sigma_ph*sigma_ph) );
			gaussian_bl0 = gaussian_ph0 * linwindowg[l] ;
			gaussian_ph1 = exp( -(pd1) / (2*sigma_ph*sigma_ph) );
			gaussian_bl1 = gaussian_ph1 * linwindowg[l+1] ;
			gaussian_ph2 = exp( -(pd2) / (2*sigma_ph*sigma_ph) );
			gaussian_bl2 = gaussian_ph2 * linwindowg[l+2] ;
			gaussian_ph3 = exp( -(pd3) / (2*sigma_ph*sigma_ph) );
			gaussian_bl3 = gaussian_ph3 * linwindowg[l+3] ;

			//if(k == filter_hw){gaussian_bl2 = 0; gaussian_bl3 = 0;}
			//if(l == filter_hw){gaussian_bl1 = 0; gaussian_bl3 = 0;}
			if(l == lwcount) {gaussian_bl1 = 0; gaussian_bl2 = 0; gaussian_bl3 = 0; continue;}
			if(l == lwcount-1) {gaussian_bl2 = 0; gaussian_bl3 = 0; continue;}
			if(l == lwcount-2) { gaussian_bl3 = 0;continue; }

			filtered_pixel0 += gaussian_bl0 * local_pixel0;	
			filtered_pixel1 += gaussian_bl1 * local_pixel1;
			filtered_pixel2 += gaussian_bl2 * local_pixel2;	
			filtered_pixel3 += gaussian_bl3 * local_pixel3;

			normal_factor += (gaussian_bl0 + gaussian_bl1 + gaussian_bl1 + gaussian_bl3);	
							
		}

		//printf("(%f, %f, %f, %f)\n", filtered_pixel0, filtered_pixel1, filtered_pixel2, filtered_pixel3);
			
		filtered_pixel = (filtered_pixel0 + filtered_pixel1 + filtered_pixel1 + filtered_pixel3) / normal_factor;

		//printf("(%d, %d) = %f\t%f\n", i, j, normal_factor, filtered_pixel);

		filteredImg(i-filter_hw, j - filter_hw) =  filtered_pixel;

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








