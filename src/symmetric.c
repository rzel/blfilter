#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>


#define HALF_WIDTH 10
#define SIGMA 25.0
#define SIGSP 50.0
#define FULLWIDTH (2*HALF_WIDTH+1)


struct uchar4{
	unsigned char x, y, z, w;
};


struct inputVariables {
	 int 
		jMax, 	//hieght of the matrix
		iMax,	//Width of the matrix
		size,
		halfWidth,
		width;
};
struct hostVariables {
	struct uchar4 pix;
	unsigned char	* pixelArray;
	float * gaussian_sp_WINDOW_ARRAY;
};
struct deviceVariables {
	unsigned char	* pixelArray;
	float
			* Array1,
			* Array2;
	float * gaussian_sp_WINDOW_ARRAY;
};
//"struct timer" and functions "start," "stop," and "startNewTimer" are for easily setting up timers.

#pragma pack(1)
struct BMPHeader{
    short type;
    int size;
    short reserved1;
    short reserved2;
    int offset;
};
#pragma pack(1)
struct BMPInfoHeader{
    unsigned int size;
    int width;
    int height;
    unsigned short planes;
    unsigned short int bitsPerPixel;
    unsigned int compression;
    unsigned int imageSize;
    int xPelsPerMeter;
    int yPelsPerMeter;
    unsigned int clrUsed;
    unsigned int clrImportant;
};
void * memset ( void * ptr, int value, size_t num );
unsigned char *LoadBMPFile(BMPHeader *hdr, BMPInfoHeader *infoHdr, struct uchar4 **dst, int *width, int *height, const char *name){

    int x, y;
    FILE *fd;
    printf("Loading %s...\n", name);
    if( !(fd = fopen(name,"rb")) ){
        printf("***BMP load error: file access denied***\n");
        exit(0);
    }

    fread(hdr, sizeof(*hdr), 1, fd);
    if((*hdr).type != 0x4D42){
        printf("***BMP load error: bad file format***\n");
        exit(0);
    }
    fread(infoHdr, sizeof(*infoHdr), 1, fd);

    if((*infoHdr).compression){
        printf("***BMP load error: compressed image***\n");
        exit(0);
    }

    *width  = abs((*infoHdr).width);
    *height = abs((*infoHdr).height);
	int malsize = (*width+HALF_WIDTH*2) * (*height+HALF_WIDTH) * 4;
    *dst    = (uchar4 *)malloc(malsize);

    printf("BMP width: %d\n", *width);
    printf("BMP height: %d\n", *height);
	*width += HALF_WIDTH*2;
	*height += HALF_WIDTH;
	
	int jMax = *width;
	int iMax = *height;
	
   for(y = 0; y < *height-HALF_WIDTH; y++){		
		for (int j = 0; j < HALF_WIDTH; j++) {
			(*dst)[y*jMax + j].x = 120;
			(*dst)[y*jMax + j].y = 120;
			(*dst)[y*jMax + j].z = 120;
			(*dst)[y*jMax + j].w = 120;
		}
        for(x = HALF_WIDTH; x < jMax-HALF_WIDTH; x++){
			(*dst)[( y * (*width) + x)].x = fgetc(fd);
            (*dst)[( y * (*width) + x)].y = fgetc(fd);
            (*dst)[( y * (*width) + x)].z = fgetc(fd);
            (*dst)[( y * (*width) + x)].w = fgetc(fd);
        }
		for (int j = jMax-HALF_WIDTH; j < jMax; j++) {
			(*dst)[y*jMax + j].x = 120;
			(*dst)[y*jMax + j].y = 120;
			(*dst)[y*jMax + j].z = 120;
			(*dst)[y*jMax + j].w = 120;
		}
	}
	for ( int i=iMax-HALF_WIDTH; i<iMax; i++) {
		for ( int j=0; j<jMax; j++) {
			(*dst)[i*jMax + j].x = 120;
			(*dst)[i*jMax + j].y = 120;
			(*dst)[i*jMax + j].z = 120;
			(*dst)[i*jMax + j].w = 120;
		}
	}	


    fclose(fd);

	fd=fopen("new.bmp","wb");//makes new .bmp file for bin writing
	
	(*infoHdr).width += (((*infoHdr).width < 0) ? - 1 : 1) * HALF_WIDTH*2;
	(*infoHdr).height += (((*infoHdr).height < 0) ? - 1 : 1) * HALF_WIDTH;
	//writes in the headers I got from the other file
	fwrite(hdr,sizeof(*hdr),1,fd);
	fwrite(infoHdr,sizeof(*infoHdr),1,fd);

	//writes in the image information from the other image
	fwrite(*dst,sizeof(uchar4 *),malsize,fd);
	fclose(fd); //closes file

	struct uchar4 * pixelArray = (*dst);
	
	struct uchar4 * pixelArray2 = (uchar4*)malloc(jMax * iMax * 4);
	unsigned char * pix = (unsigned char*)malloc(jMax * iMax * sizeof(unsigned char));
	
	
	
	
	for (int i = 0; i < iMax; i++)
		for (int j = 0; j < jMax; j++) {
			pixelArray2[i*jMax+j].x = (unsigned char)(0.299 * (float)pixelArray[i*jMax+j].x + 0.587 * (float)pixelArray[i*jMax+j].y + 0.114 * (float)pixelArray[i*jMax+j].z);
			pixelArray2[i*jMax+j].y = (unsigned char)(0.299 * (float)pixelArray[i*jMax+j].x + 0.587 * (float)pixelArray[i*jMax+j].y + 0.114 * (float)pixelArray[i*jMax+j].z);
			pixelArray2[i*jMax+j].z = (unsigned char)(0.299 * (float)pixelArray[i*jMax+j].x + 0.587 * (float)pixelArray[i*jMax+j].y + 0.114 * (float)pixelArray[i*jMax+j].z);
			pixelArray2[i*jMax+j].w = (unsigned char)(0.299 * (float)pixelArray[i*jMax+j].x + 0.587 * (float)pixelArray[i*jMax+j].y + 0.114 * (float)pixelArray[i*jMax+j].z);
			pix[i*jMax+j] = pixelArray2[i*jMax+j].x;
		}

	fd=fopen("new2.bmp","wb");//makes new .bmp file for bin writing
	fwrite(hdr,sizeof(*hdr),1,fd);
	fwrite(infoHdr,sizeof(*infoHdr),1,fd);	
	fwrite(pixelArray2,sizeof(uchar4 *),malsize,fd);
	fclose(fd); //closes file

	return pix;
};

inline float fsquare(float n) {
	return n*n;
};
int main () {

	printf("program started\n");

	inputVariables inp;
	hostVariables hst;
	inp.halfWidth = HALF_WIDTH;
	inp.width = inp.halfWidth * 2;

    BMPHeader hdr;
    BMPInfoHeader infoHdr;
	hst.pixelArray = LoadBMPFile(&hdr, &infoHdr, (uchar4 **)&hst.pix, &inp.jMax, &inp.iMax, "./me.bmp");
	printf("Finished making array on host\n");

struct timeval start_time;

struct timeval end_time;

gettimeofday(&start_time,NULL);
		
//	timer stopWatch = startNewTimer();
	
	inp.size = inp.jMax*inp.iMax;

    //gaussian kernel for the spatial distances
    int windowSize = FULLWIDTH*(HALF_WIDTH+1)*sizeof(float);
	hst.gaussian_sp_WINDOW_ARRAY = (float*)malloc(windowSize);
    for (int i = 0; i<= HALF_WIDTH ; i++) {
        for (int j = -HALF_WIDTH; j<= HALF_WIDTH; j++) {
            hst.gaussian_sp_WINDOW_ARRAY[i*FULLWIDTH+j+HALF_WIDTH] = exp( -(j*j + i*i) / (2 * SIGSP * SIGSP) );
        }
    }
	float sigmaDerived = 1/(-2*SIGMA*SIGMA);

	float * Array1 = (float*)malloc(inp.size*sizeof(float));
	float * Array2 = (float*)malloc(inp.size*sizeof(float));	
	memset(Array1, 0, inp.size*sizeof(float));
	memset(Array2, 0, inp.size*sizeof(float));	

	int inx1, inx2;
	float reg1, reg2;
	float gaussian_bl, gaussian_ph;
	for (int l = 0; l < inp.iMax-HALF_WIDTH; l++) {
		for (int k = HALF_WIDTH; k < inp.jMax-HALF_WIDTH; k++) {
			inx1 = l*inp.jMax+k;
			reg1 = 0.0f;
			reg2 = 0.0f;
			for (int i = 1; i <= HALF_WIDTH; i++) {
				for (int j = -HALF_WIDTH; j <= HALF_WIDTH; j++) {
					inx2 = inx1 + i*inp.jMax + j;
					gaussian_ph = expf(sigmaDerived * fsquare(hst.pixelArray[inx1]-hst.pixelArray[inx2]));
					gaussian_bl = gaussian_ph * hst.gaussian_sp_WINDOW_ARRAY[i*FULLWIDTH+j+HALF_WIDTH];
					//add results to center pixel
					reg1 += gaussian_bl * hst.pixelArray[inx2];
					reg2 += gaussian_bl;
					Array1[inx2] += gaussian_bl * hst.pixelArray[inx1];
					Array2[inx2] += gaussian_bl;					
				}
			}
			for (int j = -1; j >= -HALF_WIDTH; j--) {
				inx2 = inx1 + j;
				gaussian_ph = expf(sigmaDerived * fsquare(hst.pixelArray[inx1]-hst.pixelArray[inx2]));
				gaussian_bl = gaussian_ph * hst.gaussian_sp_WINDOW_ARRAY[j+HALF_WIDTH];
				//add results to center pixel
				reg1 += gaussian_bl * hst.pixelArray[inx2];
				reg2 += gaussian_bl;				
				Array1[inx2] += gaussian_bl * hst.pixelArray[inx1];
				Array2[inx2] += gaussian_bl;
			}
			Array1[inx1] += reg1;
			Array2[inx1] += reg2;			
		}
	}

	//copy results to host

//	stop(&stopWatch);


gettimeofday(&end_time,NULL);

double start_dtime=(double)start_time.tv_sec+(double)start_time.tv_usec/1000000.0;

double end_dtime=(double)end_time.tv_sec+(double)end_time.tv_usec/1000000.0; double diff=end_dtime-start_dtime;

printf("%ld %ld %f %ld %ld %f %f\n",start_time.tv_sec,start_time.tv_usec,start_dtime,end_time.tv_sec,end_time.tv_usec,end_dtime,diff);

	//////////////////////////////////////////////////////////////////////////////////////////////	
	//store final results to file
	uchar4 *test = (uchar4*)malloc(inp.size*sizeof(uchar4));
	for (int i = 0; i < inp.iMax; i++)
		for (int j = 0; j < inp.jMax; j++) {
			test[i*inp.jMax+j].x = (unsigned char)Array1[i*inp.jMax+j];
			test[i*inp.jMax+j].y = (unsigned char)Array1[i*inp.jMax+j];
			test[i*inp.jMax+j].z = (unsigned char)Array1[i*inp.jMax+j];
			test[i*inp.jMax+j].w = (unsigned char)Array1[i*inp.jMax+j];
		}
	FILE * fd=fopen("new3.bmp","wb");//makes new .bmp file for bin writing
	fwrite(&hdr,sizeof(hdr),1,fd);
	fwrite(&infoHdr,sizeof(infoHdr),1,fd);	
	fwrite(test,sizeof(uchar4 *),inp.size*sizeof(uchar4),fd);
	fclose(fd); //closes file

	printf("\n");
//	printf("\nElapsed time = %f\n",stopWatch.elapsedTime);
//	cudaFree(hst.pixelArray);
	
	return 0;
}
