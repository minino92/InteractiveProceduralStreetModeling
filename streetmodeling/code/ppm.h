/*

This file contains routines to read and write a .ppm file

Eugene Zhang  January 2005

*/

#ifndef ppm_is_defined
#define ppm_is_defined

int SSS_SavePPM(int w, int h,unsigned char* data, char* fName){
        FILE* f_p = fopen(fName, "wb");
        fprintf(f_p, "P6\n%i %i\n255\n",w,h);

        for(int i = h-1; i >=0; i--){
                fwrite( &(data[i*w*3]), sizeof( char ), w*3, f_p);
                fflush(0);
        }
        fclose(f_p);
        return 0;
}

int SSS_ReadPPM(int *w, int *h,unsigned char* data, char* fName){
	int dyn_range;  // usually 255, so no use
	unsigned char temp;

        FILE* f_p = fopen(fName, "rb");
        fscanf(f_p, "P6\n%i %i\n%i",w,h, &dyn_range);
				fread(&temp, sizeof (unsigned char), 1, f_p);

        for(int i = *h-1; i >= 0; i--){
                fread( &(data[i*(*w)*3]), sizeof( char ), *w*3, f_p);
                fflush(0);
        }
        fclose(f_p);

        return 0;
}
#endif