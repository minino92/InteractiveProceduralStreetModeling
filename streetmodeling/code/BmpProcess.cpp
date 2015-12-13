#include "stdafx.h"

#include "BmpProcess.h"

#include <GL/glut.h>

const int BIRGB = { 0 };
unsigned char _pad[] = { 0, 0, 0 };

struct BMPFILEHEADER
{
	short bfType;
	int bfSize;
	short bfReserved1;
	short bfReserved2;
	int bfOffBits;
} FileHeader;

struct BMPINFOHEADER
{
	int biSize;
	int biWidth;
	int biHeight;
	short biPlanes;
	short biBitCount;
	int biCompression;
	int biSizeImage;
	int biXPelsPerMeter;
	int biYPelsPerMeter;
	int biClrUsed;
	int biClrImportant;
} InfoHeader;

void
WriteInt( FILE *fp, int i )
{
	fputc( (i>> 0) & 0xff, fp );
	fputc( (i>> 8) & 0xff, fp );
	fputc( (i>>16) & 0xff, fp );
	fputc( (i>>24) & 0xff, fp );
}


void
WriteShort( FILE *fp, short int i )
{
	fputc( (i>>0) & 0xff, fp );
	fputc( (i>>8) & 0xff, fp );
}


int
WriteBmp( const char *filename )
{
	FILE *fp;			/* pointer to BMP file		*/
	unsigned char *array;		/* array to hold the RGB dump	*/
	unsigned char *p;		/* pointer into array[]		*/
	int dx, dy;			/* window size in pixels	*/
	int i, j;			/* counters			*/
	int numpad;			/* # if bytes to pad with	*/


	/* try to open the file:					*/

	fp = fopen( filename, "wb" );
	if( fp == NULL )
	{
		fprintf( stderr, "Cannot open BMP file '%s'\n", filename );
		return -1;
	}


	/* allocate the array to hold the RGB pixel values:		*/

	int viewport[4] = {0};

	glGetIntegerv(GL_VIEWPORT, viewport);

	//glutSetWindow( MainWindow );
	//dx = glutGet( GLUT_WINDOW_WIDTH );
	//dy = glutGet( GLUT_WINDOW_HEIGHT );
	dx = viewport[2] - viewport[0];
	dy = viewport[3] - viewport[1];

	array = new unsigned char[3*dx*dy];
	if( array == NULL )
	{
		fprintf( stderr, "WriteBmp() cannot allocate the temporary array!\n" );
		fclose( fp );
		unlink( filename );	/* remove the file we just created */
		return -1;
	}


	/* setup and do the dump:					*/

	glPixelStorei( GL_PACK_ALIGNMENT, 4 );
	//glPixelStorei( GL_PACK_ROW_LENGTH, dx );
	//glReadBuffer( GL_FRONT );
	glReadBuffer( GL_BACK );
	glReadPixels( 0, 0, dx, dy, GL_RGB, GL_UNSIGNED_BYTE, array );



	// First write out the bitmap file header:

	BITMAPFILEHEADER fh;
	fh.bfType = 0x4d42;

	int dwidth = dx;
	int dheight= dy;
	int destsize = 3 * dwidth * dheight;


	// Write the file header

	fh.bfOffBits = 14 + 40;
	fh.bfSize = fh.bfOffBits + destsize;
	fh.bfReserved1 = 0;
	fh.bfReserved2 = 0;

	WriteShort( fp, fh.bfType );
	WriteInt( fp, fh.bfSize );
	WriteShort( fp, fh.bfReserved1 );
	WriteShort( fp, fh.bfReserved2 );
	WriteInt( fp, fh.bfOffBits );


	// Write the BITMAPINFOHEADER

	BITMAPINFOHEADER bi;
	bi.biSize = 40;
	bi.biWidth = dwidth;
	bi.biHeight = dheight;
	bi.biPlanes = 1;
	bi.biBitCount = 24;
	bi.biCompression = BI_RGB;
	bi.biSizeImage = 0;
	//bi.biXPelsPerMeter = 100;
	//bi.biYPelsPerMeter = 100;
	bi.biXPelsPerMeter = 4; //use lower resolution here
	bi.biYPelsPerMeter = 4;
	bi.biClrUsed = 0;
	bi.biClrImportant = 0;

	WriteInt( fp, bi.biSize );
	WriteInt( fp, bi.biWidth );
	WriteInt( fp, bi.biHeight );
	WriteShort( fp, bi.biPlanes );
	WriteShort( fp, bi.biBitCount );
	WriteInt( fp, bi.biCompression );
	WriteInt( fp, bi.biSizeImage );
	WriteInt( fp, bi.biXPelsPerMeter );
	WriteInt( fp, bi.biYPelsPerMeter );
	WriteInt( fp, bi.biClrUsed );
	WriteInt( fp, bi.biClrImportant );

	numpad = 4*(int)((dwidth*3+3)/4) -  3*dwidth;
	




	/* write the pixels:						*/

	for( i = 0; i < dy; i++ )
	{
		p = &array[ 3*dx*i ];
		for( j = 0; j < dx; j++, p += 3 )
		{
			fputc( p[2], fp );
			fputc( p[1], fp );
			fputc( p[0], fp );
			fwrite( _pad, 1, numpad, fp );
		}
	}


	/* end game:							*/

	fflush( fp );
	fclose( fp );
	delete [] array;
	return 0;
}


int
Write_Array_To_Bmp(const char *filename,
				   unsigned char *array,
				   int width, int height)
{
	FILE *fp;			/* pointer to BMP file		*/
	int i, j;			/* counters			*/
	
	/* try to open the file:					*/
	int dx, dy;			/* window size in pixels	*/
	dx = width;
	dy = height;
	unsigned char *p;		/* pointer into array[]		*/
	int numpad;			/* # if bytes to pad with	*/

	fp = fopen( filename, "wb" );
	if( fp == NULL )
	{
		fprintf( stderr, "Cannot open BMP file '%s'\n", filename );
		return -1;
	}

	// First write out the bitmap file header:

	BITMAPFILEHEADER fh;
	fh.bfType = 0x4d42;

	int dwidth = dx;
	int dheight= dy;
	int destsize = 3 * dwidth * dheight;

	// Write the file header

	fh.bfOffBits = 14 + 40;
	fh.bfSize = fh.bfOffBits + destsize;
	fh.bfReserved1 = 0;
	fh.bfReserved2 = 0;

	WriteShort( fp, fh.bfType );
	WriteInt( fp, fh.bfSize );
	WriteShort( fp, fh.bfReserved1 );
	WriteShort( fp, fh.bfReserved2 );
	WriteInt( fp, fh.bfOffBits );


	// Write the BITMAPINFOHEADER

	BITMAPINFOHEADER bi;
	bi.biSize = 40;
	bi.biWidth = dwidth;
	bi.biHeight = dheight;
	bi.biPlanes = 1;
	bi.biBitCount = 24;
	bi.biCompression = BI_RGB;
	bi.biSizeImage = 0;
	//bi.biXPelsPerMeter = 100;
	//bi.biYPelsPerMeter = 100;
	bi.biXPelsPerMeter = 4; //use lower resolution here
	bi.biYPelsPerMeter = 4;
	bi.biClrUsed = 0;
	bi.biClrImportant = 0;

	WriteInt( fp, bi.biSize );
	WriteInt( fp, bi.biWidth );
	WriteInt( fp, bi.biHeight );
	WriteShort( fp, bi.biPlanes );
	WriteShort( fp, bi.biBitCount );
	WriteInt( fp, bi.biCompression );
	WriteInt( fp, bi.biSizeImage );
	WriteInt( fp, bi.biXPelsPerMeter );
	WriteInt( fp, bi.biYPelsPerMeter );
	WriteInt( fp, bi.biClrUsed );
	WriteInt( fp, bi.biClrImportant );

	numpad = 4*(int)((dwidth*3+3)/4) -  3*dwidth;
	




	/* write the pixels:						*/

	for( i = 0; i < dy; i++ )
	{
		p = &array[ 3*dx*i ];
		for( j = 0; j < dx; j++, p += 3 )
		{
			fputc( p[2], fp );
			fputc( p[1], fp );
			fputc( p[0], fp );
			fwrite( _pad, 1, numpad, fp );
		}
	}


	/* end game:							*/

	fflush( fp );
	fclose( fp );
	delete [] array;
	return 0;
}



/**
 ** read a BMP file into a Texture:
 **/

unsigned char *
BmpToTexture( char *filename, int *width, int *height )
{

	int s, t, e;		// counters
	int numextra;		// # extra bytes each line in the file is padded with
	FILE *fp;
	unsigned char *texture;
	int nums, numt;
	unsigned char *tp;


	fp = fopen( filename, "rb" );
	if( fp == NULL )
	{
		fprintf( stderr, "Cannot open Bmp file '%s'\n", filename );
		return NULL;
	}

	FileHeader.bfType = ReadShort( fp );


	// if bfType is not 0x4d42, the file is not a bmp:

	if( FileHeader.bfType != 0x4d42 )
	{
		fprintf( stderr, "Wrong type of file: 0x%0x\n", FileHeader.bfType );
		fclose( fp );
		return NULL;
	}


	FileHeader.bfSize = ReadInt( fp );
	FileHeader.bfReserved1 = ReadShort( fp );
	FileHeader.bfReserved2 = ReadShort( fp );
	FileHeader.bfOffBits = ReadInt( fp );


	InfoHeader.biSize = ReadInt( fp );
	InfoHeader.biWidth = ReadInt( fp );
	InfoHeader.biHeight = ReadInt( fp );

	nums = InfoHeader.biWidth;
	numt = InfoHeader.biHeight;

	InfoHeader.biPlanes = ReadShort( fp );
	InfoHeader.biBitCount = ReadShort( fp );
	InfoHeader.biCompression = ReadInt( fp );
	InfoHeader.biSizeImage = ReadInt( fp );
	InfoHeader.biXPelsPerMeter = ReadInt( fp );
	InfoHeader.biYPelsPerMeter = ReadInt( fp );
	InfoHeader.biClrUsed = ReadInt( fp );
	InfoHeader.biClrImportant = ReadInt( fp );


	// fprintf( stderr, "Image size found: %d x %d\n", ImageWidth, ImageHeight );


	texture = new unsigned char[ 3 * nums * numt ];
	if( texture == NULL )
	{
		fprintf( stderr, "Cannot allocate the texture array!\b" );
		return NULL;
	}


	// extra padding bytes:

	numextra =  4*(( (3*InfoHeader.biWidth)+3)/4) - 3*InfoHeader.biWidth;


	// we do not support compression:

	if( InfoHeader.biCompression != BI_RGB )
	{
		fprintf( stderr, "Wrong type of image compression: %d\n", InfoHeader.biCompression );
		fclose( fp );
		return NULL;
	}
	


	rewind( fp );
	fseek( fp, 14+40, SEEK_SET );

	if( InfoHeader.biBitCount == 24 )
	{
		for( t = 0, tp = texture; t < numt; t++ )
		{
			for( s = 0; s < nums; s++, tp += 3 )
			{
				*(tp+2) = fgetc( fp );		// b
				*(tp+1) = fgetc( fp );		// g
				*(tp+0) = fgetc( fp );		// r
			}

			for( e = 0; e < numextra; e++ )
			{
				fgetc( fp );
			}
		}
	}

	fclose( fp );

	*width = nums;
	*height = numt;
	return texture;
}



int
ReadInt( FILE *fp )
{
	unsigned char b3, b2, b1, b0;
	b0 = fgetc( fp );
	b1 = fgetc( fp );
	b2 = fgetc( fp );
	b3 = fgetc( fp );
	return ( b3 << 24 )  |  ( b2 << 16 )  |  ( b1 << 8 )  |  b0;
}


short
ReadShort( FILE *fp )
{
	unsigned char b1, b0;
	b0 = fgetc( fp );
	b1 = fgetc( fp );
	return ( b1 << 8 )  |  b0;
}
