
/* BmpProcess.h */



int		ReadInt( FILE * );
short	ReadShort( FILE * );

int     WriteBmp( const char *filename );
void	WriteInt( FILE *, int );
void	WriteShort( FILE *, short int );
unsigned char *
BmpToTexture( char *filename, int *width, int *height );
int
Write_Array_To_Bmp(const char *filename,
				   unsigned char *array,
				   int width, int height);
