#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

int bigEndian()
{
    union {
        uint32_t i;
        char c[4];
    } a = { 0x01020304 };

    return a.c[0] == 1; 
}

void reverseByteOrder(char *dst, size_t nbytes, char *src)
{
    for(int i=0; i<nbytes; ++i) dst[i] = src[nbytes-1-i];
}

#define EXTRACT(type, name, offset)\
    if(!bigEndian()) reverseByteOrder(bytes, sizeof(type), fdata+offset);\
    type name = *(type *)bytes;
    
char bytes[4];  
int TTF2Texture(char *dst, const char *filename)
{
    FILE *f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    size_t fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *fdata = (char*)malloc(fsize+1);
    fread(fdata, 1, fsize, f);
    fclose(f);
    
    // Check if spline outlines are saved in ttf file.
    EXTRACT(uint32_t, sfntVersion, 0);
    if(sfntVersion != 0x00010000) return -1;
    
    EXTRACT(uint16_t, numTables, 4);
    for(int i=0; i<(int)numTables; ++i)
    {
        EXTRACT(uint32_t, tableTagInt, 12 + 16 * i);
        char tableTag[5];
        for(int j=0; j<4; ++j) tableTag[j] = bytes[4-j];
        tableTag[4] = 0;
        
        EXTRACT(uint32_t, checkSum, 12 + 16 * i + 4);
        EXTRACT(uint32_t, offset, 12 + 16 * i + 8);
        EXTRACT(uint32_t, length, 12 + 16 * i + 12);
        
        printf("Table %d: %s with checksum %d, offset %d and length %d\n", i, tableTag, checkSum, offset, length);
    }
    
    return 0;
}

int main(int argc, char **args)
{
    char *texture;
    TTF2Texture(texture, "C:\\Windows\\Fonts\\arial.ttf");
    
    return 0;
}
