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
    for(int i=0; i<nbytes; ++i)
        dst[i] = src[nbytes-1-i];
}

#define ACQUIRE(type, name)\
    fread(bytes1, sizeof(char), sizeof(type), f);\
    reverseByteOrder(bytes2, sizeof(type), bytes1);\
    type name = *(type *)bytes2;

int TTF2Texture(char *dst, const char *filename)
{
    FILE *f = fopen(filename, "rb");
    char bytes1[4], bytes2[4];
    
    if(!bigEndian())
    {
        ACQUIRE(uint32_t, sfntVersion);
        if(sfntVersion != 0x00010000) return -1;
        
        ACQUIRE(uint16_t, numTables);
        
        fseek(f, 12, SEEK_SET);
        for(int i=0; i<(int)numTables; ++i)
        {
            ACQUIRE(uint32_t, tableTagInt);
            char tableTag[5];
            for(int j=0; j<4; ++j) tableTag[j] = bytes2[j];
            tableTag[4] = 0;
            
            ACQUIRE(uint32_t, checkSum);
            ACQUIRE(uint32_t, offset);
            ACQUIRE(uint32_t, length);
            
            printf("Table %d: %s with checksum %d, offset %d and length %d\n", i, tableTag, checkSum, offset, length);
        }
    }
    else 
    {
    }
    
    return 0;
}

int main(int argc, char **args)
{
    char *texture;
    TTF2Texture(texture, "C:\\Windows\\Fonts\\arial.ttf");
    
    return 0;
}
