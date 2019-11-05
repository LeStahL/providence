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

int bit(char data, int index)
{
    return (data >> index) & 0x1;
}

#define EXTRACT(type, name, offset)\
    if(!bigEndian()) reverseByteOrder(bytes, sizeof(type), fdata+offset);\
    type name = *(type *)bytes;\
    
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
    
    EXTRACT(uint32_t, sfntVersion, 0);
    if(sfntVersion != 0x00010000) return -1;
    
    uint32_t maxpOffset, glyfOffset, headOffset, locaOffset;
    
    EXTRACT(uint16_t, numTables, 4);
    for(int i=0; i<(int)numTables; ++i)
    {
        EXTRACT(uint32_t, tableTagInt, 12 + 16 * i);        
        EXTRACT(uint32_t, checkSum, 12 + 16 * i + 4);
        EXTRACT(uint32_t, offset, 12 + 16 * i + 8);
        EXTRACT(uint32_t, length, 12 + 16 * i + 12);

        if(tableTagInt == 0x6d617870) maxpOffset = offset;
        else if(tableTagInt == 0x676c7966) glyfOffset = offset;
        else if(tableTagInt == 0x68656164) headOffset = offset;
        else if(tableTagInt == 0x6c6f6361) locaOffset = offset;
    }
    
    EXTRACT(uint16_t, numGlyphs, maxpOffset + 4);
    printf("Number of glyphs saved in font: %d\n", (int)numGlyphs);
    
    EXTRACT(int16_t, indexToLocFormat, headOffset + 50);
    if(indexToLocFormat) // Only parse uint32_t based tables
    {
        for(int i=0; i<128; ++i) 
        {
            EXTRACT(uint32_t, glyphOffset, locaOffset + 4 * i);
            printf(">> Glyph offset of %d (%c) is %d\n", (int)i, (char)i, (int)glyphOffset);
            
            EXTRACT(int16_t, numberOfContours, glyfOffset + glyphOffset);
            printf("Number of contours in glyph: %d\n", (int)numberOfContours);
            
            //BEGIN FIXME: remove, unnecessary
            EXTRACT(int16_t, xMin, glyfOffset + glyphOffset + 2);
            EXTRACT(int16_t, yMin, glyfOffset + glyphOffset + 4);
            EXTRACT(int16_t, xMax, glyfOffset + glyphOffset + 6);
            EXTRACT(int16_t, yMax, glyfOffset + glyphOffset + 8);
            printf("Bounding box: (%d,%d), (%d,%d)\n", xMin, yMin, xMax, yMax);
            //END FIXME
            
            printf("---\n");
            if(glyfOffset >= 0)
            {
                for(int j=0; j<numberOfContours; ++j)
                {
                    EXTRACT(uint16_t, endPointOfContours, glyfOffset + glyphOffset + 10 + 2 * j);
                    printf("end point: %d\n", endPointOfContours);
                }
                
                EXTRACT(uint16_t, instructionLength, glyfOffset + glyphOffset + 10 + 2*numberOfContours);
                printf("Number of instrutions: %d\n", (int)instructionLength);
                
                int nPointsX = 0, nPointsY = 0;
                for(int j=0; j<instructionLength; ++j)
                {
                    EXTRACT(uint8_t, instruction, glyfOffset + glyphOffset + 10 + 2*numberOfContours + j);
                    
                    int onCurve = bit(instruction, 0);
                    size_t xWidth = bit(instruction, 1) + 1, 
                        yWidth = bit(instruction, 2) + 1;
                        
                    
                }
            }
            printf("\n");
        }
    }
    
    
    return 0;
}

int main(int argc, char **args)
{
    char *texture;
//     TTF2Texture(texture, "C:\\Windows\\Fonts\\arial.ttf");
    TTF2Texture(texture, "C:\\Users\\NR4\\Downloads\\Princess_Sofia\\PrincessSofia-Regular.ttf");
    
//     for(int i=0; i<256; ++i)
//     {
//         for(int j=0; j<8; ++j) printf("%d", bit(i,j));
//         printf("\n");
//     }
    
    
    return 0;
}
