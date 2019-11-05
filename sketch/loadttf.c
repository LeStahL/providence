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

void reverseByteOrder(char *dst, int nbytes, char *src)
{
    for(int i=0; i<nbytes; ++i) dst[i] = src[nbytes-1-i];
}

uint8_t bit(uint32_t data, size_t index)
{
    return (data >> index) & 0x1;
}

// BEGIN TODO debug only
void printBits(uint8_t data, size_t len)
{
    for(int j=len-1; j>=0; j-=1) printf("%d", bit((uint8_t)data,j));
}
// END TODO debug only

#define EXTRACT_REVERSED(type, name, offset)\
    if(!bigEndian()) reverseByteOrder(bytes, sizeof(type), fdata+offset);\
    type name = *(type *)bytes;
    
#define EXTRACT(type, name, offset)\
    type name = *(type *)bytes;
    
char bytes[4];  
int TTF2Texture(uint16_t *x, uint16_t *y, uint8_t *flags, const char *filename)
{
    FILE *f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    size_t fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *fdata = (char*)malloc(fsize+1);
    fread(fdata, 1, fsize, f);
    fclose(f);
    
//     EXTRACT_REVERSED(uint32_t, sfntVersion, 0);
//     if(sfntVersion != 0x00010000) return -1;
    
    uint32_t maxpOffset, glyfOffset, headOffset, locaOffset, cmapOffset;
    
    EXTRACT_REVERSED(uint16_t, numTables, 4);
    for(int i=0; i<(int)numTables; ++i)
    {
        EXTRACT_REVERSED(uint32_t, checkSum, 12 + 16 * i + 4);
        EXTRACT_REVERSED(uint32_t, offset, 12 + 16 * i + 8);
        EXTRACT_REVERSED(uint32_t, length, 12 + 16 * i + 12);
        EXTRACT_REVERSED(uint32_t, tableTagInt, 12 + 16 * i);        

//         printf("%4s %x\n", bytes, tableTagInt);
        
        if(tableTagInt == 0x6d617870) maxpOffset = offset;
        else if(tableTagInt == 0x676c7966) glyfOffset = offset;
        else if(tableTagInt == 0x68656164) headOffset = offset;
        else if(tableTagInt == 0x6c6f6361) locaOffset = offset;
        else if(tableTagInt == 0x636d6170) cmapOffset = offset;
    }
    
    EXTRACT_REVERSED(uint16_t, numGlyphs, maxpOffset + 4);
//     printf("Number of glyphs saved in font: %d\n", (int)numGlyphs);
    
    EXTRACT_REVERSED(uint16_t, numEncodingTables, cmapOffset + 2);
    uint32_t windowsUnicodeOffset;
    for(int i=0; i<numEncodingTables; ++i)
    {
        EXTRACT_REVERSED(uint16_t, platformID, cmapOffset + 4 + 8 * i);
        EXTRACT_REVERSED(uint16_t, encodingID, cmapOffset + 4 + 8 * i + 2);
        
        EXTRACT_REVERSED(uint32_t, encodingOffset, cmapOffset + 4 + 8 * i + 4);
        if(platformID == 3 && encodingID == 1) windowsUnicodeOffset = encodingOffset;
        
//         printf("PLAT %d ENC %d\n", platformID, encodingID);
    }
    
    char asciiIndices[128];
    for(int i=0; i<128; ++i) 
    {
        EXTRACT(uint8_t, index, windowsUnicodeOffset + i);
        asciiIndices[i] = index;
        printf("%c is at index: %d\n", (char)i, (int)index);
    }
//     printf("Unicode offset is at: %d\n", windows
    
    EXTRACT_REVERSED(int16_t, indexToLocFormat, headOffset + 50);
    if(indexToLocFormat) // Only parse uint32_t based tables
    {
        for(int i=0; i<128; ++i) 
        {
            EXTRACT_REVERSED(uint32_t, glyphOffset, locaOffset + 4 * i);
            printf(">> Glyph offset of %d (%c) is %d\n", (int)i, (char)i, (int)glyphOffset);
            
            EXTRACT_REVERSED(int16_t, numberOfContours, glyfOffset + glyphOffset);
//             printf("Number of contours in glyph: %d\n", (int)numberOfContours);
            
            //BEGIN FIXME: remove, unnecessary
//             EXTRACT_REVERSED(int16_t, xMin, glyfOffset + glyphOffset + 2);
//             EXTRACT_REVERSED(int16_t, yMin, glyfOffset + glyphOffset + 4);
//             EXTRACT_REVERSED(int16_t, xMax, glyfOffset + glyphOffset + 6);
//             EXTRACT_REVERSED(int16_t, yMax, glyfOffset + glyphOffset + 8);
//             printf("Bounding box: (%d,%d), (%d,%d)\n", xMin, yMin, xMax, yMax);
            //END FIXME
            
            printf("---\n");
            if(numberOfContours >= 0) // Glyph is not composite
            {
                for(int j=0; j<numberOfContours; ++j)
                {
                    EXTRACT_REVERSED(uint16_t, endPointOfContours, glyfOffset + glyphOffset + 10 + 2 * j);
//                     printf("end point: %d\n", endPointOfContours);
                }
                
//                 EXTRACT_REVERSED(uint16_t, instructionLength, glyfOffset + glyphOffset + 10 + 2*numberOfContours);
//                 printf("Number of instrutions: %d\n", (int)instructionLength);
                
                // Extract number of flags and data from contour end point array
                EXTRACT_REVERSED(uint16_t, numberOfPoints, glyfOffset + glyphOffset + 10 + 2 * (numberOfContours-1));
                numberOfPoints += 1;
//                 printf("Number of contained points: %d\n", numberOfPoints);
                
                // Extract size of data from weird compressed format
                uint8_t size = 0, 
                    compressedSize = 0;
                for( ; size < numberOfPoints; ++compressedSize)
                {
                    EXTRACT(uint8_t, flag, glyfOffset + glyphOffset + 10 + 2*numberOfContours + compressedSize);
//                     printBits(flag,8);
                    int8_t repeat = bit(flag, 3);
                    if(repeat)
                    {
//                         printf(" r ");
                        EXTRACT(uint8_t, nRepetitions, glyfOffset + glyphOffset + 10 + 2*numberOfContours + compressedSize + 1);
//                         printBits(nRepetitions, 8);
                        ++compressedSize;
//                         printf(" ");
                        printf("Repetition of size: %d\n", nRepetitions);
                        size += nRepetitions;
                    }
                    else 
                    {
                        size += 1;
                    }
                    
//                     printf("\n");
                }
                printf("actual size: %d, compressed size: %d\n", size, compressedSize);
                
                
                
                
                
//                 if(instructionLength == 0)
//                     printf("----------> nice, no instructions\n");
//                 int nPointsX = 0, nPointsY = 0;
//                 for(int j=0; j<instructionLength; ++j)
//                 {
//                     EXTRACT_REVERSED(uint8_t, instruction, glyfOffset + glyphOffset + 10 + 2*numberOfContours + j);
//                     
//                     int onCurve = (int)bit((uint32_t)instruction, 0);
//                     size_t xWidth = (size_t)bit((uint32_t)instruction, 1) + 1, 
//                         yWidth = (size_t)bit((uint32_t)instruction, 2) + 1;
//                     int repeat = (int)bit((uint32_t)instruction, 3);
//                     if((int)repeat) 
//                     {
//                         for(int k=0; k<8; ++k) printf("%d", bit(instruction,k));
//                         printf("\n");
//                         EXTRACT_REVERSED(uint8_t, nRepetitions, glyfOffset + glyphOffset + 10 + 2*numberOfContours + j + 1);
//                         ++j;
//                         printf("repeat %d times.\n", (int)nRepetitions);
//                     }
//                 }
            }
            printf("\n");
        }
    }
    
    
    return 0;
}

int main(int argc, char **args)
{
    char *texture;
    uint16_t *x = 0, *y = 0;
    uint8_t *flags;
    
    TTF2Texture(x, y, flags, "C:\\Windows\\Fonts\\arial.ttf");
//     TTF2Texture(texture, "C:\\Users\\NR4\\Downloads\\Princess_Sofia\\PrincessSofia-Regular.ttf");
    

    
    
    return 0;
}
