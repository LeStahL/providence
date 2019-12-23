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
    // Read file to memory
    FILE *f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    size_t fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *fdata = (char*)malloc(fsize+1);
    fread(fdata, 1, fsize, f);
    fclose(f);
    
    // Extract relevant table offsets
    uint32_t maxpOffset, glyfOffset, headOffset, locaOffset, cmapOffset;
    EXTRACT_REVERSED(uint16_t, numTables, 4);
    for(int i=0; i<(int)numTables; ++i)
    {
        EXTRACT_REVERSED(uint32_t, checkSum, 12 + 16 * i + 4);
        EXTRACT_REVERSED(uint32_t, offset, 12 + 16 * i + 8);
        EXTRACT_REVERSED(uint32_t, length, 12 + 16 * i + 12);
        EXTRACT_REVERSED(uint32_t, tableTagInt, 12 + 16 * i);        

        if(tableTagInt == 0x6d617870) maxpOffset = offset;
        else if(tableTagInt == 0x676c7966) glyfOffset = offset;
        else if(tableTagInt == 0x68656164) headOffset = offset;
        else if(tableTagInt == 0x6c6f6361) locaOffset = offset;
        else if(tableTagInt == 0x636d6170) cmapOffset = offset;
    }
    
    // Extract number of glyphs
    EXTRACT_REVERSED(uint16_t, numGlyphs, maxpOffset + 4);

    // Extract encoding and platform information
    EXTRACT_REVERSED(uint16_t, numEncodingTables, cmapOffset + 2);
    uint32_t windowsUnicodeOffset;
    for(int i=0; i<numEncodingTables; ++i)
    {
        EXTRACT_REVERSED(uint16_t, platformID, cmapOffset + 4 + 8 * i);
        EXTRACT_REVERSED(uint16_t, encodingID, cmapOffset + 4 + 8 * i + 2);
        
        EXTRACT_REVERSED(uint32_t, encodingOffset, cmapOffset + 4 + 8 * i + 4);
        
        printf("Contained platform ID: %u, encoding ID: %u\n", platformID, encodingID);
        
        // Windows + Unicode is the way to go
        if(platformID == 3 && encodingID == 1) windowsUnicodeOffset = encodingOffset;
    }
    
    // Extract CMAP subtable information
    EXTRACT_REVERSED(uint16_t, segCountX2, cmapOffset + windowsUnicodeOffset + 6);
    uint16_t segCount = segCountX2 / 2;
    printf("Number of segments: %u\n", segCount);
    
    // Extract range intervals
    uint16_t *startCode = (uint16_t*) malloc(segCount*sizeof(uint16_t)),
        *endCode = (uint16_t*) malloc(segCount*sizeof(uint16_t)), 
        *idRangeOffset = (uint16_t*) malloc(segCount*sizeof(uint16_t));
    for(int i=0; i<segCount; ++i)
    {
        EXTRACT_REVERSED(uint16_t, endCodeI, cmapOffset + windowsUnicodeOffset + 14 + 2 * i);
        endCode[i] = endCodeI;
        EXTRACT_REVERSED(uint16_t, startCodeI, cmapOffset + windowsUnicodeOffset + 16 + segCountX2 + 2 * i);
        startCode[i] = startCodeI;
        EXTRACT_REVERSED(uint16_t, idRangeOffsetI, cmapOffset + windowsUnicodeOffset + 16 + 3 * segCountX2 + 2 * i);
        idRangeOffset[i] = idRangeOffsetI;
    }
    //14 :: end code
    
//     EXTRACT_REVERSED(uint16_t, cmapFormat, cmapOffset + windowsUnicodeOffset);
//     printf("CMAP subtable 3/1 format is %u\n", cmapFormat);
//     
//     EXTRACT_REVERSED(uint16_t, subtableSize, cmapOffset + windowsUnicodeOffset + 2);
//     printf("CMAP subtable has length %u bytes.\n", subtableSize);
    
    // Determine the glyph indices of all relevant ascii chars from CMAP
    uint16_t asciiIndices[128];
    for(int i=32; i<127; ++i) 
    {
        EXTRACT_REVERSED(uint16_t, index, cmapOffset + windowsUnicodeOffset + 24 + 2*i);
//         asciiIndices[i] = index;
        int segment = -1;
        for(int j=0; j<segCount; ++j)
        {
            if(i > startCode[j] && i < endCode[j]) 
            {
                segment = j;
                printf("Found segment index: %d\n", j);
                break;
            }
        }
        if(segment == -1) printf("glyph data for %c not present in font file.\n", (char)i);
        
        asciiIndices[i] = *(idRangeOffset[segment]/2
            + (i - startCode[segment])
            + &idRangeOffset[segment]);
        
        printf("%c / %d -> index %u\n", (char)i, i, (unsigned int)asciiIndices[i]);
    }
    
    EXTRACT_REVERSED(int16_t, indexToLocFormat, headOffset + 50);
    if(indexToLocFormat == 1) // Only parse uint32_t based tables
    {
        for(int i=0; i<128; ++i) 
        {
            printf("\nChar: %d (%c)\n", i, (char)(i=='\n'?'n':i));
            
            // Convert glyph index to GLYF offset with LOCA table
            EXTRACT_REVERSED(uint32_t, glyphOffset, locaOffset + 4*asciiIndices[i]);
            printf("Glyph offset: %u\n", glyphOffset);
            
            // Extract actual glyph data from GLYF
            // First determine the number of contours
            EXTRACT_REVERSED(int16_t, numberOfContours, glyfOffset + glyphOffset);
            printf("Contours: %d @ ", numberOfContours);
            if(numberOfContours >= 0) // Glyph is not composite
            {
                int *endpoints = (int *)malloc(numberOfContours*sizeof(int));
                for(int j=0; j<numberOfContours; ++j)
                {
                    EXTRACT_REVERSED(uint16_t, endPointOfContours, glyfOffset + glyphOffset + 10 + 2 * j);
                    endpoints[j] = endPointOfContours;
                    printf("%d ", endpoints[j]);
                }
                printf("\n");
                
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
//                         printf("Repetition of size: %d\n", nRepetitions);
                        size += nRepetitions;
                    }
                    else 
                    {
                        size += 1;
                    }
                    
//                     printf("\n");
                }
//                 printf("actual size: %d, compressed size: %d\n", size, compressedSize);
                
                
                
                
                
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
//             printf("\n");
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
