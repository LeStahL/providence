#ifndef MIDI_H
#define MIDI_H

#include <Windows.h>
#include <math.h>

#define max(a,b) (a)>(b)?(a):(b)
#define min(a,b) (a)<(b)?(a):(b)

UINT nDevices, nOutputDevices;
HMIDIOUT hAPC40mk2, hnanoKONTROL2;
double nanokontrol2_faders[8], 
    nanokontrol2_dials[8], 
    apc40mk2_faders[8],
    apc40mk2_top_dials[8],
    apc40mk2_right_dials[8];

const double apcColors[] = {
    0.00,0.00,0.00,
    0.36,0.18,0.26,
    0.68,0.44,0.58,
    0.78,0.60,0.73,
    0.96,0.02,0.45,
    0.96,0.00,0.00,
    0.92,0.00,0.02,
    0.54,0.00,0.03,
    0.86,0.43,0.46,
    0.98,0.38,0.01,
    0.83,0.06,0.00,
    0.53,0.21,0.00,
    0.99,0.67,0.01,
    0.92,0.64,0.00,
    0.70,0.35,0.00,
    0.41,0.21,0.00,
    0.82,0.78,0.60,
    0.58,0.74,0.00,
    0.51,0.55,0.00,
    0.38,0.36,0.07,
    0.56,0.82,0.60,
    0.02,0.84,0.01,
    0.00,0.56,0.05,
    0.00,0.37,0.10,
    0.53,0.81,0.67,
    0.00,0.85,0.34,
    0.01,0.53,0.27,
    0.00,0.34,0.11,
    0.63,0.74,0.69,
    0.00,0.85,0.60,
    0.00,0.56,0.38,
    0.01,0.33,0.35,
    0.45,0.72,0.79,
    0.00,0.74,0.75,
    0.00,0.40,0.46,
    0.05,0.29,0.49,
    0.74,0.83,1.00,
    0.00,0.72,0.89,
    0.00,0.36,0.54,
    0.00,0.32,0.44,
    0.61,0.62,0.91,
    0.00,0.61,0.95,
    0.00,0.32,0.66,
    0.04,0.11,0.53,
    0.60,0.41,0.84,
    0.16,0.65,0.96,
    0.03,0.06,0.45,
    0.16,0.07,0.39,
    0.51,0.33,0.73,
    0.50,0.09,0.73,
    0.28,0.02,0.39,
    0.34,0.00,0.23,
    0.68,0.17,0.61,
    0.70,0.00,0.69,
    0.65,0.00,0.57,
    0.30,0.00,0.18,
    0.87,0.08,0.41,
    0.97,0.01,0.38,
    0.57,0.00,0.24,
    0.39,0.01,0.20,
    0.98,0.00,0.00,
    0.76,0.04,0.00,
    0.59,0.15,0.00,
    0.59,0.24,0.02,
    0.02,0.31,0.02,
    0.00,0.35,0.30,
    0.00,0.32,0.42,
    0.11,0.25,0.81,
    0.00,0.29,0.40,
    0.39,0.19,0.69,
    0.41,0.20,0.30,
    0.33,0.09,0.25,
    1.00,0.00,0.01,
    0.64,0.51,0.16,
    0.60,0.49,0.00,
    0.35,0.58,0.00,
    0.00,0.49,0.01,
    0.00,0.62,0.65,
    0.00,0.51,0.81,
    0.08,0.35,0.81,
    0.42,0.16,0.75,
    0.73,0.16,0.85,
    0.75,0.00,0.45,
    0.51,0.05,0.00,
    0.95,0.00,0.02,
    0.59,0.48,0.01,
    0.54,0.47,0.07,
    0.00,0.67,0.04,
    0.33,0.65,0.25,
    0.42,0.44,0.31,
    0.33,0.77,0.86,
    0.35,0.33,0.65,
    0.36,0.37,0.62,
    0.47,0.30,0.62,
    0.75,0.06,0.74,
    0.87,0.01,0.46,
    0.81,0.32,0.00,
    0.60,0.23,0.03,
    0.52,0.49,0.01,
    0.72,0.22,0.00,
    0.33,0.05,0.20,
    0.15,0.31,0.24,
    0.00,0.26,0.30,
    0.29,0.14,0.31,
    0.19,0.16,0.46,
    0.46,0.18,0.16,
    0.79,0.00,0.16,
    0.70,0.23,0.32,
    0.95,0.45,0.43,
    0.83,0.36,0.26,
    0.66,0.44,0.22,
    0.55,0.40,0.09,
    0.34,0.16,0.35,
    0.67,0.43,0.36,
    0.51,0.52,0.57,
    0.55,0.32,0.53,
    0.56,0.32,0.60,
    0.43,0.18,0.37,
    0.46,0.18,0.24,
    0.62,0.35,0.58,
    0.79,0.00,0.00,
    0.51,0.01,0.01,
    0.01,0.61,0.09,
    0.00,0.33,0.01,
    0.58,0.26,0.01,
    0.51,0.14,0.01,
    0.99,0.25,0.00,
    0.48,0.04,0.00
};

int apcIndexFromRGB(double r, double g, double b)
{
    int minimalIndex = 0;
    double minimalDiff = 1.;
    
    for(int i=0; i<=0x7f; ++i)
    {
        double diff = max(fabs(r-apcColors[3*i]), fabs(g-apcColors[3*i+1]));
        diff = max(diff, fabs(b-apcColors[3*i+2]));
        
        if(diff < minimalDiff)
        {
            minimalDiff = diff;
            minimalIndex = i;
        }
    }
    
    return minimalIndex;
}

void apcButtonOn(int index, int color)
{
    DWORD msg = 0x9 << 4 | index << 8 | color << 16;
    midiOutShortMsg(hAPC40mk2, msg);
}

void nanokontrol2ButtonOn(int index, int color)
{
    DWORD msg = 0x9 << 4 | index << 8 | color << 16;
    midiOutShortMsg(hnanoKONTROL2, msg);
}

void (*apc40mk2_fader_notifier)(int, double) = 0, 
    (*apc40mk2_dial_top_notifier)(int, double) = 0,
    (*apc40mk2_dial_right_notifier)(int, double) = 0;
void CALLBACK MidiInProc_apc40mk2(HMIDIIN hMidiIn, UINT wMsg, DWORD dwInstance, DWORD dwParam1, DWORD dwParam2)
{
    if(wMsg == MIM_DATA)
    {
        BYTE b1 = (dwParam1 >> 24) & 0xFF,
            b2 = (dwParam1 >> 16) & 0xFF,
            b3 = (dwParam1 >> 8) & 0xFF,
            b4 = dwParam1 & 0xFF;
        BYTE b3lo = b3 & 0xF,
            b3hi = (b3 >> 4) & 0xF,
            b4lo = b4 & 0xF,
            b4hi = (b4 >> 4) & 0xF;
        
        BYTE channel = b4lo,
            button = b3;
        
        printf("Akai APC40 mkII: wMsg=MIM_DATA, dwParam1=%08x, byte=%02x %02x h_%01x l_%01x %02x, dwParam2=%08x\n", dwParam1, b1, b2, b3hi, b3lo, b4, dwParam2);
        
        if(b4 == 0xb0 && b3hi == 0x3) // top Dial
        {
            apc40mk2_top_dials[b3lo] = (double)b2/(double)0x7f;
            if(apc40mk2_dial_top_notifier != 0) (*apc40mk2_dial_top_notifier)(b3lo, apc40mk2_top_dials[b3lo]);
        }
        else if(b4 == 0xb0 && b3hi == 0x1) // right Dial
        {
            apc40mk2_right_dials[b3lo] = (double)b2/(double)0x7f;
            if(apc40mk2_dial_top_notifier != 0) (*apc40mk2_dial_right_notifier)(b3lo, apc40mk2_right_dials[b3lo]);
        }
        else if(b4hi == 0xb && b3 == 0x07)
        {
            apc40mk2_faders[b4lo] = (double)b2/(double)0x7f;
            if(apc40mk2_fader_notifier != 0) (*apc40mk2_fader_notifier)(b4lo, apc40mk2_faders[b4lo]);
//             printf("%d\n", b3lo
        }
    }
    
	return;
}

void (*nanokontrol2_fader_notifier)(int, double) = 0, (*nanokontrol2_dial_notifier)(int, double) = 0;
void CALLBACK MidiInProc_nanoKONTROL2(HMIDIIN hMidiIn, UINT wMsg, DWORD dwInstance, DWORD dwParam1, DWORD dwParam2)
{
    if(wMsg == MIM_DATA)
    {
        BYTE b1 = (dwParam1 >> 24) & 0xFF,
            b2 = (dwParam1 >> 16) & 0xFF,
            b3 = (dwParam1 >> 8) & 0xFF,
            b4 = dwParam1 & 0xFF;
        BYTE b3lo = b3 & 0xF,
            b3hi = (b3 >> 4) & 0xF,
            b4lo = b4 & 0xF,
            b4hi = (b4 >> 4) & 0xF;
        
        BYTE channel = b4lo,
            button = b3;
        
        printf("KORG nanoKONTROL2: wMsg=MIM_DATA, dwParam1=%08x, byte=%02x %02x h_%01x l_%01x %02x, dwParam2=%08x\n", dwParam1, b1, b2, b3hi, b3lo, b4, dwParam2);
        
        if(b4 == 0xb0) // Fader or dial
        {
            if(b3hi == 0) // Fader
            {
                nanokontrol2_faders[b3lo] = (double)b2/(double)0x7f;
                if(nanokontrol2_fader_notifier != 0) (*nanokontrol2_fader_notifier)(b3lo, nanokontrol2_faders[b3lo]);
            }
            else if(b3hi == 1) // Dial
            {
                nanokontrol2_dials[b3lo] = (double)b2/(double)0x7f;
                if(nanokontrol2_dial_notifier != 0) (*nanokontrol2_dial_notifier)(b3lo, nanokontrol2_dials[b3lo]);
            }
        }
    }
    
	return;
}

void initKorgNanoKontrol2Input(void *callback)
{
    // Number of devices
    nDevices = midiInGetNumDevs();
    if(nDevices == 0) return;
    
    // Inputs
    MIDIINCAPS capabilities;
    for(int i=0; i<nDevices; ++i)
    {
        midiInGetDevCaps(i, &capabilities, sizeof(MIDIINCAPS));
        
        HMIDIIN hMidiDevice;
        if(!strcmp(capabilities.szPname, "nanoKONTROL2"))
            midiInOpen(&hMidiDevice, i, (DWORD_PTR)callback, 0, CALLBACK_FUNCTION);
        midiInStart(hMidiDevice);
    }
}

void initApc40Mk2Input(void *callback)
{
    // Number of devices
    nDevices = midiInGetNumDevs();
    if(nDevices == 0) return;
    
    // Inputs
    MIDIINCAPS capabilities;
    for(int i=0; i<nDevices; ++i)
    {
        midiInGetDevCaps(i, &capabilities, sizeof(MIDIINCAPS));
        
        HMIDIIN hMidiDevice;
        if(!strcmp(capabilities.szPname, "APC40 mkII"))
            midiInOpen(&hMidiDevice, i, (DWORD_PTR)callback, 0, CALLBACK_FUNCTION);
        midiInStart(hMidiDevice);
    }
}

void initControllers()
{
    // Number of devices
    nDevices = midiInGetNumDevs();
    if(nDevices == 0) return;
    
    // Inputs
    MIDIINCAPS capabilities;
    for(int i=0; i<nDevices; ++i)
    {
        midiInGetDevCaps(i, &capabilities, sizeof(MIDIINCAPS));
        
        HMIDIIN hMidiDevice;
        if(!strcmp(capabilities.szPname, "APC40 mkII"))
            midiInOpen(&hMidiDevice, i, (DWORD)(void*)MidiInProc_apc40mk2, 0, CALLBACK_FUNCTION);
        else if(!strcmp(capabilities.szPname, "nanoKONTROL2"))
            midiInOpen(&hMidiDevice, i, (DWORD)(void*)MidiInProc_nanoKONTROL2, 0, CALLBACK_FUNCTION);
        midiInStart(hMidiDevice);
    }
    
    nOutputDevices = midiOutGetNumDevs();
    if(nOutputDevices == 0) return;
    
    // Outputs
    MIDIOUTCAPS outCapabilites;
    for(int i=0; i<nOutputDevices; ++i)
    {
        midiOutGetDevCaps(i, &outCapabilites, sizeof(MIDIOUTCAPS));
        
        if(!strcmp(outCapabilites.szPname, "APC40 mkII"))
        {
            midiOutOpen(&hAPC40mk2, i, 0, 0, CALLBACK_NULL);
        }
        else if(!strcmp(outCapabilites.szPname, "nanoKONTROL2"))
        {
            midiOutOpen(&hnanoKONTROL2, i, 0, 0, CALLBACK_NULL);
        }
    }
}

#endif
