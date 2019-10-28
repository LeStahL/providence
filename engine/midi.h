#ifndef MIDI_H
#define MIDI_H

#include "common.h"
#include <Windows.h>

UINT nDevices, nOutputDevices;

double faders[8], dials[8];

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
    }
    
	return;
}

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
    }
    
	return;
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
        {
            MMRESULT rv = midiInOpen(&hMidiDevice, i, (DWORD)(void*)MidiInProc_apc40mk2, 0, CALLBACK_FUNCTION);
        }
        else if(!strcmp(capabilities.szPname, "nanoKONTROL2"))
        {
            MMRESULT rv = midiInOpen(&hMidiDevice, i, (DWORD)(void*)MidiInProc_nanoKONTROL2, 0, CALLBACK_FUNCTION);
        }
        midiInStart(hMidiDevice);
    }
    
    nOutputDevices = midiOutGetNumDevs();
    if(nOutputDevices == 0) return;
    
    // Outputs
    MIDIOUTCAPS outCapabilites;
    for(int i=0; i<nOutputDevices; ++i)
    {
        midiOutGetDevCaps(i, &outCapabilites, sizeof(MIDIOUTCAPS));
        
        HMIDIOUT hMidiDevice;
        if(!strcmp(outCapabilites.szPname, "APC40 mkII") || !strcmp(outCapabilites.szPname, "nanoKONTROL2"))
        {
            MMRESULT rv = midiOutOpen(&hMidiDevice, i, 0, 0, CALLBACK_NULL);
        }
    }
}

#endif
