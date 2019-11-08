#include "apc40mk2widget.h"

apc40Mk2Widget *apc40_instance;

#ifdef WIN32
#define STUDIO
extern "C"
{
    #include "engine/midi.h"
}
#endif

#include <QDebug>

apc40Mk2Widget::apc40Mk2Widget(QWidget* parent)
{
    m_ui.setupUi(this);
    
#ifdef WIN32
    apc40mk2_fader_notifier = &apc40_instance->setFader;
    apc40mk2_dial_top_notifier = &apc40_instance->setDial;
    apc40mk2_dial_right_notifier = &apc40_instance->setRightDial;
    initApc40Mk2Input((void*)&MidiInProc_apc40mk2);
#endif
    
    apc40_instance = this;
}

apc40Mk2Widget::~apc40Mk2Widget()
{
}

void apc40Mk2Widget::setFader(int index, double value)
{
    emit apc40_instance->faderChanged(index, value);
}

void apc40Mk2Widget::faderChanged(int index, double value)
{
    if(index == 0) m_ui.verticalSlider_1->setValue((int)(value*1000.));
    else if(index == 1) m_ui.verticalSlider_2->setValue((int)(value*1000.));
    else if(index == 2) m_ui.verticalSlider_3->setValue((int)(value*1000.));
    else if(index == 3) m_ui.verticalSlider_4->setValue((int)(value*1000.));
    else if(index == 4) m_ui.verticalSlider_5->setValue((int)(value*1000.));
    else if(index == 5) m_ui.verticalSlider_6->setValue((int)(value*1000.));
    else if(index == 6) m_ui.verticalSlider_7->setValue((int)(value*1000.));
    else if(index == 7) m_ui.verticalSlider_8->setValue((int)(value*1000.));
    
    m_ui.horizontalLayout_2->update();
}

void apc40Mk2Widget::setDial(int index, double value)
{
    emit apc40_instance->dialChanged(index, value);
}

void apc40Mk2Widget::setRightDial(int index, double value)
{
    emit apc40_instance->dialChanged(index+8, value);
}

void apc40Mk2Widget::dialChanged(int index, double value)
{
    if(index == 0) m_ui.dial_1->setValue((int)(value*1000.));
    else if(index == 1) m_ui.dial_2->setValue((int)(value*1000.));
    else if(index == 2) m_ui.dial_3->setValue((int)(value*1000.));
    else if(index == 3) m_ui.dial_4->setValue((int)(value*1000.));
    else if(index == 4) m_ui.dial_5->setValue((int)(value*1000.));
    else if(index == 5) m_ui.dial_6->setValue((int)(value*1000.));
    else if(index == 6) m_ui.dial_7->setValue((int)(value*1000.));
    else if(index == 7) m_ui.dial_8->setValue((int)(value*1000.));
    else if(index == 8) m_ui.dial_9->setValue((int)(value*1000.));
    else if(index == 9) m_ui.dial_10->setValue((int)(value*1000.));
    else if(index == 10) m_ui.dial_11->setValue((int)(value*1000.));
    else if(index == 11) m_ui.dial_12->setValue((int)(value*1000.));
    else if(index == 12) m_ui.dial_13->setValue((int)(value*1000.));
    else if(index == 13) m_ui.dial_14->setValue((int)(value*1000.));
    else if(index == 14) m_ui.dial_15->setValue((int)(value*1000.));
    else if(index == 15) m_ui.dial_16->setValue((int)(value*1000.));
    m_ui.horizontalLayout_2->update();
}
