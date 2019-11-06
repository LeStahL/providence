#include "nanokontrol2Widget.h"

nanoKontrol2Widget *instance;

#ifdef WIN32
#define STUDIO
#include "engine/midi.h"
#endif

#include <QDebug>


nanoKontrol2Widget::nanoKontrol2Widget(QWidget* parent)
{
    m_ui.setupUi(this);
    
#ifdef WIN32
    fader_notifier = &instance->setFader;
    dial_notifier = &instance->setDial;
    initKorgNanoKontrol2Input((void*)&MidiInProc_nanoKONTROL2);
#endif
    
    instance = this;
}

nanoKontrol2Widget::~nanoKontrol2Widget()
{
}

void nanoKontrol2Widget::setFader(int index, double value)
{
    emit instance->faderChanged(index, value);
}

void nanoKontrol2Widget::faderChanged(int index, double value)
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

void nanoKontrol2Widget::setDial(int index, double value)
{
    emit instance->dialChanged(index, value);
}

void nanoKontrol2Widget::dialChanged(int index, double value)
{
    if(index == 0) m_ui.dial_1->setValue((int)(value*1000.));
    else if(index == 1) m_ui.dial_2->setValue((int)(value*1000.));
    else if(index == 2) m_ui.dial_3->setValue((int)(value*1000.));
    else if(index == 3) m_ui.dial_4->setValue((int)(value*1000.));
    else if(index == 4) m_ui.dial_5->setValue((int)(value*1000.));
    else if(index == 5) m_ui.dial_6->setValue((int)(value*1000.));
    else if(index == 6) m_ui.dial_7->setValue((int)(value*1000.));
    else if(index == 7) m_ui.dial_8->setValue((int)(value*1000.));
    
    m_ui.horizontalLayout_2->update();
}
