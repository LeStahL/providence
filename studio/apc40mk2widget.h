#ifndef APC40MK2WIDGET_H
#define APC40MK2WIDGET_H

#include "ui_apc40mk2widget.h"

class apc40Mk2Widget : public QWidget
{
public:
    apc40Mk2Widget(QWidget *parent);
    virtual ~apc40Mk2Widget();

    static void setFader(int index, double value);
    static void setDial(int index, double value);
    static void setRightDial(int index, double value);
    
signals:
    void faderChanged(int index, double value);
    void dialChanged(int index, double value);
    
private:
    Ui::apc40Mk2Widget m_ui;
};

#endif
