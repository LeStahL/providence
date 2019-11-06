#ifndef NANOKONTROL2WIDGET_H
#define NANOKONTROL2WIDGET_H

#include "ui_nanokontrol2.h"

class nanoKontrol2Widget : public QWidget
{
public:
    nanoKontrol2Widget(QWidget *parent);
    virtual ~nanoKontrol2Widget();

    static void setFader(int index, double value);
    static void setDial(int index, double value);
    
signals:
    void faderChanged(int index, double value);
    void dialChanged(int index, double value);
    
private:
    Ui::nanoKontrol2Widget m_ui;
};

#endif
