#ifndef NANOKONTROL2WIDGET_H
#define NANOKONTROL2WIDGET_H

#include "ui_nanokontrol2.h"

class nanoKontrol2Widget : public QWidget
{
public:
    nanoKontrol2Widget(QWidget *parent);
    virtual ~nanoKontrol2Widget();
    
private:
    Ui::nanoKontrol2Widget m_ui;
};

#endif
