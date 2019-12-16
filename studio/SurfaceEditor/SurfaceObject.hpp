#ifndef SURFACE_OBJECT_HPP
#define SURFACE_OBJECT_HPP

#include "SurfaceNode.hpp"

class SurfaceObject : public SurfaceNode
{
public:
    SurfaceObject(QPointF position = QPointF(0.,0.));
    
    virtual ~SurfaceObject();
    
    void setPosition(const QPointF position);
    const QPointF position() const;
    
private:
    QPointF m_position;
};

#endif
