#ifndef SURFACE_SCENE_HPP
#define SURFACE_SCENE_HPP

#include <QList>

#include "SurfaceNode.hpp"

class SurfaceScene 
{
public:
    SurfaceScene();
    virtual ~SurfaceScene();
    
private:
    QList<SurfaceNode *> m_contained_nodes;
};

#endif
