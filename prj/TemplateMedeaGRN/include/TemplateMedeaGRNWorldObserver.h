/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef TEMPLATEMEDEAGRNWORLDOBSERVER_H
#define TEMPLATEMEDEAGRNWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNSharedData.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNSharedData.h"

//class World;

class TemplateMedeaGRNWorldObserver : public WorldObserver
{
private:
    void updateEnvironment();
    void updateMonitoring();
    void monitorPopulation( bool localVerbose = true );
    std::vector<std::vector<double> > _robotDistances;
    
protected:
    int _generationCount;
    int _generationItCount;
    
public:
    TemplateMedeaGRNWorldObserver(World *world);
    ~TemplateMedeaGRNWorldObserver();
    
    std::vector<std::vector<double> > getRobotDistances(){return _robotDistances;}
    void reset();
    void step();
    
    int getGenerationItCount() { return _generationItCount; }

};

#endif
