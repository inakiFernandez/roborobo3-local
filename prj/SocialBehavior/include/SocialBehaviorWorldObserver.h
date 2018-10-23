/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef SOCIALBEHAVIORWORLDOBSERVER_H
#define SOCIALBEHAVIORWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "SocialBehavior/include/SocialBehaviorSharedData.h"
#include "SocialBehavior/include/SocialBehaviorSharedData.h"

//class World;

class SocialBehaviorWorldObserver : public WorldObserver
{
private:
    void updateEnvironment();
    void updateMonitoring();
    void monitorPopulation( bool localVerbose = true );
    
protected:
    int _generationCount;
    int _generationItCount;
    
public:
    SocialBehaviorWorldObserver(World *world);
    ~SocialBehaviorWorldObserver();
    
    void reset();
    void step();
    
    int getGenerationItCount() { return _generationItCount; }

};

#endif
