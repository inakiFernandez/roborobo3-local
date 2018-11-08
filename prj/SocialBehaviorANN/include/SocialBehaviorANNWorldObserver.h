/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef SOCIALBEHAVIORANNWORLDOBSERVER_H
#define SOCIALBEHAVIORANNWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNSharedData.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNSharedData.h"

//class World;

class SocialBehaviorANNWorldObserver : public WorldObserver
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
    SocialBehaviorANNWorldObserver(World *world);
    ~SocialBehaviorANNWorldObserver();
    
    std::vector<std::vector<double> > getRobotDistances(){return _robotDistances;}
    void reset();
    void step();
    
    int getGenerationItCount() { return _generationItCount; }
    double  computeGlobalDiversity();
};

#endif
