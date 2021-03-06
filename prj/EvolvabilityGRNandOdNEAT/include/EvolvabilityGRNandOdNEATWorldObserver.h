/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef EVOLVABILITYGRNANDODNEATWORLDOBSERVER_H
#define EVOLVABILITYGRNANDODNEATWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATSharedData.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATSharedData.h"
#include "EvolvabilityGRNandOdNEATController.h"

//class World;

class EvolvabilityGRNandOdNEATWorldObserver : public WorldObserver
{
private:
    void updateEnvironment();
    void updateMonitoring();
    void monitorPopulation( bool localVerbose = true );
    std::vector<std::vector<double> > _robotDistances;
    std::vector<std::map<GCIndividual,Stats > > _offspringStats;
    
protected:
    int _generationCount;
    int _generationItCount;
    
public:
    EvolvabilityGRNandOdNEATWorldObserver(World *world);
    ~EvolvabilityGRNandOdNEATWorldObserver();
    
    std::vector<std::vector<double> > getRobotDistances(){return _robotDistances;}

    std::vector<std::map<GCIndividual,Stats > > getOffspringStats() {return _offspringStats;}
    std::map<GCIndividual,Stats> computeNbOffspringAndFitness();
    void reset();
    void step();
    
    int getGenerationItCount() { return _generationItCount; }
    double  computeGlobalDiversity();
};

#endif
