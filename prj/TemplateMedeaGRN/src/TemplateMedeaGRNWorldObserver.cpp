/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNWorldObserver.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNController.h"
#include "World/World.h"

TemplateMedeaGRNWorldObserver::TemplateMedeaGRNWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;
    
    // ==== loading project-specific properties
    
    gProperties.checkAndGetPropertyValue("gSigmaRef",&TemplateMedeaGRNSharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gSigmaMin",&TemplateMedeaGRNSharedData::gSigmaMin,true);
    gProperties.checkAndGetPropertyValue("gSigmaMax",&TemplateMedeaGRNSharedData::gSigmaMax,true);
    
    gProperties.checkAndGetPropertyValue("gProbaMutation",&TemplateMedeaGRNSharedData::gProbaMutation,true);
    gProperties.checkAndGetPropertyValue("gUpdateSigmaStep",&TemplateMedeaGRNSharedData::gUpdateSigmaStep,true);
    gProperties.checkAndGetPropertyValue("gEvaluationTime",&TemplateMedeaGRNSharedData::gEvaluationTime,true);
    gProperties.checkAndGetPropertyValue("gSynchronization",&TemplateMedeaGRNSharedData::gSynchronization,true);
    
    gProperties.checkAndGetPropertyValue("gEnergyRequestOutput",&TemplateMedeaGRNSharedData::gEnergyRequestOutput,false);
    
    gProperties.checkAndGetPropertyValue("gMonitorPositions",&TemplateMedeaGRNSharedData::gMonitorPositions,true);
    
    
    gProperties.checkAndGetPropertyValue("gSnapshots",&TemplateMedeaGRNSharedData::gSnapshots,false);
    gProperties.checkAndGetPropertyValue("gSnapshotsFrequency",&TemplateMedeaGRNSharedData::gSnapshotsFrequency,false);
    
    gProperties.checkAndGetPropertyValue("gControllerType",&TemplateMedeaGRNSharedData::gControllerType,true);
    
    gProperties.checkAndGetPropertyValue("gMaxNbGenomeTransmission",&TemplateMedeaGRNSharedData::gMaxNbGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gLimitGenomeTransmission",&TemplateMedeaGRNSharedData::gLimitGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gSelectionMethod",&TemplateMedeaGRNSharedData::gSelectionMethod,true);
    
    
    gProperties.checkAndGetPropertyValue("gLogGenome",&TemplateMedeaGRNSharedData::gLogGenome,false);
    
    gProperties.checkAndGetPropertyValue("gIndividualMutationRate",&TemplateMedeaGRNSharedData::gIndividualMutationRate,false);

    gProperties.checkAndGetPropertyValue("gMutationOperator",&TemplateMedeaGRNSharedData::gMutationOperator,false);
    gProperties.checkAndGetPropertyValue("gFitness",&TemplateMedeaGRNSharedData::gFitness,true);

    gProperties.checkAndGetPropertyValue("gSelPressure",&TemplateMedeaGRNSharedData::gSelPressure,true);


    gProperties.checkAndGetPropertyValue("gMutationOperator",&TemplateMedeaGRNSharedData::gMutateProb,true);
    gProperties.checkAndGetPropertyValue("gFitness",&TemplateMedeaGRNSharedData::gCrossoverProb,true);

    gProperties.checkAndGetPropertyValue("gCommunicationRange",&TemplateMedeaGRNSharedData::gCommunicationRange,true);
    gProperties.checkAndGetPropertyValue("gCommunicationOnRadius",&TemplateMedeaGRNSharedData::gCommunicationOnRadius,true);

    gProperties.checkAndGetPropertyValue("gMatingOperator",&TemplateMedeaGRNSharedData::gMatingOperator,true);
    gProperties.checkAndGetPropertyValue("gBroadcastTime",&TemplateMedeaGRNSharedData::gBroadcastTime,true);

    if ( !gRadioNetwork)
    {
        std::cout << "Error : gRadioNetwork must be true." << std::endl;
        exit(-1);
    }
    
    // * iteration and generation counters
    _generationItCount = -1;
    _generationCount = -1;
}

TemplateMedeaGRNWorldObserver::~TemplateMedeaGRNWorldObserver()
{
    // nothing to do.
}

void TemplateMedeaGRNWorldObserver::reset()
{
}

void TemplateMedeaGRNWorldObserver::step()
{
    _generationItCount++;
    
    if( _generationItCount == TemplateMedeaGRNSharedData::gEvaluationTime+1 ) // switch to next generation.
    {
        // update iterations and generations counters
        _generationItCount = 0;
        _generationCount++;
    }
    //Local broadcast with distance: compute inter-robot distance matrix
    if(TemplateMedeaGRNSharedData::gCommunicationOnRadius)
    {
    _robotDistances = std::vector<std::vector<double> >(gNumberOfRobots, std::vector<double>(gNumberOfRobots,-1));
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {
            TemplateMedeaGRNController* c = (dynamic_cast<TemplateMedeaGRNController*>(gWorld->getRobot(i)->getController()));
            double x = c->getWorldModel()->getXReal();
            double y = c->getWorldModel()->getYReal();

            for ( int j = i ; j != gNumberOfRobots ; j++ )
            {
                TemplateMedeaGRNController* c2 = (dynamic_cast<TemplateMedeaGRNController*>(gWorld->getRobot(j)->getController()));
                double x2 = c2->getWorldModel()->getXReal();
                double y2 = c2->getWorldModel()->getYReal();
                _robotDistances[i][j] = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)); //.push_back(sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)));
                _robotDistances[j][i] = _robotDistances[i][j];
            }

        }
    }
    updateMonitoring();
    
    updateEnvironment();
    
}


void TemplateMedeaGRNWorldObserver::updateEnvironment()
{
    // example: moving landmarks
    /*
    if ( gWorld->getIterations() % 2000 == 0 )
        for ( int i = 0 ; i != gLandmarks.size() ; i++ )
        {
            Point2d* position = new Point2d( 200+rand()%(gAreaWidth-400) , 200+rand()%(gAreaHeight-400) );
            gLandmarks[i].setPosition(*position);
        }
    */
}

void TemplateMedeaGRNWorldObserver::updateMonitoring()
{
    // * Log at end of each generation

    if( //gWorld->getIterations() % TemplateMedeaGRNSharedData::gEvaluationTime == 1 ||
            gWorld->getIterations() % TemplateMedeaGRNSharedData::gEvaluationTime == TemplateMedeaGRNSharedData::gEvaluationTime-1 ) // beginning(+1) *and* end of generation. ("==1" is required to monitor the outcome of the first iteration)
    {
        monitorPopulation();
    }
    
    // * Every N generations, take a video (duration: one generation time)
    
    if ( TemplateMedeaGRNSharedData::gSnapshots )
    {
        if ( ( gWorld->getIterations() ) % ( TemplateMedeaGRNSharedData::gEvaluationTime * TemplateMedeaGRNSharedData::gSnapshotsFrequency ) == 0 )
        {
            if ( gVerbose )
                std::cout << "[START] Video recording: generation #" << (gWorld->getIterations() / TemplateMedeaGRNSharedData::gEvaluationTime ) << ".\n";
            gTrajectoryMonitorMode = 0;
            initTrajectoriesMonitor();
        }
        else
            if ( ( gWorld->getIterations() ) % ( TemplateMedeaGRNSharedData::gEvaluationTime * TemplateMedeaGRNSharedData::gSnapshotsFrequency ) == TemplateMedeaGRNSharedData::gEvaluationTime - 1 )
            {
                if ( gVerbose )
                    std::cout << "[STOP]  Video recording: generation #" << (gWorld->getIterations() / TemplateMedeaGRNSharedData::gEvaluationTime ) << ".\n";
                saveTrajectoryImage();
            }
    }    
}

void TemplateMedeaGRNWorldObserver::monitorPopulation( bool localVerbose )
{
    // * monitoring: count number of active agents.
    
    int activeCount = 0;
    for ( int i = 0 ; i != gNumberOfRobots ; i++ )
    {
        if ( (dynamic_cast<TemplateMedeaGRNController*>(gWorld->getRobot(i)->getController()))->getWorldModel()->isAlive() == true )
            activeCount++;
    }
    
    if ( gVerbose && localVerbose )
    {
        std::cout << "[gen:" << (gWorld->getIterations()/TemplateMedeaGRNSharedData::gEvaluationTime) << ";it:" << gWorld->getIterations() << ";pop:" << activeCount << "]\n";
    }
    // Logging here
    double sumFitness = 0.0;
    double sumItems = 0.0;
    double sumCollisions = 0.0;
    int gatheredGenomes = 0;
    for ( int i = 0 ; i != gNumberOfRobots ; i++ )
    {

        TemplateMedeaGRNController* c =(dynamic_cast<TemplateMedeaGRNController*>(gWorld->getRobot(i)->getController()));
         sumFitness += c->getFitness();
         sumItems += c->getCollectedItems();
         sumCollisions += c->getNbCollisions();
         gatheredGenomes += c->getGenomesList().size();
    }
    switch (TemplateMedeaGRNSharedData::gFitness)
    {
        case 1:
        {
            sumFitness *= TemplateMedeaGRNSharedData::gEvaluationTime;
            break;
        }
    }

    std::cout //<< "[gen:"
            << (gWorld->getIterations()/TemplateMedeaGRNSharedData::gEvaluationTime)
            //<< ", fitness: "
            << " "
            << sumFitness / gNumberOfRobots
            //<< ", gatheredGenomes: "
            << " "
            << sumItems / gNumberOfRobots
            << " "
            << sumCollisions / gNumberOfRobots
            << " "
            << gatheredGenomes
            << std::endl;
    // Logging, population-level: alive
    std::string sLog = std::string("") + std::to_string(gWorld->getIterations()) + ",pop,alive," + std::to_string(activeCount) + "\n";
    gLogManager->write(sLog);
    gLogManager->flush();
}
