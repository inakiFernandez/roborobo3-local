/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "SocialBehavior/include/SocialBehaviorWorldObserver.h"
#include "SocialBehavior/include/SocialBehaviorController.h"
#include "World/World.h"

SocialBehaviorWorldObserver::SocialBehaviorWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;
    
    // ==== loading project-specific properties
    
    gProperties.checkAndGetPropertyValue("gSigmaRef",&SocialBehaviorSharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gSigmaMin",&SocialBehaviorSharedData::gSigmaMin,true);
    gProperties.checkAndGetPropertyValue("gSigmaMax",&SocialBehaviorSharedData::gSigmaMax,true);
    
    gProperties.checkAndGetPropertyValue("gProbaMutation",&SocialBehaviorSharedData::gProbaMutation,true);
    gProperties.checkAndGetPropertyValue("gUpdateSigmaStep",&SocialBehaviorSharedData::gUpdateSigmaStep,true);
    gProperties.checkAndGetPropertyValue("gEvaluationTime",&SocialBehaviorSharedData::gEvaluationTime,true);
    gProperties.checkAndGetPropertyValue("gSynchronization",&SocialBehaviorSharedData::gSynchronization,true);
    
    gProperties.checkAndGetPropertyValue("gEnergyRequestOutput",&SocialBehaviorSharedData::gEnergyRequestOutput,false);
    
    gProperties.checkAndGetPropertyValue("gMonitorPositions",&SocialBehaviorSharedData::gMonitorPositions,true);
    
    gProperties.checkAndGetPropertyValue("gNbHiddenLayers",&SocialBehaviorSharedData::gNbHiddenLayers,true);
    gProperties.checkAndGetPropertyValue("gNbNeuronsPerHiddenLayer",&SocialBehaviorSharedData::gNbNeuronsPerHiddenLayer,true);
    gProperties.checkAndGetPropertyValue("gNeuronWeightRange",&SocialBehaviorSharedData::gNeuronWeightRange,true);
    
    gProperties.checkAndGetPropertyValue("gSnapshots",&SocialBehaviorSharedData::gSnapshots,false);
    gProperties.checkAndGetPropertyValue("gSnapshotsFrequency",&SocialBehaviorSharedData::gSnapshotsFrequency,false);
    
    gProperties.checkAndGetPropertyValue("gControllerType",&SocialBehaviorSharedData::gControllerType,true);
    
    gProperties.checkAndGetPropertyValue("gMaxNbGenomeTransmission",&SocialBehaviorSharedData::gMaxNbGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gLimitGenomeTransmission",&SocialBehaviorSharedData::gLimitGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gSelectionMethod",&SocialBehaviorSharedData::gSelectionMethod,true);
    
    gProperties.checkAndGetPropertyValue("gNotListeningStateDelay",&SocialBehaviorSharedData::gNotListeningStateDelay,true);
    gProperties.checkAndGetPropertyValue("gListeningStateDelay",&SocialBehaviorSharedData::gListeningStateDelay,true);
    
    gProperties.checkAndGetPropertyValue("gLogGenome",&SocialBehaviorSharedData::gLogGenome,false);
    
    gProperties.checkAndGetPropertyValue("gIndividualMutationRate",&SocialBehaviorSharedData::gIndividualMutationRate,false);

    gProperties.checkAndGetPropertyValue("gMutationOperator",&SocialBehaviorSharedData::gMutationOperator,false);
    
    gProperties.checkAndGetPropertyValue("gSigma",&SocialBehaviorSharedData::gSigma,false);
    
    // ====
    
    if ( !gRadioNetwork)
    {
        std::cout << "Error : gRadioNetwork must be true." << std::endl;
        exit(-1);
    }
    
    // * iteration and generation counters
    
    _generationItCount = -1;
    _generationCount = -1;
}

SocialBehaviorWorldObserver::~SocialBehaviorWorldObserver()
{
    // nothing to do.
}

void SocialBehaviorWorldObserver::reset()
{
}

void SocialBehaviorWorldObserver::step()
{
    _generationItCount++;
    
    if( _generationItCount == SocialBehaviorSharedData::gEvaluationTime+1 ) // switch to next generation.
    {
        // update iterations and generations counters
        _generationItCount = 0;
        _generationCount++;
    }
    
    updateMonitoring();
    
    updateEnvironment();
    
}


void SocialBehaviorWorldObserver::updateEnvironment()
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

void SocialBehaviorWorldObserver::updateMonitoring()
{
    // * Log at end of each generation

    if( gWorld->getIterations() % SocialBehaviorSharedData::gEvaluationTime == 1 || gWorld->getIterations() % SocialBehaviorSharedData::gEvaluationTime == SocialBehaviorSharedData::gEvaluationTime-1 ) // beginning(+1) *and* end of generation. ("==1" is required to monitor the outcome of the first iteration)
    {
        monitorPopulation();
    }
    
    // * Every N generations, take a video (duration: one generation time)
    
    if ( SocialBehaviorSharedData::gSnapshots )
    {
        if ( ( gWorld->getIterations() ) % ( SocialBehaviorSharedData::gEvaluationTime * SocialBehaviorSharedData::gSnapshotsFrequency ) == 0 )
        {
            if ( gVerbose )
                std::cout << "[START] Video recording: generation #" << (gWorld->getIterations() / SocialBehaviorSharedData::gEvaluationTime ) << ".\n";
            gTrajectoryMonitorMode = 0;
            initTrajectoriesMonitor();
        }
        else
            if ( ( gWorld->getIterations() ) % ( SocialBehaviorSharedData::gEvaluationTime * SocialBehaviorSharedData::gSnapshotsFrequency ) == SocialBehaviorSharedData::gEvaluationTime - 1 )
            {
                if ( gVerbose )
                    std::cout << "[STOP]  Video recording: generation #" << (gWorld->getIterations() / SocialBehaviorSharedData::gEvaluationTime ) << ".\n";
                saveTrajectoryImage();
            }
    }    
}

void SocialBehaviorWorldObserver::monitorPopulation( bool localVerbose )
{
    // * monitoring: count number of active agents.
    
    int activeCount = 0;
    for ( int i = 0 ; i != gNumberOfRobots ; i++ )
    {
        if ( (dynamic_cast<SocialBehaviorController*>(gWorld->getRobot(i)->getController()))->getWorldModel()->isAlive() == true )
            activeCount++;
    }
    
    if ( gVerbose && localVerbose )
    {
        std::cout << "[gen:" << (gWorld->getIterations()/SocialBehaviorSharedData::gEvaluationTime) << ";it:" << gWorld->getIterations() << ";pop:" << activeCount << "]\n";
    }
    
    // Logging, population-level: alive
    std::string sLog = std::string("") + std::to_string(gWorld->getIterations()) + ",pop,alive," + std::to_string(activeCount) + "\n";
    gLogManager->write(sLog);
    gLogManager->flush();
}