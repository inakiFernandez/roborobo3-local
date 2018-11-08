/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNWorldObserver.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNController.h"
#include "World/World.h"

SocialBehaviorANNWorldObserver::SocialBehaviorANNWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;
    
    // ==== loading project-specific properties
    
    gProperties.checkAndGetPropertyValue("gSigmaRef",&SocialBehaviorANNSharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gSigmaMin",&SocialBehaviorANNSharedData::gSigmaMin,true);
    gProperties.checkAndGetPropertyValue("gSigmaMax",&SocialBehaviorANNSharedData::gSigmaMax,true);
    
    gProperties.checkAndGetPropertyValue("gUpdateSigmaStep",&SocialBehaviorANNSharedData::gUpdateSigmaStep,true);
    gProperties.checkAndGetPropertyValue("gEvaluationTime",&SocialBehaviorANNSharedData::gEvaluationTime,true);
    gProperties.checkAndGetPropertyValue("gSynchronization",&SocialBehaviorANNSharedData::gSynchronization,true);
    
    gProperties.checkAndGetPropertyValue("gEnergyRequestOutput",&SocialBehaviorANNSharedData::gEnergyRequestOutput,false);
    
    gProperties.checkAndGetPropertyValue("gMonitorPositions",&SocialBehaviorANNSharedData::gMonitorPositions,true);
    
    
    gProperties.checkAndGetPropertyValue("gSnapshots",&SocialBehaviorANNSharedData::gSnapshots,false);
    gProperties.checkAndGetPropertyValue("gSnapshotsFrequency",&SocialBehaviorANNSharedData::gSnapshotsFrequency,false);
    
    gProperties.checkAndGetPropertyValue("gControllerType",&SocialBehaviorANNSharedData::gControllerType,true);
    
    gProperties.checkAndGetPropertyValue("gMaxNbGenomeTransmission",&SocialBehaviorANNSharedData::gMaxNbGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gLimitGenomeTransmission",&SocialBehaviorANNSharedData::gLimitGenomeTransmission,true);
    gProperties.checkAndGetPropertyValue("gSelectionMethod",&SocialBehaviorANNSharedData::gSelectionMethod,true);
    
    
    gProperties.checkAndGetPropertyValue("gLogGenome",&SocialBehaviorANNSharedData::gLogGenome,false);
    
    gProperties.checkAndGetPropertyValue("gIndividualMutationRate",&SocialBehaviorANNSharedData::gIndividualMutationRate,false);

    gProperties.checkAndGetPropertyValue("gMutationOperator",&SocialBehaviorANNSharedData::gMutationOperator,false);
    gProperties.checkAndGetPropertyValue("gFitness",&SocialBehaviorANNSharedData::gFitness,true);

    gProperties.checkAndGetPropertyValue("gSelPressure",&SocialBehaviorANNSharedData::gSelPressure,true);
    gProperties.checkAndGetPropertyValue("gPopSize",&SocialBehaviorANNSharedData::gPopSize,true);


    gProperties.checkAndGetPropertyValue("gMutateProb",&SocialBehaviorANNSharedData::gMutateProb,true);
    gProperties.checkAndGetPropertyValue("gCrossoverProb",&SocialBehaviorANNSharedData::gCrossoverProb,true);
    gProperties.checkAndGetPropertyValue("gNbRegulatory",&SocialBehaviorANNSharedData::gNbRegulatory,true);


    gProperties.checkAndGetPropertyValue("gCommunicationRange",&SocialBehaviorANNSharedData::gCommunicationRange,true);
    gProperties.checkAndGetPropertyValue("gCommunicationOnRadius",&SocialBehaviorANNSharedData::gCommunicationOnRadius,true);

    gProperties.checkAndGetPropertyValue("gMatingOperator",&SocialBehaviorANNSharedData::gMatingOperator,true);
    gProperties.checkAndGetPropertyValue("gBroadcastTime",&SocialBehaviorANNSharedData::gBroadcastTime,true);    

    gProperties.checkAndGetPropertyValue("mutate_link_weights_prob",&Helper::mutateLinkWeightsProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_node_prob",&Helper::mutateAddNodeProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_link_prob",&Helper::mutateAddLinkProb,true);
    gProperties.checkAndGetPropertyValue("recur_only_prob",&Helper::recurOnlyProb,true);
    gProperties.checkAndGetPropertyValue("newstructure_tries",&Helper::newStructureTries,true);


    //TODO read other properties (cf shareddata)
    if ( !gRadioNetwork)
    {
        std::cout << "Error : gRadioNetwork must be true." << std::endl;
        exit(-1);
    }
    
    // * iteration and generation counters
    _generationItCount = -1;
    _generationCount = -1;
    int nbIn = -1;
    switch (SocialBehaviorANNSharedData::gFitness) {
    case 0:
        nbIn = 8;
        break;
    case 1:
        nbIn = 16;
        break;
    case 2:
        nbIn = 16;
        break;
    default:
        std::cerr << "[ERROR] Wrong or non implemented task" << std::endl;
        exit(-1);
        break;
    }
SocialBehaviorANNSharedData::initInputsBehavior(SocialBehaviorANNSharedData::gNbInputsBehavior, nbIn);
}

SocialBehaviorANNWorldObserver::~SocialBehaviorANNWorldObserver()
{
    // nothing to do.
}

void SocialBehaviorANNWorldObserver::reset()
{
}

void SocialBehaviorANNWorldObserver::step()
{
    _generationItCount++;
    
    if( _generationItCount == SocialBehaviorANNSharedData::gEvaluationTime+1 ) // switch to next generation.
    {
        // update iterations and generations counters
        _generationItCount = 0;
        _generationCount++;
    }
    //Local broadcast with distance: compute inter-robot distance matrix
    if(SocialBehaviorANNSharedData::gCommunicationOnRadius)
    {
    _robotDistances = std::vector<std::vector<double> >(gNumberOfRobots, std::vector<double>(gNumberOfRobots,-1));
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {
            SocialBehaviorANNController* c = (dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(i)->getController()));
            double x = c->getWorldModel()->getXReal();
            double y = c->getWorldModel()->getYReal();

            for ( int j = i ; j != gNumberOfRobots ; j++ )
            {
                SocialBehaviorANNController* c2 = (dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(j)->getController()));
                double x2 = c2->getWorldModel()->getXReal();
                double y2 = c2->getWorldModel()->getYReal();
                _robotDistances[i][j] = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)); //.push_back(sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)));
                _robotDistances[j][i] = _robotDistances[i][j];
                //std::cout << _robotDistances[i][j] << "      ";
            }
            //std::cout << std::endl;

        }
        //std::cout << SocialBehaviorANNSharedData::gCommunicationRange << std::endl << std::endl << std::endl;
    }
    updateMonitoring();
    
    updateEnvironment();
    
}


void SocialBehaviorANNWorldObserver::updateEnvironment()
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

void SocialBehaviorANNWorldObserver::updateMonitoring()
{
    // * Log at end of each generation

    if( //gWorld->getIterations() % SocialBehaviorANNSharedData::gEvaluationTime == 1 ||
            gWorld->getIterations() % SocialBehaviorANNSharedData::gEvaluationTime == SocialBehaviorANNSharedData::gEvaluationTime-1 ) // beginning(+1) *and* end of generation. ("==1" is required to monitor the outcome of the first iteration)
    {
        monitorPopulation();
    }
    //Simulate centralized algorithm by maintaining panmictic communication
    if(SocialBehaviorANNSharedData::gIsCentralized)
    {
        //Send around the end of the evaluation
        if((gWorld->getIterations() + 10) % SocialBehaviorANNSharedData::gEvaluationTime == 0)
        //if( _lifeIterationCount == SocialBehaviorANNSharedData::gEvaluationTime-10)
        {
            for (int i = 0 ; i != gNumberOfRobots ; i++ )
            {
                SocialBehaviorANNController* c1 = (dynamic_cast<SocialBehaviorANNController*>
                                          (gWorld->getRobot(i)->getController()));

                for (int j = 0 ; j != gNumberOfRobots ; j++ )
                {
                    SocialBehaviorANNController* c2 = (dynamic_cast<SocialBehaviorANNController*>
                                              (gWorld->getRobot(j)->getController()));
                    c2->storeGenomeHelper(c1->getCurrentGenome());//,c1->getCurrentId(), c1->getFitness(), c1->getBehavior());
                }
            }
        }
    }
    // * Every N generations, take a video (duration: one generation time)

    if ( SocialBehaviorANNSharedData::gSnapshots )
    {
        if ( ( gWorld->getIterations() ) % ( SocialBehaviorANNSharedData::gEvaluationTime * SocialBehaviorANNSharedData::gSnapshotsFrequency ) == 0 )
        {
            if ( gVerbose )
                std::cout << "[START] Video recording: generation #" << (gWorld->getIterations() / SocialBehaviorANNSharedData::gEvaluationTime ) << ".\n";
            gTrajectoryMonitorMode = 0;
            initTrajectoriesMonitor();
        }
        else
            if ( ( gWorld->getIterations() ) % ( SocialBehaviorANNSharedData::gEvaluationTime * SocialBehaviorANNSharedData::gSnapshotsFrequency ) == SocialBehaviorANNSharedData::gEvaluationTime - 1 )
            {
                if ( gVerbose )
                    std::cout << "[STOP]  Video recording: generation #" << (gWorld->getIterations() / SocialBehaviorANNSharedData::gEvaluationTime ) << ".\n";
                saveTrajectoryImage();
            }
    }    
}

void SocialBehaviorANNWorldObserver::monitorPopulation( bool localVerbose )
{
    // * monitoring: count number of active agents.
    
    int activeCount = 0;
    for ( int i = 0 ; i != gNumberOfRobots ; i++ )
    {
        if ( (dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(i)->getController()))->getWorldModel()->isAlive() == true )
            activeCount++;
    }
    
    if ( gVerbose && localVerbose )
    {
        std::cout << "[gen:" << (gWorld->getIterations()/SocialBehaviorANNSharedData::gEvaluationTime) << ";it:" << gWorld->getIterations() << ";pop:" << activeCount << "]\n";
    }
    // Logging here
    double sumFitness = 0.0, sumAvgLocalPopFitness = 0.0;
    double sumItems = 0.0;
    double sumCollisions = 0.0;
    double gatheredGenomes = 0.0;
    double distance = 0.0;
    double nbUnits = 0.0;

    for ( int i = 0 ; i != gNumberOfRobots ; i++ )
    {

        SocialBehaviorANNController* c =(dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(i)->getController()));
         sumFitness += c->getFitness();
         sumAvgLocalPopFitness += c-> getAvgPopFitness();
         sumItems += c->getCollectedItems();
         sumCollisions += c->getNbCollisions() - c->getCollectedItems();
         gatheredGenomes += c->getGenomesList().size();
         distance += c->getDSumTravelled();
         nbUnits += c->getNbUnits();
    }
    switch (SocialBehaviorANNSharedData::gFitness)
    {
        case 1:
        {
            sumFitness *= SocialBehaviorANNSharedData::gEvaluationTime;
            break;
        }
    }
    std::cout << gWorld->getIterations() << " ";
    //<< (sumFitness  / gNumberOfRobots) / SocialBehaviorANNSharedData::gEvaluationTime

    // divided by two because each item gives 1 fitness point to both agents
    //std::cout << sumFitness / 2 << " " << sumAvgLocalPopFitness / 2 << std::endl;

    std::cout << sumFitness/(double) gNumberOfRobots << " ";
                 //" " << sumAvgLocalPopFitness <<
    if(SocialBehaviorANNSharedData::gDoMeasureDiv)
    {
        std::cout << computeGlobalDiversity();
    }
                        //" " << computeInterRobotDiversity();
    std::cout << " "
    << sumItems / gNumberOfRobots
    << " "
    << sumCollisions / gNumberOfRobots
    << " "
    << gatheredGenomes / gNumberOfRobots

                 //  " " << fullLaps <<
              //   " " << towardDoor0 <<
              //   " " << towardDoor1
    << "       "
    << distance / gNumberOfRobots
    << " "
    << nbUnits / gNumberOfRobots //nb neurones ou proteines (répartition des tailles de réseau par quartiles de fitness)
    << " "; // nbOffspring + fitness pour mesurer le kendall-tau-b pression de sélection

   /* if(SocialBehaviorANNSharedData::gDoMeasureDiv)
   {

        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

            SocialBehaviorANNController* c =(dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(i)->getController()));
             //std::cout << c-> getDoorPassages() << " ";
             std::cout << " " << c->computeIntraRobotDiversity() << " ";
        }
    }*/
    std::cout << std::endl;

    //std::cout //<< "[gen:"
    //        << (gWorld->getIterations()/SocialBehaviorANNSharedData::gEvaluationTime)
            //<< ", fitness: "
           // << " "
   //         << sumFitness / gNumberOfRobots
            //<< ", gatheredGenomes: "


    // Logging, population-level: alive
    //std::string sLog = std::string("") + std::to_string(gWorld->getIterations()) + ",pop,alive," + std::to_string(activeCount) + "\n";
    //gLogManager->write(sLog);
    //gLogManager->flush();
}
double SocialBehaviorANNWorldObserver::computeGlobalDiversity()
{
    double result = 0.0;
    std::vector<double> pairwiseDistances;
    if(gNumberOfRobots == 1)
        return 0.0;
    for (int i = 0 ; i < gNumberOfRobots ; i++ )
    {
        SocialBehaviorANNController* c1 = (dynamic_cast<SocialBehaviorANNController*>
                                  (gWorld->getRobot(i)->getController()));

        for (int j = i + 1 ; j < gNumberOfRobots ; j++ )
        {
            SocialBehaviorANNController* c2 = (dynamic_cast<SocialBehaviorANNController*>
                                      (gWorld->getRobot(j)->getController()));
            pairwiseDistances.push_back(c1->computeBehavDistance(c1->getBehavior(),c2->getBehavior()));
        }
    }
    for(auto d : pairwiseDistances)
    {
        result += d;
    }
    result = result/pairwiseDistances.size();
    //Average of pairwise distances between active behaviors
    return result;

}
