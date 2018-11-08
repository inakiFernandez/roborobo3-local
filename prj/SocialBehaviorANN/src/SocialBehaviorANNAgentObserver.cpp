/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include "SocialBehaviorANN/include/SocialBehaviorANNAgentObserver.h"
#include "World/World.h"
#include "Utilities/Misc.h"
#include "RoboroboMain/roborobo.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNController.h"
#include <cmath>
#include "SocialBehaviorANN/include/SocialBehaviorANNWorldObserver.h"
#include <string>

SocialBehaviorANNAgentObserver::SocialBehaviorANNAgentObserver( RobotWorldModel *wm )
{
    _wm = (RobotWorldModel*)wm;
    
}

SocialBehaviorANNAgentObserver::~SocialBehaviorANNAgentObserver()
{
    // nothing to do.
}

void SocialBehaviorANNAgentObserver::reset()
{
    // nothing to do.
}

void SocialBehaviorANNAgentObserver::step()
{
    
    // * send callback messages to objects touched or walked upon.
    
    // through distance sensors
    for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
    {
        int targetIndex = _wm->getObjectIdFromCameraSensor(i);
        
        if ( PhysicalObject::isInstanceOf(targetIndex) )   // sensor ray bumped into a physical object
        {
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            //std::cout << "[DEBUG] Robot #" << _wm->getId() << " touched " << targetIndex << "\n";
            gPhysicalObjects[targetIndex]->isTouched(_wm->getId());
        }
    }
    //TODO if(collision) _fitness = _fitness - 0.5;


    // through floor sensor
    int targetIndex = _wm->getGroundSensorValue();
    // ground sensor is upon a physical object (OR: on a place marked with
    //this physical object footprint, cf. groundsensorvalues image)
    if(PhysicalObject::isInstanceOf(targetIndex))
    {
        switch (SocialBehaviorANNSharedData::gFitness) {
        case 0:
        {
            //Floreano's increment: nothing to do
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
            break;
        }
        case 1:
        {
            //Counting items
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
            SocialBehaviorANNController* c = dynamic_cast<SocialBehaviorANNController*>(gWorld->getRobot(_wm -> _id)->getController());
            c->collectItem(); //updateFitness(1.0);
            break;
        }
        case 2:
        {
            //Recording the id of the robot for the item stepped upon
            //targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            //dynamic_cast<SocialBehaviorANNWorldObserver*>(gWorld->getWorldObserver())
            //                    ->listCollected[targetIndex].push_back(_wm->getId());
            std::cerr << "TODO collab. foraging in AgentObsever"<< std::endl;
            break;
        }
        default:
            break;
        }
    }
}
