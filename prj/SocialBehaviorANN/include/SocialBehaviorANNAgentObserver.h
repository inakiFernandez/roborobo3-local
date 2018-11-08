/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef SOCIALBEHAVIORANNAGENTOBSERVER_H
#define SOCIALBEHAVIORANNAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNSharedData.h"

#include <iomanip>

class SocialBehaviorANNAgentObserver : public AgentObserver
{
	public:
		SocialBehaviorANNAgentObserver(RobotWorldModel *wm);
		~SocialBehaviorANNAgentObserver();

		void reset();
		void step();

};

#endif

