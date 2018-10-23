/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef SOCIALBEHAVIORAGENTOBSERVER_H
#define SOCIALBEHAVIORAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "SocialBehavior/include/SocialBehaviorSharedData.h"

#include <iomanip>

class SocialBehaviorAgentObserver : public AgentObserver
{
	public:
		SocialBehaviorAgentObserver(RobotWorldModel *wm);
		~SocialBehaviorAgentObserver();

		void reset();
		void step();

};

#endif

