/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef EVOLVABILITYGRNANDODNEATAGENTOBSERVER_H
#define EVOLVABILITYGRNANDODNEATAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATSharedData.h"

#include <iomanip>

class EvolvabilityGRNandOdNEATAgentObserver : public AgentObserver
{
	public:
		EvolvabilityGRNandOdNEATAgentObserver(RobotWorldModel *wm);
		~EvolvabilityGRNandOdNEATAgentObserver();

		void reset();
		void step();

};

#endif

