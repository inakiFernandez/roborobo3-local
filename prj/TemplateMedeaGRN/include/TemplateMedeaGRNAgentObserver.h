/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef TEMPLATEMEDEAGRNAGENTOBSERVER_H
#define TEMPLATEMEDEAGRNAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNSharedData.h"

#include <iomanip>

class TemplateMedeaGRNAgentObserver : public AgentObserver
{
	public:
		TemplateMedeaGRNAgentObserver(RobotWorldModel *wm);
		~TemplateMedeaGRNAgentObserver();

		void reset();
		void step();

};

#endif

