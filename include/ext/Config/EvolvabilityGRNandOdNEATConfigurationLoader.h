/*
 * MedeaConfigurationLoader.h
 */

#ifndef EVOLVABILITYGRNANDODNEATCONFIGURATIONLOADER_H
#define EVOLVABILITYGRNANDODNEATCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class EvolvabilityGRNandOdNEATConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		EvolvabilityGRNandOdNEATConfigurationLoader();
		~EvolvabilityGRNandOdNEATConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
