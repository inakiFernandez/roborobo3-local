/*
 * MedeaConfigurationLoader.h
 */

#ifndef SOCIALBEHAVIORANNCONFIGURATIONLOADER_H
#define SOCIALBEHAVIORANNCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class SocialBehaviorANNConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		SocialBehaviorANNConfigurationLoader();
		~SocialBehaviorANNConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
