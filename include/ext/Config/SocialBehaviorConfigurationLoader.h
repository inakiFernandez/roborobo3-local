/*
 * MedeaConfigurationLoader.h
 */

#ifndef SOCIALBEHAVIORCONFIGURATIONLOADER_H
#define SOCIALBEHAVIORCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class SocialBehaviorConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		SocialBehaviorConfigurationLoader();
		~SocialBehaviorConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
