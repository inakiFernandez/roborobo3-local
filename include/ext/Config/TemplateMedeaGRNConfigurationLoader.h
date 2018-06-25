/*
 * MedeaConfigurationLoader.h
 */

#ifndef TEMPLATEMEDEAGRNCONFIGURATIONLOADER_H
#define TEMPLATEMEDEAGRNCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class TemplateMedeaGRNConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		TemplateMedeaGRNConfigurationLoader();
		~TemplateMedeaGRNConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
