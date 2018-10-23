#if defined PRJ_SOCIALBEHAVIOR || !defined MODULAR

#include "Config/SocialBehaviorConfigurationLoader.h"

#include "SocialBehavior/include/SocialBehaviorWorldObserver.h"
#include "SocialBehavior/include/SocialBehaviorAgentObserver.h"
#include "SocialBehavior/include/SocialBehaviorController.h"

#include "WorldModels/RobotWorldModel.h"

SocialBehaviorConfigurationLoader::SocialBehaviorConfigurationLoader()
{
}

SocialBehaviorConfigurationLoader::~SocialBehaviorConfigurationLoader()
{
	//nothing to do
}

WorldObserver* SocialBehaviorConfigurationLoader::make_WorldObserver(World* wm)
{
	return new SocialBehaviorWorldObserver(wm);
}

RobotWorldModel* SocialBehaviorConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* SocialBehaviorConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new SocialBehaviorAgentObserver(wm);
}

Controller* SocialBehaviorConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new SocialBehaviorController(wm);
}

#endif
