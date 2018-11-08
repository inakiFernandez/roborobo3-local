#if defined PRJ_SOCIALBEHAVIORANN || !defined MODULAR

#include "Config/SocialBehaviorANNConfigurationLoader.h"

#include "SocialBehaviorANN/include/SocialBehaviorANNWorldObserver.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNAgentObserver.h"
#include "SocialBehaviorANN/include/SocialBehaviorANNController.h"

#include "WorldModels/RobotWorldModel.h"

SocialBehaviorANNConfigurationLoader::SocialBehaviorANNConfigurationLoader()
{
}

SocialBehaviorANNConfigurationLoader::~SocialBehaviorANNConfigurationLoader()
{
	//nothing to do
}

WorldObserver* SocialBehaviorANNConfigurationLoader::make_WorldObserver(World* wm)
{
	return new SocialBehaviorANNWorldObserver(wm);
}

RobotWorldModel* SocialBehaviorANNConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* SocialBehaviorANNConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new SocialBehaviorANNAgentObserver(wm);
}

Controller* SocialBehaviorANNConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new SocialBehaviorANNController(wm);
}

#endif
