#if defined PRJ_TEMPLATEMEDEAGRN || !defined MODULAR

#include "Config/TemplateMedeaGRNConfigurationLoader.h"

#include "TemplateMedeaGRN/include/TemplateMedeaGRNWorldObserver.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNAgentObserver.h"
#include "TemplateMedeaGRN/include/TemplateMedeaGRNController.h"

#include "WorldModels/RobotWorldModel.h"

TemplateMedeaGRNConfigurationLoader::TemplateMedeaGRNConfigurationLoader()
{
}

TemplateMedeaGRNConfigurationLoader::~TemplateMedeaGRNConfigurationLoader()
{
	//nothing to do
}

WorldObserver* TemplateMedeaGRNConfigurationLoader::make_WorldObserver(World* wm)
{
	return new TemplateMedeaGRNWorldObserver(wm);
}

RobotWorldModel* TemplateMedeaGRNConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* TemplateMedeaGRNConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new TemplateMedeaGRNAgentObserver(wm);
}

Controller* TemplateMedeaGRNConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new TemplateMedeaGRNController(wm);
}

#endif
