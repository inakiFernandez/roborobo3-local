#if defined PRJ_EVOLVABILITYGRNANDODNEAT || !defined MODULAR

#include "Config/EvolvabilityGRNandOdNEATConfigurationLoader.h"

#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATWorldObserver.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATAgentObserver.h"
#include "EvolvabilityGRNandOdNEAT/include/EvolvabilityGRNandOdNEATController.h"

#include "WorldModels/RobotWorldModel.h"

EvolvabilityGRNandOdNEATConfigurationLoader::EvolvabilityGRNandOdNEATConfigurationLoader()
{
}

EvolvabilityGRNandOdNEATConfigurationLoader::~EvolvabilityGRNandOdNEATConfigurationLoader()
{
	//nothing to do
}

WorldObserver* EvolvabilityGRNandOdNEATConfigurationLoader::make_WorldObserver(World* wm)
{
	return new EvolvabilityGRNandOdNEATWorldObserver(wm);
}

RobotWorldModel* EvolvabilityGRNandOdNEATConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* EvolvabilityGRNandOdNEATConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new EvolvabilityGRNandOdNEATAgentObserver(wm);
}

Controller* EvolvabilityGRNandOdNEATConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new EvolvabilityGRNandOdNEATController(wm);
}

#endif
