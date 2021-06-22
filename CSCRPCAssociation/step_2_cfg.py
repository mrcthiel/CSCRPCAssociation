# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --filein file:step1.root --fileout file:step2.root --eventcontent AODSIM --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 --datatier AODSIM --step RAW2DIGI,L1Reco,RECO --conditions auto:phase2_realistic_T21 --geometry Extended2026D76 --era Phase2C11I13M9 --python_filename step_2_cfg.py --nThreads 8
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9

process = cms.Process('RECO',Phase2C11I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1),	
	output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:step1.root'),
	secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
	FailPath = cms.untracked.vstring(),
	IgnoreCompletely = cms.untracked.vstring(),
	Rethrow = cms.untracked.vstring(),
	SkipEvent = cms.untracked.vstring(),
	allowUnscheduled = cms.obsolete.untracked.bool,
	canDeleteEarly = cms.untracked.vstring(),
	deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
	emptyRunLumiMode = cms.obsolete.untracked.string,	
	eventSetup = cms.untracked.PSet(
		forceNumberOfConcurrentIOVs = cms.untracked.PSet(
			allowAnyLabel_=cms.required.untracked.uint32
		),
		numberOfConcurrentIOVs = cms.untracked.uint32(1)
	),	
	fileMode = cms.untracked.string('FULLMERGE'),
	forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
	makeTriggerResults = cms.obsolete.untracked.bool,
	numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),	
	numberOfConcurrentRuns = cms.untracked.uint32(1),
	numberOfStreams = cms.untracked.uint32(0),
	numberOfThreads = cms.untracked.uint32(8),
	printDependencies = cms.untracked.bool(False),
	sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,	
	throwIfIllegalParameter = cms.untracked.bool(True),
	wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
	annotation = cms.untracked.string('step2 nevts:1'),
	name = cms.untracked.string('Applications'),
	version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
	compressionAlgorithm = cms.untracked.string('LZMA'),
	compressionLevel = cms.untracked.int32(4),
	dataset = cms.untracked.PSet(
		dataTier = cms.untracked.string('AODSIM'),
		filterName = cms.untracked.string('')
	),
	eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
	fileName = cms.untracked.string('file:step2.root'),
	outputCommands = process.AODSIMEventContent.outputCommands
)
# Joao
process.AODSIMoutput.outputCommands = cms.untracked.vstring('drop *')
process.AODSIMoutput.outputCommands = cms.untracked.vstring('keep *_*simMuonRPCDigis*_*_*','keep *_*csc*_*_*','keep *_*rpcRecHits*_*_*','keep *_*simCscTriggerPrimitiveDigis*_*_*')
# Additional output definition
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 8
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 2
process.options.eventSetup.numberOfConcurrentIOVs = 1
if hasattr(process, 'DQMStore'): process.DQMStore.assertLegacySafe=cms.untracked.bool(False)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# End of customisation functions


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
