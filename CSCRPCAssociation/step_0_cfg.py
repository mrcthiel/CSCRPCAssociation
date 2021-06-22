# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleMuPt100_Eta2p85_cfi --beamspot HLLHC --conditions auto:phase2_realistic_T21 --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 --datatier GEN-SIM --era Phase2C11I13M9 --eventcontent RAWSIM --fileout file:step0.root --geometry Extended2026D76 --nThreads 8 --number 10 --python_filename step_0_cfg.py --relval 9000,100 --step GEN,SIM
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9

process = cms.Process('SIM',Phase2C11I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D76_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(1000),
	output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)
# Input source
process.source = cms.Source("EmptySource")

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
	annotation = cms.untracked.string('SingleMuPt100_Eta2p85_cfi nevts:10'),
	name = cms.untracked.string('Applications'),
	version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
	SelectEvents = cms.untracked.PSet(
		SelectEvents = cms.vstring('generation_step')
	),
	compressionAlgorithm = cms.untracked.string('LZMA'),
	compressionLevel = cms.untracked.int32(1),
	dataset = cms.untracked.PSet(
		dataTier = cms.untracked.string('GEN-SIM'),
		filterName = cms.untracked.string('')
	),
	eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
	fileName = cms.untracked.string('file:step0.root'),
	outputCommands = process.RAWSIMEventContent.outputCommands,
	splitLevel = cms.untracked.int32(0)
)

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.generator = cms.EDFilter("Pythia8PtGun",
	PGunParameters = cms.PSet(
		AddAntiParticle = cms.bool(True),
		MaxEta = cms.double(2.85),
		MaxPhi = cms.double(3.14159265359),
		MaxPt = cms.double(100.01),
		MinEta = cms.double(-2.85),
		MinPhi = cms.double(-3.14159265359),
		MinPt = cms.double(99.99),
		ParticleID = cms.vint32(-13)
	),
	PythiaParameters = cms.PSet(
		parameterSets = cms.vstring()
	),
	Verbosity = cms.untracked.int32(0),
	firstRun = cms.untracked.uint32(1),
	psethack = cms.string('single mu pt 100')
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)
# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 8
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 2
process.options.eventSetup.numberOfConcurrentIOVs = 1
if hasattr(process, 'DQMStore'): process.DQMStore.assertLegacySafe=cms.untracked.bool(False)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# End of customisation functions
# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
                                                                                                                            154,1         Bot
