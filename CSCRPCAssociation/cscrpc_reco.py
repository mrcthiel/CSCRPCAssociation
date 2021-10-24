import FWCore.ParameterSet.Config as cms
#process = cms.Process("DUMPCSC")
process = cms.Process("demo")


#Geometry
from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")


## Global tag for 10_6 phase2 mc
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
#   fileNames = cms.untracked.vstring('file:step2.root')
   fileNames = cms.untracked.vstring('file:/eos/home-j/jgomespi/workspace/MC_ForSummerStudents/step2.root')
   ,duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string("output_rumi.root"),
	closeFileFast = cms.untracked.bool(True))



process.demo = cms.EDAnalyzer('CSCRPCAssociation',
   cscSegTag = cms.InputTag("cscSegments")
   , rpcRecHitTag = cms.InputTag("rpcRecHits")
   , cscCorrDigisTag = cms.InputTag("simCscTriggerPrimitiveDigis")
   , rpcDigisTag = cms.InputTag("simMuonRPCDigis")
   , wiresDigisTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
)

process.p = cms.Path(process.demo)
