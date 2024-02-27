import FWCore.ParameterSet.Config as cms

# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'year',
    '2022',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "add year file"
)


# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
)

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(params.inputFiles),
                        )


# Load the standard set of configuration modules
#process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')

#Modules to get L1 bits
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

#Global tag based on year
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag.globaltag = params.GlobalTagData

if params.year == '2022':
    process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v10', '')
elif params.year == '2023':
    process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_v1', '')


# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("output_scouting.root")
                               )

#L1Info = ['L1_DoubleMu_12_5', 'L1_DoubleMu_15_7', 'L1_HTT200er', 'L1_HTT255er', 'L1_HTT280er', 'L1_HTT320er', 'L1_HTT360er', 'L1_ETT2000', 'L1_HTT400er', 'L1_HTT450er', 'L1_SingleJet180', 'L1_SingleJet200', 'L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5', 'L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5', 'L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu4_SQ_OS_dR_Max1p2', 'L1_SingleEG36er2p5', 'L1_SingleLooseIsoEG28er2p1', 'L1_SingleEG8er2p5', 'L1_SingleEG10er2p5', 'L1_SingleEG15er2p5', 'L1_SingleEG26er2p5', 'L1_SingleEG28_FWD2p5', 'L1_DoubleEG4_er1p2_dR_Max0p9', 'L1_DoubleEG4p5_er1p2_dR_Max0p9', 'L1_DoubleEG5_er1p2_dR_Max0p9', 'L1_DoubleEG5p5_er1p2_dR_Max0p8', 'L1_DoubleEG7_er1p2_dR_Max0p8', 'L1_DoubleEG7p5_er1p2_dR_Max0p7', 'L1_DoubleEG_15_10_er2p5', 'L1_DoubleEG_20_10_er2p5', 'L1_DoubleEG_22_10_er2p5', 'L1_DoubleEG_25_12_er2p5', 'L1_DoubleEG_25_14_er2p5', 'L1_DoubleEG_27_14_er2p5', 'L1_DoubleEG_LooseIso22_12_er2p5', 'L1_DoubleEG_LooseIso25_12_er2p5', 'L1_TripleEG_18_17_8_er2p5', 'L1_TripleEG_18_18_12_er2p5', 'L1_TripleEG16er2p5', 'L1_DoubleEG8er2p5_HTT300er']
#Only muon seeds


#L1 triggers
L1Info = ['L1_DoubleMu_12_5', 'L1_DoubleMu_15_7', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18', 'L1_DoubleMu4_SQ_OS_dR_Max1p2', 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2']

# Make tree
process.mmtree = cms.EDAnalyzer('ScoutingAnalyzer',
                                triggerresults = cms.InputTag("TriggerResults::HLT"),
                                primaryVtx = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx", "HLT"),
                                secondaryVtx = cms.InputTag("hltScoutingMuonPacker", "displacedVtx", "HLT"),
                                pfparticles = cms.InputTag("hltScoutingPFPacker", "" ,"HLT"),
                                electrons = cms.InputTag("hltScoutingEgammaPacker", "" ,"HLT"),
                                muons = cms.InputTag("hltScoutingMuonPacker", "" ,"HLT"),
                                rho = cms.InputTag("hltScoutingPFPacker", "rho", "HLT"),
                                AlgInputTag = cms.InputTag("gtStage2Digis"),
                                l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                l1Seeds = cms.vstring(L1Info),
                                doL1 = cms.bool(False)
                            )

#Following process is to get Stage 2 Digis for L1 info
#process.p = cms.Path( process.gtStage2Digis*process.mmtree )
process.p = cms.Path( process.mmtree )