import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_tools import *
from os import path

mvaTag = "Summer16IdIso"

weightFileDir = "RecoEgamma/ElectronIdentification/data/MVAWeightFiles/Summer_16_ID_ISO"

mvaWeightFiles = cms.vstring(
     path.join(weightFileDir, "EB1_5.weights.xml.gz"), # EB1_5
     path.join(weightFileDir, "EB2_5.weights.xml.gz"), # EB2_5
     path.join(weightFileDir, "EE_5.weights.xml.gz"), # EE_5
     path.join(weightFileDir, "EB1_10.weights.xml.gz"), # EB1_10
     path.join(weightFileDir, "EB2_10.weights.xml.gz"), # EB2_10
     path.join(weightFileDir, "EE_10.weights.xml.gz"), # EE_10
     )

categoryCuts = cms.vstring(
     "pt < 10. & abs(superCluster.eta) < 0.800", # EB1_5
     "pt < 10. & abs(superCluster.eta) >= 0.800 & abs(superCluster.eta) < 1.479", # EB2_5
     "pt < 10. & abs(superCluster.eta) >= 1.479", # EE_5
     "pt >= 10. & abs(superCluster.eta) < 0.800", # EB1_10
     "pt >= 10. & abs(superCluster.eta) >= 0.800 & abs(superCluster.eta) < 1.479", # EB2_10
     "pt >= 10. & abs(superCluster.eta) >= 1.479", # EE_10
     )

mvaEleID_Summer16_ID_ISO_HZZ_container = EleMVARaw_WP(
    idName = "mvaEleID-Summer16-ID-ISO-HZZ", mvaTag = mvaTag,
    cutCategory0 = "1.83536653843", # EB1_5
    cutCategory1 = "1.79287061619", # EB2_5
    cutCategory2 = "1.72721183566", # EE_5
    cutCategory3 = "0.397999302311", # EB1_10
    cutCategory4 = "0.375408289311", # EB2_10
    cutCategory5 = "-0.654236429469", # EE_10
    )


mvaEleID_Summer16_ID_ISO_producer_config = cms.PSet(
    mvaName             = cms.string(mvaClassName),
    mvaTag              = cms.string(mvaTag),
    nCategories         = cms.int32(6),
    categoryCuts        = categoryCuts,
    weightFileNames     = mvaWeightFiles,
    variableDefinition  = cms.string(mvaVariablesFile)
    )

mvaEleID_Summer16_ID_ISO_HZZ = configureVIDMVAEleID( mvaEleID_Summer16_ID_ISO_HZZ_container )

mvaEleID_Summer16_ID_ISO_HZZ.isPOGApproved = cms.untracked.bool(True)
