from ROOT import *
from array import array
import os

SamplesData = [
'Data_2016.root',
]

SamplesMC = [
'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1',
'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2',
'TTTo2L2Nu_13TeV-powheg_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1',
'WWTo2L2Nu_13TeV-powheg_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'ZZTo2L2Nu_13TeV_powheg_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'GluGluHToZZTo2L2Q_M750_13TeV_powheg2_JHUgenV698_pythia8_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1',
'VBF_HToZZTo2L2Q_M750_13TeV_powheg2_JHUgenV698_pythia8_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1'

]

###################################################### 
RootFile = {} 
Tree = {} 
nEvents = {} 
sumw = {}

def LoadData():
    for i in range(0,len(SamplesMC)):

        print i,SamplesMC[i]
        sample = SamplesMC[i].rstrip('.root')
        
        RootFile[sample] = TFile('ZjetsNtuples/'+sample+'.root',"READ")
        Tree[sample]  = RootFile[sample].Get("passedEvents")
        
        h_nevents = RootFile[sample].Get("Ana/nEvents")
        h_sumw = RootFile[sample].Get("Ana/sumWeights")
        
        if (h_nevents): nEvents[sample] = h_nevents.Integral()
        else: nEvents[sample] = 0.

        if (h_sumw): sumw[sample] = h_sumw.Integral()
        else: sumw[sample] = 0.

        if (not Tree[sample]): print sample+' has no passedEvents tree'
        else:
            print sample,"nevents",nEvents[sample],"sumw",sumw[sample]

    for i in range(0,len(SamplesData)):
        
        sample = SamplesData[i].rstrip('.root')
        
        RootFile[sample] = TFile('ZjetsNtuples/'+sample+'.root',"READ")
        Tree[sample]  = RootFile[sample].Get("passedEvents")
        
        
