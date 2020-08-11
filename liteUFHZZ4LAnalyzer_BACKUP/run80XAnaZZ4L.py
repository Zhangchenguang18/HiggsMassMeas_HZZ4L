import os, time

#dirMC = '/cms/data/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Run1Fid_20160222'
#dirMC = '/cms/data/store/user/t2/users/dsperka/dsperka/rootfiles_MC80X_20160716_MuCalib/'
#dirMC = '/cms/data/store/user/t2/users/mhl/'
#dirMC = '/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/HiggsMass_HZZ4L/packages/doMeasurement/Fit_PereventMerr_tmp/'
#dirMC = '/cms/data/store/user/hmei/UFHZZAnalysisRun2/myTask_2015MC_ggH_m125_v2/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/crab_GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2/161006_100207/0000/'

#2016MC
#dirMC = '/cms/data/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_4lskim_M17_Jan26/'
#dirMC = '/cms/data/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_4lskim_M17_Feb11/'
#dirMC = '/cms/data/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_4lskim_M17_Feb21/'
#dirMC = '/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/HiggsMass_HZZ4L/packages/liteUFHZZ4LAnalyzer/skimTrees/'
#dirMC = '/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/CMSSW_8_0_26_patch1/src/'
#dirMC = '/cms/data/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_hzz4lskim_M17_Feb21/'
#dirMC = '/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/test/CMSSW_8_0_26_patch1/src/'
#dirMC = '/raid/raid8/mhl/HZZ4L_Run2_post2016ICHEP/HiggsMass_HZZ4L/packages/inputRoots_2016MC/'

#dirMC = '/raid/raid9/chenguan/Mass/input/'
#dirMC = '/raid/raid9/qguo/Pro_NT_2016_STXS1P1/DataMC_2016_new1p1/'


#2016MC from lucien
#dirMC = '/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_4l_Feb21/'

#Run2MC from ferrico
#dirMC = '/store/user/t2/users/ferrico/Full_RunII/ggH/'
#dirMC = '/store/user/t2/users/ferrico/Full_RunII/ZH/'
#Run2MC made by myself
dirMC = '/raid/raid9/chenguan/input/AfterSync_200424/'
#dirMC = '/store/user/t2/users/ferrico/Qianying/2017/DATA/'
samplesMC  = [
'2017GGH_125'
#'2017_noDuplicates'
#'slimmed_2017qqZZ'
#'2017ggZZ_4mu'
#'ZH_HToZZ_4LFilter_M125_2018'	
#'2017GGH_125_vtx'
#'GluGluHToZZTo4L_M125_2017'#from ferrico
#'GluGluHToZZTo4L_M120_2017'
#'GluGluHToZZTo4L_M124_2017'
#'GluGluHToZZTo4L_M126_2017'
#'GluGluHToZZTo4L_M130_2017'
#NTuple made by my self
#'GGH_M120_2016'
#'GGH_M124_2016'
#'GGH_M125_2017'
#'GGH_M126_2016'
#'GGH_M130_2016'

#these samples are in /cms/data/store/user/t2/users/klo/Higgs/HZZ4l/NTuple/Run2/MC80X_M17_4l_Feb21/
#'ZZTo4L_13TeV_powheg_pythia8'
#'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8'
#'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8'
#'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8'
#'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8'
#'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8'
#'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8'
#'GluGluHToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8'
#'GluGluHToZZTo4L_M124_13TeV_powheg2_JHUgenV6_pythia8'
#'GluGluHToZZTo4L_M126_13TeV_powheg2_JHUgenV6_pythia8'
#'GluGluHToZZTo4L_M130_13TeV_powheg2_JHUgenV6_pythia8'
#'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8'
#'WH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8'
#'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8'
#'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8'
#'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8'
#'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8'




#76x
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2_20170627_v1'
#'ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2_20170627_v1'
#'Sync_80X_Moriond'
#'Data_Run2016-03Feb2017_hzz4l'
#'Sync_80X_Moriond'
#'skimTreesqqZZ_lowMass'
#'skimTreesggZZ_4mu_lowMass'
#'skimTreesggZZ_4e_lowMass'
#'skimTreesggZZ_2e2mu_lowMass'
#'ZZTo4L_13TeV_powheg_pythia8_ext1_RunIISummer16MiniAODv2'
#'ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2_tchan_ext'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2'
#'GluGluHToZZTo4L_M124_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2'
#'GluGluHToZZTo4L_M126_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2'
#'GluGluHToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2'
#'GluGluHToZZTo4L_M130_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2'
#'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2'
#'ZZTo4L_13TeV_powheg_pythia8_RunIISummer16MiniAODv2'
#'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2'
#'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2'
#'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISummer16MiniAODv2'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1'
#'GluGluHToZZTo4L_M124_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3'
#'GluGluHToZZTo4L_M126_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1'
#'GluGluHToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1'
#'GluGluHToZZTo4L_M130_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1'
##80x
#'ggZZ'
#'qqZZ'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2'
#'ggH_2015MC_mH125'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2_3'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1'
#'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2'
#'ZZTo4L_13TeV_powheg_pythia8_RunIISpring16MiniAODv2'
#'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv1'
#'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2'
#'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2'
]

njobs = 10
for job in range(1,njobs+1):
#  if ( job <= 5 and job > 10): continue
#  if ( job <= 10 or job > 20): continue
#  if ( job <= 20 or job > 30): continue
#  if ( job <= 30 or job > 40): continue
  for sample in samplesMC:
#    cmd = './ZZ4L_Ana.exe '+dirMC+'/'+sample+' Ntuples/'+sample+' 0 '+str(job)+' '+str(njobs)+' &'
    cmd = 'nohup ./ZZ4L_Ana.exe '+dirMC+'/'+sample+' Ntuples/'+sample+' 0 '+str(job)+' '+str(njobs)+' >& Dump/'+sample+'_'+str(job)+'.log &'
#    cmd = './ZZ4L_Ana.exe '+dirMC+'/'+sample+' Ntuples/'+sample+' 0 '+str(job)+' '+str(njobs)#+' >& Dump/'+sample+'_'+str(job)+'.log &'

    print cmd
    os.system(cmd)

