# CSCRPCAssociation

This is a simple analyzer created with the command *mkedanlzr* and modified to print the CSC-RPC correlations on the screen. The code was adapted from the module:
https://github.com/cms-sw/cmssw/blob/master/RecoLocalMuon/RPCRecHit/src/CSCSegtoRPC.cc

To run, follow the commands.
To get the code:
> cmsrel CMSSW_10_6_0_patch2\
> cmsenv\
> cd CMSSW_10_6_0_patch2/src/\
> git clone git@github.com:mrcthiel/CSCRPCAssociation.git\
> scram b -j8\
> cd CSCRPCAssociation/CSCRPCAssociation

To run over a RAW sample and create a root file with the CSC and RPC collection:
>cmsRun L1T-PhaseIITDRSpring19GS-00007-fragment_py_GEN_SIM.py\
>cmsRun step1_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT.py\
>cmsRun step2_RAW2DIGI_L1Reco_RECO.py 

To run the analyzer and print the CSC-RPC correlations:
>cmsRun cscrpc_reco.py
