# CSCRPCAssociation

To run, follow the commands.
To get the code:
> cmsrel CMSSW_12_0_0_pre1\
> cd CMSSW_12_0_0_pre1/src\
> cmsenv\
> git clone git@github.com:jgomespi/CSCRPCAssociation.git\
> cd CSCRPCAssociation\
> git checkout new_branch\
> scram b -j8\
> cd CSCRPCAssociation\

To run over a RAW sample and create a root file with the CSC and RPC collection:
>cmsRun L1T-PhaseIITDRSpring19GS-00007-fragment_py_GEN_SIM.py\
>cmsRun step1_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT.py\
>cmsRun step2_RAW2DIGI_L1Reco_RECO.py 

To run the analyzer and print the CSC-RPC correlations:
>cmsRun cscrpc_reco.py
