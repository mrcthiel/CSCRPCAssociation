# CSCRPCAssociation

This is a simple analyzer created with the command *mkedanlzr* and modified to print the CSC-RPC correlations in the screem. The code was adapted the module:
https://github.com/cms-sw/cmssw/blob/master/RecoLocalMuon/RPCRecHit/src/CSCSegtoRPC.cc

To run, follow the commands.
To get the code:
> cmsrel CMSSW_10_6_0_patch2\
> cmsenv\
> cd CMSSW_10_6_0_patch2/src/\
> git clone git@github.com:mrcthiel/CSCRPCAssociation.git\
> scram b -j8\
> cd CSCRPCAssociation/CSCRPCAssociation\

To run:
>cmsRun run_CSC_and_RPC_rechits_and_segments_cfg.py\

To run:
>cmsRun cscrpc_reco.py\



