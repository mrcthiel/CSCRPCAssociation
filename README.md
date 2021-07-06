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
> cd CSCRPCAssociation

To run the analyzer and try to print the RE4 information:
>cmsRun cscrpc_reco.py

The information that access the RE4 information is between lines 289 and 384
