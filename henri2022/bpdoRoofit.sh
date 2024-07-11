DOANALYSISPbPb_FULL_BP=0

DOANALYSISPbPb_BINNED_PT_BP=0
DOANALYSISPbPb_BINNED_Y_BP=0
DOANALYSISPbPb_BINNED_MULTI_BP=0

DOANALYSISPbPb_BINNED_PT_BP_TRK=0
DOANALYSISPbPb_BINNED_Y_BP_TRK=0
DOANALYSISPbPb_BINNED_MULTI_BP_TRK=0

DOANALYSISPbPb_MATCHING_Bs_BINS_BP_pt=0
DOANALYSISPbPb_MATCHING_Bs_BINS_BP_pt_TRK=0
DOANALYSISPbPb_MATCHING_Bs_BINS_BP_y=0
DOANALYSISPbPb_MATCHING_Bs_BINS_BP_y_TRK=0

#Data and MC Samples
MC_BP="/lstore/cms/henrique/Hadronization/data_mc/BPMC_nom.root"
Data_BP="/lstore/cms/henrique/Hadronization/data_mc/BPData_nom.root"
#MC_BP="/data3/smcosta/data/BPMC_nom.root"
#Data_BP="/data3/smcosta/data/BPData_nom.root"

INPUTJPSI="//lstore/cms/henrique/Hadronization/data_mc/jpsinp_nom.root"
#INPUTJPSI="/data3/smcosta/data/jpsinp_nom.root"
#Data and MC Samples

#NEW NMB from https://twiki.cern.ch/twiki/pub/CMS/HINUpsilonRaa2016/Jason_MinBiasCounting_2017-02-02.pdf

CUTPbPb="Bpt>0"
cut_trk_tight="(track>1)"

mkdir -p ROOTfiles/
OutputFile_BP_FULL="ROOTfiles/yields_BP_full"
OutputFile_BP_BINNED_Y="ROOTfiles/yields_BP_binned_y"
OutputFile_BP_BINNED_Mult="ROOTfiles/yields_BP_binned_Mult"
OutputFile_BP_BINNED_PT="ROOTfiles/yields_BP_binned_pt"
OutputFile_BP_BINNED_PT_trk="ROOTfiles/yields_BP_binned_pt_trk"
OutputFile_BP_BINNED_Y_trk="ROOTfiles/yields_BP_binned_y_trk"
OutputFile_BP_BINNED_Mult_trk="ROOTfiles/yields_BP_binned_Mult_trk"
OutputFile_BP_MatchingBINS_PT="ROOTfiles/yields_BP_binned_pt_BsBPBINS"
OutputFile_BP_MatchingBINS_PT_trk="ROOTfiles/yields_BP_binned_pt_trk_BsBPBINS"
OutputFile_BP_MatchingBINS_y="ROOTfiles/yields_BP_binned_y_BsBPBINS"
OutputFile_BP_MatchingBINS_y_trk="ROOTfiles/yields_BP_binned_y_trk_BsBPBINS"


#The Function to be called:
#
#void roofitB(TString tree = "ntphi", int full = 0, TString inputdata = "", TString inputmc = "", TString varExp = "", TString cut = "", TString outputfile = "", TString outplotf = "", TString jpsiFile = "", int BsBPBins = 0){
#



if [ $DOANALYSISPbPb_FULL_BP -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','1','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BP_FULL\"','\"results/BP\"','\"$INPUTJPSI\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_PT_BP -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BP_BINNED_PT\"','\"results/BP/Bpt\"','\"$INPUTJPSI\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_PT_BP_TRK -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$cut_trk_tight\"','\"$OutputFile_BP_BINNED_PT_trk\"','\"results/BP/Bpt\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_Y_BP  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"By\"','\"$CUTPbPb\"','\"$OutputFile_BP_BINNED_Y\"','\"results/BP/By\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_Y_BP_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"By\"','\"$cut_trk_tight\"','\"$OutputFile_BP_BINNED_Y_trk\"','\"results/BP/By\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_MULTI_BP  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"nMult\"','\"$CUTPbPb\"','\"$OutputFile_BP_BINNED_Mult\"','\"results/BP/nMult\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_MULTI_BP_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"nMult\"','\"$cut_trk_tight\"','\"$OutputFile_BP_BINNED_Mult_trk\"','\"results/BP/nMult\"','\"$INPUTJPSI\"')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_MATCHING_Bs_BINS_BP_pt_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$cut_trk_tight\"','\"$OutputFile_BP_MatchingBINS_PT_trk\"','\"results/BP/Bpt\"','\"$INPUTJPSI\"','1')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_MATCHING_Bs_BINS_BP_pt  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BP_MatchingBINS_PT\"','\"results/BP/Bpt\"','\"$INPUTJPSI\"','1')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_MATCHING_Bs_BINS_BP_y_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"By\"','\"$cut_trk_tight\"','\"$OutputFile_BP_MatchingBINS_y_trk\"','\"results/BP/By\"','\"$INPUTJPSI\"','1')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_MATCHING_Bs_BINS_BP_y  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntKp\"','0','\"$Data_BP\"','\"$MC_BP\"','\"By\"','\"$CUTPbPb\"','\"$OutputFile_BP_MatchingBINS_y\"','\"results/BP/By\"','\"$INPUTJPSI\"','1')' 

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi



