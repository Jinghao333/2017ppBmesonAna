DOANALYSISPbPb_FULL_BS=0
DOANALYSISPbPb_BINNED_PT_BS=1
DOANALYSISPbPb_BINNED_Y_BS=0
DOANALYSISPbPb_BINNED_MULT_BS=0

DOANALYSISPbPb_BINNED_PT_BS_TRK=0
DOANALYSISPbPb_BINNED_Y_BS_TRK=0
DOANALYSISPbPb_BINNED_MULT_BS_TRK=0

#Data and MC Samples
MC_Bs="~/data_mc/BsMC_nom.root"
#lstore/cms/henrique/Hadronization/data_mc/BsMC_nom.root"
Data_Bs="~/data_mc/BsData_nom.root"
#/lstore/cms/henrique/Hadronization/data_mc/BsData_nom.root"
#Data and MC Samples

#NEW NMB from https://twiki.cern.ch/twiki/pub/CMS/HINUpsilonRaa2016/Jason_MinBiasCounting_2017-02-02.pdf

CUTPbPb="Bpt>0"
cut_trk_tight="(track>1)"

mkdir -p ROOTfiles/
OutputFile_BS_FULL="ROOTfiles/yields_Bs_full"
OutputFile_BS_BINNED_Y="ROOTfiles/yields_Bs_binned_y"
OutputFile_BS_BINNED_Y_trk="ROOTfiles/yields_Bs_binned_y_trk"
OutputFile_BS_BINNED_PT="ROOTfiles/yields_Bs_binned_pt"
OutputFile_BS_BINNED_PT_trk="ROOTfiles/yields_Bs_binned_pt_trk"
OutputFile_BS_BINNED_MULT="ROOTfiles/yields_Bs_binned_Mult"
OutputFile_BS_BINNED_MULT_trk="ROOTfiles/yields_Bs_binned_Mult_trk"


#The Function to be called:
#
#void roofitB(TString tree = "ntphi", int full = 0, TString inputdata = "", TString inputmc = "", TString varExp = "", TString cut = "", TString outputfile = "", TString outplotf = "", TString jpsiFile = "", int BsBPBins = 0){
#



if [ $DOANALYSISPbPb_FULL_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('\"ntphi\"','1','\"$Data_Bs\"','\"$MC_Bs\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BS_FULL\"','\"results/Bs/Bpt\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_PT_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BS_BINNED_PT\"','\"results/Bs/Bpt\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_PT_BS_TRK -eq 1 ]; then
root -b  -q 'roofitB.C('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"Bpt\"','\"$cut_trk_tight\"','\"$OutputFile_BS_BINNED_PT_trk\"','\"results/Bs/trk_tight\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_Y_BS  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"By\"','\"$CUTPbPb\"','\"$OutputFile_BS_BINNED_Y\"','\"results/Bs/By\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_Y_BS_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"By\"','\"$cut_trk_tight\"','\"$OutputFile_BS_BINNED_Y_trk\"','\"results/Bs/By\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_MULT_BS  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"nMult\"','\"$CUTPbPb\"','\"$OutputFile_BS_BINNED_MULT\"','\"results/Bs/nMult\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_BINNED_MULT_BS_TRK  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"nMult\"','\"$cut_trk_tight\"','\"$OutputFile_BS_BINNED_MULT_trk\"','\"results/Bs/nMult\"')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi



