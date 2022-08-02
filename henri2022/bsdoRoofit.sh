DOANALYSISPbPb_ROOFIT_BINNED_PT_BS=1
DOANALYSISPbPb_ROOFIT_BINNED_PT_BS_TRK=0
DOANALYSISPbPb_ROOFIT_BINNED_MULT_BS=0
DOANALYSISPbPb_ROOFIT_BINNED_Y_BS=0
DOANALYSISPbPb_ROOFIT_FULL_BS=0


#INPUTDATAPbPbCANDWISE_BS="../../SkimmedSamples/BsData.root"
#INPUTMCPbPbCANDWISE_BS="../../SkimmedSamples/OfficialMC/BsMC.root"

INPUTDATAPbPbCANDWISE_BS="~/data/BsData.root"
INPUTMCPbPbCANDWISE_BS="~/data/BsMC.root"
#INPUTDATAPbPbCANDWISE_BS="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/Bs/NewCut/BsData.root"
#INPUTMCPbPbCANDWISE_BS="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/Bs/NewCut/BsMC.root"

LABEL=""
#LUMIPbPb=13.1983052423 #paper 20170227
LUMIPbPb=56.564165324
#NEW NMB from https://twiki.cern.ch/twiki/pub/CMS/HINUpsilonRaa2016/Jason_MinBiasCounting_2017-02-02.pdf
#2339717742. 2339717742*5.607/1e9 = 13.1187973794
#13.1187973794/0.99*0.996 = 13.1983052423
NMBEVT=2329685794.627413
ISMCPbPb=0
ISDOWEIGHTPbPb=1
SELGENPbPb="TMath::Abs(Gy)<2.4&&TMath::Abs(GpdgId)==531&&GisSignal>0"
SELGENPbPbACCPbPb="TMath::Abs(Gy)<2.4&&abs(GpdgId)==531&&GisSignal>0&&((TMath::Abs(Gmu1eta)<1.2&&Gmu1pt>3.5)||(TMath::Abs(Gmu1eta)>1.2&&TMath::Abs(Gmu1eta)<2.1&&Gmu1pt>(5.77-1.8*TMath::Abs(Gmu1eta)))||(TMath::Abs(Gmu1eta)>2.1&&TMath::Abs(Gmu1eta)<2.4&&Gmu1pt>1.8))&&((TMath::Abs(Gmu2eta)<1.2&&Gmu2pt>3.5)||(TMath::Abs(Gmu2eta)>1.2&&TMath::Abs(Gmu2eta)<2.1&&Gmu2pt>(5.77-1.8*TMath::Abs(Gmu2eta)))||(TMath::Abs(Gmu2eta)>2.1&&TMath::Abs(Gmu2eta)<2.4&&Gmu2pt>1.8))&&Gtk1pt>0.&&Gtk2pt>0.&&TMath::Abs(Gtk1eta)<2.4&&TMath::Abs(Gtk2eta)<2.4"

#BASECUTPbPb="(hiBin<181)&&Btrk1Pt>0.9&&Btrk2Pt>0.9&&Bchi2cl>0.05&&BsvpvDistance/BsvpvDisErr>2&&Bpt>5&&abs(Btrk1Eta-0.0)<2.4&&abs(Btrk2Eta-0.0)<2.4&&(TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.77-1.8*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.8))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.77-1.8*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.8))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>1.&&Btrk2Pt>1.&&abs(Btktkmass-1.019455)<0.015)&&(abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter)&&(Btrk1PixelHit+Btrk1StripHit>10)&&(Btrk2PixelHit+Btrk2StripHit>10)&&(Btrk1PtErr/Btrk1Pt<0.1)&&(Btrk2PtErr/Btrk2Pt<0.1)&&Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer)<0.18&&Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer)<0
BASECUTPbPb="(hiBin<181)&&Btrk1Pt>1.0&&Btrk2Pt>1.0&&Bchi2cl>0.05&&BsvpvDistance/BsvpvDisErr>2.2&&Bpt>5&&abs(Btrk1Eta-0.0)<2.4&&abs(Btrk2Eta-0.0)<2.4&&(TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.77-1.8*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.8))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.77-1.8*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.8))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>1.&&Btrk2Pt>1.&&abs(Btktkmass-1.019455)<0.015)&&(abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter)&&(Btrk1PixelHit+Btrk1StripHit>10)&&(Btrk2PixelHit+Btrk2StripHit>10)&&(Btrk1PtErr/Btrk1Pt<0.1)&&(Btrk2PtErr/Btrk2Pt<0.1)&&Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer)<0.18&&Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer)<0.18"

CUTPbPb=${BASECUTPbPb}"&&((Bpt>5&&Bpt<10&&BDT_pt_5_10>0.17)||(Bpt>10&&Bpt<15&&BDT_pt_10_15>0.17)||(Bpt>15&&Bpt<20&&BDT_pt_15_20>0.26)||(Bpt>20&&Bpt<50&&BDT_pt_20_50>0.25))"

CUTPbPb=${CUTPbPb}"&&abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter"


CUTPbPb="1"

#TRGPbPb="(HLT_HIL1DoubleMu0_v1||HLT_HIL1DoubleMu0_part1_v1||HLT_HIL1DoubleMu0_part2_v1||HLT_HIL1DoubleMu0_part3_v1)"
#TRGPbPbMC="(HLT_HIL1DoubleMu0_v1||HLT_HIL1DoubleMu0_part1_v1||HLT_HIL1DoubleMu0_part2_v1||HLT_HIL1DoubleMu0_part3_v1)"
TRGPbPb="(Bpt>0)"
TRGPbPbMC="(Bpt>0)"

echo "TRGPbPb="$TRGPbPb

mkdir -p ROOTfiles/
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_FULL="ROOTfiles/yields_Bs_full_${CENTPbPbMIN}_${CENTPbPbMAX}.root"

OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_1ST="ROOTfiles/yields_Bs_binned_cent_1st_pt.root"
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_2ND="ROOTfiles/yields_Bs_binned_cent_2nd_pt.root"
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_1ST_Y="ROOTfiles/yields_Bs_binned_cent_1st_y.root"
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_2ND_Y="ROOTfiles/yields_Bs_binned_cent_2nd_y.root"
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_Y="ROOTfiles/yields_Bs_binned_y.root"
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_PT="ROOTfiles/yields_Bs_binned_pt.root"
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_PT_trk="ROOTfiles/yields_Bs_binned_pt_trk.root"
OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_MULT="ROOTfiles/yields_Bs_binned_Mult.root"

NPROOFIT_PbPb="1"
NPROOFIT_PbPb_BP="701.019629*TMath::Erf((Bmass-5.140349)/-0.035471)+701.019629+16.946432*TMath::Gaus(Bmass,5.343914,0.040000)/(sqrt(2*3.14159)*0.040000)"

#NPROOFIT_PbPb="1.299998*TMath::Gaus(BmassNewNew,6.099828,-0.242801)/(sqrt(2*3.14159)*-0.242801)+8.186179*TMath::Erf((BmassNew-5.000000)/-0.205218)+8.186179+1.263652*(0.426611*TMath::Gaus(BmassNew,5.383307,0.249980)/(sqrt(2*3.14159)*0.249980)+(1-0.426611)*TMath::Gaus(Bmass,5.383307,0.037233)/(sqrt(2*3.14159)*0.037233))"


if [ $DOANALYSISPbPb_ROOFIT_BINNED_DOUBLE_BS_1ST  -eq 1  ]; then
root -b  -q 'roofitB.C+('1','\"ntphi\"','0','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"Bpt\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_1ST\"','\"plotFits/doubly_roofit\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_ROOFIT_BINNED_DOUBLE_BS_2ND  -eq 1  ]; then
root -b  -q 'roofitB.C+('2','\"ntphi\"','0','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"Bpt\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_2ND\"','\"plotFits/doubly_roofit\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_ROOFIT_BINNED_DOUBLE_BS_1ST_Y  -eq 1  ]; then
root -b  -q 'roofitB.C+('1','\"ntphi\"','0','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"abs\(By\)\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_1ST_Y\"','\"plotFits/doubly_roofit\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_ROOFIT_BINNED_DOUBLE_BS_2ND_Y  -eq 1  ]; then
root -b  -q 'roofitB.C+('2','\"ntphi\"','0','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"abs\(By\)\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_DOUBLE_2ND_Y\"','\"plotFits/doubly_roofit\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_ROOFIT_BINNED_PT_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('0','\"ntphi\"','0','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"Bpt\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_PT\"','\"results/Bs/Bpt\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

if [ $DOANALYSISPbPb_ROOFIT_BINNED_PT_BS_TRK -eq 1  ]; then
    cut_trk_tight="(track>1)"
    root -b  -q 'roofitB.C('0','\"ntphi\"','0','1','0', \
'\"$INPUTDATAPbPbCANDWISE_BS\"', \
'\"$INPUTMCPbPbCANDWISE_BS\"','\"Bpt\"', \
'\"$TRGPbPb\"','\"$cut_trk_tight\"','\"$SELGENPbPb\"', \
'$ISMCPbPb','1','$ISDOWEIGHTPbPb', \
'\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_PT_trk\"', \
'\"results/Bs/trk_tight_roofit\"','\"$NPROOFIT_PbPb\"', \
'0')' |& tee binned_pt_trk.log
fi

if [ $DOANALYSISPbPb_ROOFIT_BINNED_Y_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('0','\"ntphi\"','0','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"By\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_Y\"','\"results/Bs/By\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi


if [ $DOANALYSISPbPb_ROOFIT_BINNED_MULT_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('1','\"ntphi\"','0','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"nMult\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_BINNED_MULT\"','\"results/Bs/nMult\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi


if [ $DOANALYSISPbPb_ROOFIT_FULL_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('0','\"ntphi\"','1','1','0','\"$INPUTDATAPbPbCANDWISE_BS\"','\"$INPUTMCPbPbCANDWISE_BS\"','\"Bpt\"','\"$TRGPbPb\"','\"$CUTPbPb\"','\"$SELGENPbPb\"','$ISMCPbPb','1','$ISDOWEIGHTPbPb','\"$OUTPUTFILEPbPbSAVEHIST_ROOFIT_BS_FULL\"','\"results/Bs/Bpt\"','\"$NPROOFIT_PbPb\"','0')'

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
fi

