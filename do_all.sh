#! /bin/bash

maketnp () {
    ## run TnP and attach them to the unskimmed MC files
    pushd MakeMCTnP/
    root -q -b -l TnPWeight.C'(0)' &
    root -q -b -l TnPWeight.C'(1)' &
    wait
    popd
}

# new bpt reweighting uses FONLL and doesn't depend on acceptance
# function >> MCEff.C

bpYield () {
    ## yield extraction
    pushd BP/RawYieldFits
    # roofitB.C contains 2 By cuts
    source doRoofit.sh
    popd
}

bpEff () {
    ## eff correction (2d map)
    # this is everything from yield to syst, other than MC syst
    pushd BP/EffAna
    # about 1hr
    # takes unskimmed MC sample as input
    echo "Takes BPw.root as input. Make sure it is up to date"
    ls -l BDTWeights/BPw.root
    root -b -l -q MCEff.C'(1,0)' > eff.log # >> bpsyst2d.root
    # CrossSectionAna.C contains 1 By cuts
    root -b -l -q CrossSectionAna.C'(1)'
    root -b -l -q CrossSectionAna.C'(1, 1)'
    root -b -l -q CrossSectionAna.C'(1, 2)'
    # root -b -l -q CrossSectionAnaMult.C'(1)'
    # >> BP/EffAna/FinalFiles/BPPPCorrYieldPT.root
    popd
}


bsYield () {
    cd Bs/RawYieldFits
    # roofitB.C contains 2 By cuts
    source doRoofit.sh
    cd ../..
}

bsEff () {
    pushd Bs/EffAna
    # about 1hr
    echo "Takes Bsw.root as input"
    root -b -l -q MCEff.C'(1,0)' > eff.log
    ls -l BDTWeights/Bsw.root
    root -b -l -q CrossSectionAna.C'(1)'
    # root -b -l -q CrossSectionAnaMult.C'(1)'
    # >> Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root
    popd
}


nominal () {
    # central value
    cd BsBPFinalResults/BsBPRatio
    root -b -l -q PlotBsBPRatio.C'(1)'
    cd ../RAA
    root -b -l -q BPRAA.C
    root -b -l -q BsRAA.C
    cd ..

    cd Comparisons/Fiducial
    root -b -l -q BPComparison.C
    root -b -l -q BsComparison.C
    root -b -l -q BPNewFidNoScale.C
    root -b -l -q BsNewFidNoScale.C
    cd ../../..
}

syst () {
    pushd 2DMapSyst
    root -b -l -q CalEffSystBP.C # >> outfiles/bpsyst2d.root
    root -b -l -q CalEffSystBs.C
    root -b -l -q PlotEffSyst2D.C'(0)'
    root -b -l -q PlotEffSyst2D.C'(1)'
    # need to pass the above to cross calculation
    popd
}

## MC Stat Systematics
# takes 2D map eff as input
bpStat () {
    pushd MCStatSyst/BP
    root -b -l -q Generate2DMaps.C
    # more than 2hr
    root -b -l -q MCStatCal.C > mcstat.log
    popd
}

bsStat() {
    cd MCStatSyst/Bs
    root -b -l -q Generate2DMaps.C
    # ~1hr
    root -b -l -q MCStatCal.C > mcstat.log
    cd ../..
}

comp () {
    # get pdf variation errors
    python master.py
    # comparison plot again
    pushd BsBPFinalResults/Comparisons/Fiducial/
    # << BP/EffAna/FinalFiles/BPPPCorrYieldPT.root
    root -b -l -q BPComparison.C
    root -b -l -q BsComparison.C
    python syst_table.py
    cd ../../RAA/
    root -b -l -q BPRAA.C
    root -b -l -q BsRAA.C
    popd
}

paperPlots () {
    # input
    pushd MakeFinalPlots/NominalPlots/CrossSection
    root -b -l -q plotPt.C'(1,1,0,1,1)'
    cd ../RAA
    root -b -l -q plotPt.C'(1,1,0,1,1)'
    popd
}
# maketnp

# bpYield &
# bsYield &
# wait

bpEff &
bsEff &
wait

nominal
syst

# bpStat&
# bsStat&
# wait

comp
paperPlots
