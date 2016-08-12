#!/bin/bash


eosfolder=/store/group/dpg_hcal/comm_hcal/RecoAlgos/Summer16Method2Update/HcalRecoTesting/HighPtJet80/0000/
#eosfolder=/store/group/dpg_hcal/comm_hcal/RecoAlgos/Summer16Method2Update/HcalRecoTesting/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/0000/
#eosfolder=/store/group/dpg_hcal/comm_hcal/RecoAlgos/Summer16Method2Update/HcalRecoTesting/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_new/0000/


outfolder=/afs/cern.ch/user/h/hum/work/public/CMSSW_8_0_1/src/PulseShapes/Datatrans/
#outfolder=/afs/cern.ch/work/j/jlawhorn/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/

for file in `/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls ${eosfolder} | grep root`
do
    #if [[ -e ${outfolder}/${file} ]]; then
    #echo "Output file exists. Not submitting."
    #else
    #echo run.sh `pwd` $file $eosfolder $outfolder
    bsub -o out.%J -q 8nm run.sh ${CMSSW_BASE} $file $eosfolder $outfolder
    #fi
done
