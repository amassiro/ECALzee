Run
====

    cmsRun dumpMC.py  inputFiles=file:/eos/cms/store/relval/CMSSW_9_2_14_patch1/RelValTTbar_13/GEN-SIM-RECO/92X_upgrade2017_realistic_v14_tst-v1/10000/BEF26E0E-C8C5-E711-8877-0025905B85FC.root   outputFile=test.root
    cmsRun dumpMC.py  inputFiles=file:/tmp/amassiro/BEF26E0E-C8C5-E711-8877-0025905B85FC.root   outputFile=test.root
    
    
    /DoubleEG/Run2017C-ZElectron-PromptReco-v1/RAW-RECO
    cmsRun dumpDATA.py  inputFiles=file:/eos/cms/store/data/Run2017C/DoubleEG/RAW-RECO/ZElectron-PromptReco-v1/000/299/368/00000/1E1F5306-7C6D-E711-A4F7-02163E01A1D9.root  outputFile=test.root
    cmsRun dumpDATA.py  inputFiles=file:/tmp/amassiro/1E1F5306-7C6D-E711-A4F7-02163E01A1D9.root   outputFile=test.data.root  maxEvents=-1
    
    
    cp /eos/cms/store/data/Run2017C/DoubleEG/RAW-RECO/ZElectron-PromptReco-v1/000/299/368/00000/566CA2BF-8E6D-E711-A37F-02163E012531.root    /tmp/amassiro/
    cmsRun dumpDATA.py  inputFiles=file:/tmp/amassiro/566CA2BF-8E6D-E711-A37F-02163E012531.root   outputFile=test.data.2.root  maxEvents=-1

    cp /eos/cms/store/data/Run2017C/DoubleEG/RAW-RECO/ZElectron-PromptReco-v1/000/299/368/00000/64C2C610-926D-E711-A8E6-02163E0128FE.root     /tmp/amassiro/
    cmsRun dumpDATA.py  inputFiles=file:/tmp/amassiro/64C2C610-926D-E711-A8E6-02163E0128FE.root   outputFile=test.data.3.root  maxEvents=-1

    
    

Plot
====


    tree = (TTree*) _file0->Get("TreeProducer/tree")
    tree->Draw("mll", "mll>10")
    tree->Draw("mll")
    
    tree->Draw("@std_vector_Ele_pt.size()", "")
    
    
    tree->Draw("std_vector_Ele_pt.at(0)", "")
    tree->Draw("std_vector_Ele_pt>10", "")

    
    tree->Draw("mll", "mll>10 && (@std_vector_Ele_pt.size()>=2) && (std_vector_Ele_pt.at(0)>30)")
    