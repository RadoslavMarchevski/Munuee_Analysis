void merge() { 
   Target = TFile::Open("kplus_pipie_merged.root","RECREATE");
   FileList = new TList();
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie1/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie2/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie3/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie4/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie5/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie6/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie7/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie8/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie9/histos.root"));
   FileList->Add(TFile::Open("/afs/cern.ch/user/r/rmarchev/compact/reader/kplus_pipie10/histos.root"));
   MergeRootfile(Target,FileList);
}
