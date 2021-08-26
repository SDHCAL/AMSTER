/*
 * example:
 * root
 * .L rootstats.cpp
 * TNtuple * t=readevent("d20evt1000stats.txt")
 * t->SetMarkerStyle(20)
 * t->Draw("efficiency:purity")
 * t->SetMarkerColor(2)
 * t->Draw("efficiency:purity","leakedenergy>8","SAME")
 */

TNtuple* readevent(const char* filename)
{
   ifstream in;
   in.open(filename);
   TString rootfilename(filename);
   rootfilename.ReplaceAll(".txt",".root");
   TFile *f = new TFile(rootfilename,"RECREATE");
   Float_t eventnumber,efficiency,purity,Nhits,Nedgesstart,Nedgesend,depositedenergy,leakedenergy,Erecotot,Erecocomp0,Erecopion,Nhitspion,Nhitskaon, timetot;
   TNtuple *ntuple = new TNtuple
   (
      "ntuple",
      "event stat",
      "eventnumber:efficiency:purity:Nhits:Nedgesstart:Nedgesend:depositedenergy:leakedenergy:Erecotot:Erecocomp0:Erecopion:Nhitspion:Nhitskaon:timetot"
   );
   while(1)
   {
      in >> eventnumber
         >> efficiency
         >> purity
         >> Nhits
         >> Nedgesstart
         >> Nedgesend
         >> depositedenergy
         >> leakedenergy
         >> Erecotot
         >> Erecocomp0
         >> Erecopion
         >> Nhitspion
         >> Nhitskaon
         >> timetot;
      if (not in.good()) break;
      ntuple->Fill(eventnumber,efficiency,purity,Nhits,Nedgesstart,Nedgesend,depositedenergy,leakedenergy,Erecotot,Erecocomp0,Erecopion,Nhitspion,Nhitskaon,timetot);
   }
   in.close();
   f->Write();
   return ntuple;
}
