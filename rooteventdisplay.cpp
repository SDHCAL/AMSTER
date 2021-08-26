/*
 * example:
 * root
 * .L rooteventdisplay.cpp
 * TNtuple* t=read("graph.txt")
 * Draw(t)
 */

TNtuple* read(const char* filename)
{
   ifstream in;
   in.open(filename);
   TString rootfilename(filename);
   rootfilename.ReplaceAll(".txt",".root");
   TFile *f = new TFile(rootfilename,"RECREATE");
   Float_t id,parent,component,energy,x,y,z;
   TNtuple *ntuple = new TNtuple("ntuple","event display","id:parent:component:energy:x:y:z");

   while(1)
     {
       in >> id >> parent >> component >> energy >> x >> y >> z;
       if (!in.good()) break;
       ntuple->Fill(id,parent,component,energy,x,y,z);
     }

   in.close();
   f->Write();
   return ntuple;
}

TCanvas* Draw(TNtuple* t,bool withParent=true,bool withComponent=true)
{
  TCanvas* c=new TCanvas();
  if (withParent)
    {
      t->SetMarkerStyle(20); /* big dot   */
      t->SetMarkerColor(2);  /* 2 for red */
      t->Draw("x:y:z");
      t->SetMarkerColor(3);  /* 3 for green */
      t->Draw("x:y:z","parent!=-211","SAME"); /* pion in red and kaon in green */
    }
  if (withComponent)
    {
      t->SetMarkerStyle(4);
      t->SetMarkerSize(2);
      t->SetMarkerColor(4); /* 4 for blue */
      if (! withParent) t->Draw("x:y:z"); 
      else t->Draw("x:y:z","component==0","SAME"); 
      t->SetMarkerColor(5); /* 5 for yellow */
      t->Draw("x:y:z","component==1","SAME"); 
    }
  t->SetMarkerStyle(1);
  t->SetMarkerSize(1);
  t->SetMarkerColor(0); /* 0 for black */
  t->SetTitle("");
  t->SetName("");
  return c;
}
