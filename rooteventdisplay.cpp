

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
      t->SetMarkerStyle(20);
      t->SetMarkerColor(2); //red 2
      t->Draw("x:y:z");
      t->SetMarkerColor(3); //green 3
      t->Draw("x:y:z","parent!=-211","SAME"); //pion en rouge et kaon en vert
    }
  if (withComponent)
    {
      t->SetMarkerStyle(4);
      t->SetMarkerSize(2);
      t->SetMarkerColor(4); //blue 4
      if (! withParent) t->Draw("x:y:z"); 
      else t->Draw("x:y:z","component==0","SAME"); 
      t->SetMarkerColor(5); //yellow 5
      t->Draw("x:y:z","component==1","SAME"); 
    }
  t->SetMarkerStyle(1);
  t->SetMarkerSize(1);
  t->SetMarkerColor(0); //black
  t->SetTitle("");
  t->SetName("");
  return c;
}
