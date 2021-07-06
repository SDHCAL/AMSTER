#include "hamsters.hpp"   /* headers                         */
#include "paramsters.hpp" /* parameters and global variables */
#include "amster.hpp"     /* data input and MST              */
#include "amsterithm.hpp" /* clustering algorithm            */

int
main(int argc, char** argv)
{
   createdirectories<E,V>();
   lcio::LCReader* lcReader = lcio::LCFactory::getInstance()->createLCReader();
   lcReader->open(argv[1]);
   lcio::LCEvent* evt;
   int eventnumber = 0;
   float leakedenergy = 0.;
   float depositedenergy = 0.;
   std::vector<E> mean(10,0.0);
   while((evt = lcReader->readNextEvent()) != 0)
   {
      leakedenergy = evt->getParameters().getFloatVal("LeakedEnergy");
      depositedenergy = evt->getParameters().getFloatVal("DepositedEnergy");
      //if (leakedenergy > 3) continue;

      t0 = std::chrono::high_resolution_clock::now();
         Graph<E,V>                           g;
         Graph<E,V>                           gmstkruskal;
         std::vector<EDesc<E,V>>              vecmstkruskal;
         std::vector<std::pair<EDesc<E,V>,E>> edgeslistkruskal;
      t1 = std::chrono::high_resolution_clock::now();
         slcioinput<E,V>(evt, g);
      t2 = std::chrono::high_resolution_clock::now();
         boost::kruskal_minimum_spanning_tree(g, std::back_inserter(vecmstkruskal));
      t3 = std::chrono::high_resolution_clock::now();
         fillmstkruskal<E,V>(g, vecmstkruskal, gmstkruskal);
      t4 = std::chrono::high_resolution_clock::now();
         sortedges<E,V>(gmstkruskal, edgeslistkruskal);
      t5 = std::chrono::high_resolution_clock::now();
         cluster<E,V>(gmstkruskal, edgeslistkruskal, 2);
      t6 = std::chrono::high_resolution_clock::now();

      std::cout << "/********** Event " << eventnumber << " **********/"      << std::endl;
      std::cout << "N_V             = " << boost::num_vertices(g) << std::endl;
      std::cout << "N_Estart        = " << boost::num_edges(g) << std::endl;
      std::cout << "N_Eend          = " << boost::num_edges(gmstkruskal) << std::endl;
      std::cout << "E_tot           = " << hcalenergy<E,V>(gmstkruskal)  << std::endl;
      std::cout << "E_0 | E_MC0     = " << hcalenergycomponent<E,V>(gmstkruskal,0) << " | " << gmstkruskal[0].energy << std::endl;
      std::cout << "E_1 | E_MC1     = " << hcalenergycomponent<E,V>(gmstkruskal,1) << " | " << gmstkruskal[1].energy << std::endl;
      std::cout << "E_deposited     = " << depositedenergy << std::endl;
      std::cout << "E_leaked        = " << leakedenergy << std::endl;
      std::cout << "Efficiency      = " << efficiency(gmstkruskal) << std::endl;
      std::cout << "Purity          = " << purity(gmstkruskal) << std::endl;
      std::cout << "Parents c0      = " << countparents<E,V>(gmstkruskal,0)[0] << " x (-211) "
                                        << countparents<E,V>(gmstkruskal,0)[1] << " x 311 "
                                        << countparents<E,V>(gmstkruskal,0)[2] << " others"
                                        << std::endl;
      std::cout << "Parents c1      = " << countparents<E,V>(gmstkruskal,1)[0] << " x (-211) "
                                        << countparents<E,V>(gmstkruskal,1)[1] << " x 311 "
                                        << countparents<E,V>(gmstkruskal,1)[2] << " others"
                                        << std::endl;
      std::cout << std::endl;
      exportdata<E,V>("graph", gmstkruskal, 0, eventnumber);

      mean[0] += boost::num_vertices(g);
      mean[1] += boost::num_edges(g);
      mean[2] += boost::num_edges(gmstkruskal);
      mean[3] += hcalenergy<E,V>(gmstkruskal);
      mean[4] += hcalenergycomponent<E,V>(gmstkruskal,0);
      mean[5] += hcalenergycomponent<E,V>(gmstkruskal,1);
      mean[6] += gmstkruskal[0].energy;
      mean[7] += gmstkruskal[1].energy;
      mean[8] += efficiency(gmstkruskal);
      mean[9] += purity(gmstkruskal);
      t[0] += std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t0).count()*1e-9;
      t[1] += std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()*1e-9;
      t[2] += std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()*1e-9;
      t[3] += std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()*1e-9;
      t[4] += std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()*1e-9;
      t[5] += std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()*1e-9;
      t[6] += std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()*1e-9;
      ++eventnumber;
   }

   std::cout << "/********** End **********/"      << std::endl;
   std::cout << "N_events        = " << std::fixed << eventnumber              << std::endl;
   std::cout << "<N_V>           = " << std::fixed << mean[0]/eventnumber      << std::endl;
   std::cout << "<N_Estart>      = " << std::fixed << mean[1]/eventnumber      << std::endl;
   std::cout << "<N_Eend>        = " << std::fixed << mean[2]/eventnumber      << std::endl;
   std::cout << "<E_tot>         = " << std::fixed << mean[3]/eventnumber      << std::endl;
   std::cout << "<E_0> | <E_MC0> = " << std::fixed << mean[4]/eventnumber << " | " << mean[6]/eventnumber << std::endl;
   std::cout << "<E_1> | <E_MC1> = " << std::fixed << mean[5]/eventnumber << " | " << mean[7]/eventnumber << std::endl;
   std::cout << "<Efficiency>    = " << std::fixed << mean[8]/eventnumber      << std::endl;
   std::cout << "<Purity>        = " << std::fixed << mean[9]/eventnumber      << std::endl;
   std::cout << "t_tot           = " << std::fixed << t[0]             << " s" << std::endl;
   std::cout << "<t/event>       = " << std::fixed << t[0]/eventnumber << " s" << std::endl;
   std::cout << "<t/variables>   = " << std::fixed << t[1]/eventnumber << " s" << std::endl;
   std::cout << "<t/input>       = " << std::fixed << t[2]/eventnumber << " s" << std::endl;
   std::cout << "<t/kruskal>     = " << std::fixed << t[3]/eventnumber << " s" << std::endl;
   std::cout << "<t/fillmst>     = " << std::fixed << t[4]/eventnumber << " s" << std::endl;
   std::cout << "<t/sortedges>   = " << std::fixed << t[5]/eventnumber << " s" << std::endl;
   std::cout << "<t/cluster>     = " << std::fixed << t[6]/eventnumber << " s" << std::endl;
   
   lcReader->close();
   delete lcReader;

   return EXIT_SUCCESS;
}
