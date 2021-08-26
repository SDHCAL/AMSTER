#include "hamsters.hpp"   /* headers                         */
#include "paramsters.hpp" /* parameters and global variables */
#include "amster.hpp"     /* data input and MST              */
#include "amsterithm.hpp" /* clustering algorithm            */

int
main(int argc, char** argv)
{
   lcio::LCReader* lcReader = lcio::LCFactory::getInstance()->createLCReader();
   lcio::LCEvent*  evt;
   lcReader->open(argv[1]);

   createdirectories<E,V>();
   std::ofstream fout("d20evt100statskruskaldir.txt");

   while((evt = lcReader->readNextEvent()) != 0 && eventnumber < 100)
   //while((evt = lcReader->readNextEvent()) != 0)
   {
      leakedenergy    = evt->getParameters().getFloatVal("LeakedEnergy");
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
         cluster<E,V>(gmstkruskal, edgeslistkruskal, hcalenergy<E,V>(gmstkruskal)*0.01);
      t6 = std::chrono::high_resolution_clock::now();

//      std::cout << "/********** Event " << eventnumber << " **********/"      << std::endl;
/*      std::cout << "N_V                  = " << boost::num_vertices(g) << std::endl;
      std::cout << "N_Estart             = " << boost::num_edges(g) << std::endl;
      std::cout << "N_Eend               = " << boost::num_edges(gmstkruskal) << std::endl;
      std::cout << "E_tot                = " << hcalenergy<E,V>(gmstkruskal)  << std::endl;
      std::cout << "E_0 | E_0max | E_MC0 = " << hcalenergycomponent<E,V>(gmstkruskal,0) << " | " << hcalenergyparent<E,V>(gmstkruskal, -211) << " | " << gmstkruskal[0].energy << std::endl;
      std::cout << "E_1 | E_1max | E_MC0 = " << hcalenergycomponent<E,V>(gmstkruskal,1) << " | " << hcalenergyparent<E,V>(gmstkruskal, 311) << " | " << gmstkruskal[1].energy << std::endl;
      std::cout << "E_deposited          = " << depositedenergy << std::endl;
      std::cout << "E_leaked             = " << leakedenergy << std::endl;
      std::cout << "Efficiency           = " << efficiency<E,V>(gmstkruskal) << std::endl;
      std::cout << "Purity | Puritymin   = " << purity<E,V>(gmstkruskal) << " | " << minimumpurity<E,V>(gmstkruskal) << std::endl;
      std::cout << "Parents c0           = " << numberparentscomponent<E,V>(gmstkruskal,0)[0] << " x (-211) "
                                             << numberparentscomponent<E,V>(gmstkruskal,0)[1] << " x 311 "
                                             << numberparentscomponent<E,V>(gmstkruskal,0)[2] << " others"
                                             << std::endl;
      std::cout << "Parents c1           = " << numberparentscomponent<E,V>(gmstkruskal,1)[0] << " x (-211) "
                                             << numberparentscomponent<E,V>(gmstkruskal,1)[1] << " x 311 "
                                             << numberparentscomponent<E,V>(gmstkruskal,1)[2] << " others"
                                             << std::endl;
      std::cout << std::endl;
*/
      exportdata<E,V>("graph", gmstkruskal, eventnumber);
      fout << eventnumber                              << " "
           << efficiency<E,V>(gmstkruskal)             << " "
           << purity<E,V>(gmstkruskal)                 << " "
           << boost::num_vertices(g)                   << " "
           << boost::num_edges(g)                      << " "
           << boost::num_edges(gmstkruskal)            << " "
           << depositedenergy                          << " "
           << leakedenergy                             << " "
           << hcalenergy<E,V>(gmstkruskal)             << " "
           << hcalenergycomponent<E,V>(gmstkruskal,0)  << " "
           << hcalenergyparent<E,V>(gmstkruskal, -211) << " "
           << numberparents<E,V>(gmstkruskal,-211)      << " "
           << numberparents<E,V>(gmstkruskal,311)       << " "
           << std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t0).count()*1e-9
           << '\n';

      m[0]  += boost::num_vertices(g);
      m[1]  += boost::num_edges(g);
      m[2]  += boost::num_edges(gmstkruskal);
      m[3]  += hcalenergy<E,V>(gmstkruskal);
      m[4]  += hcalenergycomponent<E,V>(gmstkruskal,0);
      m[5]  += hcalenergycomponent<E,V>(gmstkruskal,1);
      m[6]  += hcalenergyparent<E,V>(gmstkruskal, -211);
      m[7]  += hcalenergyparent<E,V>(gmstkruskal, 311);
      m[8]  += gmstkruskal[0].energy;
      m[9]  += gmstkruskal[1].energy;
      m[10] += efficiency(gmstkruskal);
      m[11] += purity(gmstkruskal);
      t[0]  += std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t0).count()*1e-9;
      t[1]  += std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count()*1e-9;
      t[2]  += std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()*1e-9;
      t[3]  += std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count()*1e-9;
      t[4]  += std::chrono::duration_cast<std::chrono::nanoseconds>(t4-t3).count()*1e-9;
      t[5]  += std::chrono::duration_cast<std::chrono::nanoseconds>(t5-t4).count()*1e-9;
      t[6]  += std::chrono::duration_cast<std::chrono::nanoseconds>(t6-t5).count()*1e-9;
      ++eventnumber;
   }

   fout.close();
   std::cout << "/********** End **********/"      << std::endl;
   std::cout << "N_events                   = " << std::fixed << eventnumber              << std::endl;
   std::cout << "<N_V>                      = " << std::fixed << m[0]/eventnumber         << std::endl;
   std::cout << "<N_Estart>                 = " << std::fixed << m[1]/eventnumber         << std::endl;
   std::cout << "<N_Eend>                   = " << std::fixed << m[2]/eventnumber         << std::endl;
   std::cout << "<E_tot>                    = " << std::fixed << m[3]/eventnumber         << std::endl;
   std::cout << "<E_0> | <E_0max> | <E_MC0> = " << std::fixed << m[4]/eventnumber << " | " << m[6]/eventnumber << " | " << m[8]/eventnumber << std::endl;
   std::cout << "<E_1> | <E_1max> | <E_MC1> = " << std::fixed << m[5]/eventnumber << " | " << m[7]/eventnumber << " | " << m[9]/eventnumber << std::endl;
   std::cout << "<Efficiency>               = " << std::fixed << m[10]/eventnumber        << std::endl;
   std::cout << "<Purity>                   = " << std::fixed << m[11]/eventnumber        << std::endl;
   std::cout << "t_tot                      = " << std::fixed << t[0]             << " s" << std::endl;
   std::cout << "<t/event>                  = " << std::fixed << t[0]/eventnumber << " s" << std::endl;
   std::cout << "<t/variables>              = " << std::fixed << t[1]/eventnumber << " s" << std::endl;
   std::cout << "<t/input>                  = " << std::fixed << t[2]/eventnumber << " s" << std::endl;
   std::cout << "<t/kruskal>                = " << std::fixed << t[3]/eventnumber << " s" << std::endl;
   std::cout << "<t/fillmst>                = " << std::fixed << t[4]/eventnumber << " s" << std::endl;
   std::cout << "<t/sortedges>              = " << std::fixed << t[5]/eventnumber << " s" << std::endl;
   std::cout << "<t/cluster>                = " << std::fixed << t[6]/eventnumber << " s" << std::endl;
   
   lcReader->close();
   delete lcReader;

   return EXIT_SUCCESS;
}
