#include "hamsters.hpp"   /* headers                         */
#include "paramsters.hpp" /* parameters and global variables */
#include "amster.hpp"     /* data input and MST              */
#include "amsterithm.hpp" /* clustering algorithm            */
#include "example.hpp"    /* example functions               */

int
main()
{
   Graph<E,V> g;
   Graph<E,V> gmstkruskal;
   Graph<E,V> gmstprim;
   std::vector<EDesc<E,V>> vecmstkruskal;
   std::vector<VDesc<E,V>> vecmstprim;
   std::vector<std::pair<EDesc<E,V>,E>> edgeslistkruskal;
   std::vector<std::pair<EDesc<E,V>,E>> edgeslistprim;

   exampleinput<E,V>(hits, trackenergy, g);

   boost::kruskal_minimum_spanning_tree(g, std::back_inserter(vecmstkruskal));
   fillmstkruskal<E,V>(g, vecmstkruskal, gmstkruskal);
   
   vecmstprim.resize(boost::num_vertices(g));
   boost::prim_minimum_spanning_tree(g, &vecmstprim[0]);
   fillmstprim<E,V>(g, vecmstprim, gmstprim);

   sortedges<E,V>(gmstkruskal, edgeslistkruskal);
   sortedges<E,V>(gmstprim,    edgeslistprim);

   if (hits < max)
   {
      createdirectories<E,V>();
      printgraph<E,V>("figuresexample/1g",            format, g,           -1);
      printgraph<E,V>("figuresexample/21gmstkruskal", format, gmstkruskal, -1);
      printgraph<E,V>("figuresexample/31gmstprim",    format, gmstprim,    -1);
   } 
 
   examplecluster<E,V>(gmstkruskal, edgeslistkruskal, deltaE);
   examplecluster<E,V>(gmstprim,    edgeslistprim,    deltaE);

   if (hits < max)
   {
      printgraph<E,V>("figuresexample/22gmstkruskalcluster", format, gmstkruskal, -1);
      printgraph<E,V>("figuresexample/32gmstprimcluster",    format, gmstprim,    -1);
   }

   std::cout << "/********** Kruskal **********/"                       << std::endl;
   std::cout << "N_V         = " << boost::num_vertices(g)              << std::endl;
   std::cout << "N_Estart    = " << boost::num_edges(g)                 << std::endl;
   std::cout << "N_Eend      = " << boost::num_edges(gmstkruskal)       << std::endl;
   std::cout << "E_tot       = " << examplehcalenergy<E,V>(gmstkruskal) << std::endl;
   std::cout << "E_0 | E_MC0 = " << examplehcalenergycomponent<E,V>(gmstkruskal,0) << " | " << gmstkruskal[0].energy << std::endl;
   std::cout << "E_1 | E_MC1 = " << examplehcalenergycomponent<E,V>(gmstkruskal,1) << " | " << gmstkruskal[1].energy << std::endl;
   std::cout << "Efficiency  = " << exampleefficiency(gmstkruskal)      << std::endl;
   std::cout << "Purity      = " << examplepurity(gmstkruskal)          << std::endl;
   std::cout << "Parents c0  = " << examplenumberparentscomponent<E,V>(gmstkruskal,0,0) << " x (0) "
                                 << examplenumberparentscomponent<E,V>(gmstkruskal,1,0) << " x (1) "
                                 << std::endl;
   std::cout << "Parents c1  = " << examplenumberparentscomponent<E,V>(gmstkruskal,0,1) << " x (0) "
                                 << examplenumberparentscomponent<E,V>(gmstkruskal,1,1) << " x (1) "
                                 << std::endl;
   std::cout << std::endl;
   std::cout << "/********** Prim **********/"                       << std::endl;
   std::cout << "N_V         = " << boost::num_vertices(g)           << std::endl;
   std::cout << "N_Estart    = " << boost::num_edges(g)              << std::endl;
   std::cout << "N_Eend      = " << boost::num_edges(gmstprim)       << std::endl;
   std::cout << "E_tot       = " << examplehcalenergy<E,V>(gmstprim) << std::endl;
   std::cout << "E_0 | E_MC0 = " << examplehcalenergycomponent<E,V>(gmstprim,0) << " | " << gmstprim[0].energy << std::endl;
   std::cout << "E_1 | E_MC1 = " << examplehcalenergycomponent<E,V>(gmstprim,1) << " | " << gmstprim[1].energy << std::endl;
   std::cout << "Efficiency  = " << exampleefficiency(gmstprim)      << std::endl;
   std::cout << "Purity      = " << examplepurity(gmstprim)          << std::endl;
   std::cout << "Parents c0  = " << examplenumberparentscomponent<E,V>(gmstprim,0,0) << " x (0) "
                                 << examplenumberparentscomponent<E,V>(gmstprim,1,0) << " x (1) "
                                 << std::endl;
   std::cout << "Parents c1  = " << examplenumberparentscomponent<E,V>(gmstprim,0,1) << " x (0) "
                                 << examplenumberparentscomponent<E,V>(gmstprim,1,1) << " x (1) "
                                 << std::endl;
   return EXIT_SUCCESS;
}
