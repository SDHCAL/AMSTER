/*
 * functions for the example
 */

template <typename E, typename V>
void
printgraph(std::string graphname, std::string format, Graph<E,V> graph, V componentnumber)
{
   EIter<E,V> eiter, eiterend;
   std::string graphdot = graphname + ".dot";
   std::string output   = graphname + "." + format;
   std::ofstream fout(graphdot);
   fout << "graph A {\n"
        << " rankdir=LR\n"
        << " size=\"10,10\"\n"
        << " ratio=\"filled\"\n"
        << " edge[style=\"bold\"]\n"
        << " node[shape=\"circle\"]\n";
   for (boost::tie(eiter, eiterend) = boost::edges(graph); eiter != eiterend; ++eiter)
   {
      if (componentnumber == -1)
      {
         fout << boost::source(*eiter, graph) << " -- " << boost::target(*eiter, graph);
         fout << "[color=\"gray\", label=\""
              << boost::get(boost::edge_weight, graph, *eiter)
              << "\"];\n";
      } 
      if (graph[boost::source(*eiter, graph)].component == componentnumber)
      {
         fout << boost::source(*eiter, graph) << " -- " << boost::target(*eiter, graph);
         fout << "[color=\"gray\", label=\""
              << boost::get(boost::edge_weight, graph, *eiter)
              << "\"];\n";
      }
   }
   fout << "}\n";
   fout.close();
   graphname = "dot -T" + format + " " + graphdot + " > " + output;
   const char* command = graphname.c_str();
   system(command);
}

template <typename E, typename V>
void
exampleinput(int numberofhits, E trackenergy, Graph<E,V> &graph)
{
   boost::add_vertex(graph);
   graph[0].track     = true;
   graph[0].parent    = 0;
   graph[0].component = 0;
   graph[0].energy    = trackenergy;
   graph[0].xpos      = -100.0;
   graph[0].ypos      = 0.0;
   graph[0].zpos      = 0.0;

   boost::add_vertex(graph);
   graph[1].track     = true;
   graph[1].parent    = 1;
   graph[1].component = 1;
   graph[1].energy    = trackenergy*0.85;
   graph[1].xpos      = 100.0;
   graph[1].ypos      = 0.0;
   graph[1].zpos      = 0.0;

   for (int i=2; i<numberofhits; ++i)
   {
      boost::add_vertex(graph);
      graph[i].track     = false;
      if (i%2 == 0)
      {
         graph[i].component = 0;
         graph[i].parent    = 0;
         graph[i].energy    = i*0.618;
         graph[i].xpos      = -100-i;
         graph[i].ypos      = -i*i*0.618;
         graph[i].zpos      = -i*i*0.618;
      }
      else
      {  
         graph[i].component = 0;
         graph[i].parent    = 1;
         graph[i].energy    = i*0.618;
         graph[i].xpos      = 100+i;
         graph[i].ypos      = i*i*0.618;
         graph[i].zpos      = i*i*0.618;
      }
   }

   E distance;
   for (int i=0; i!=numberofhits; ++i)
   {
      for (int j=0; j<i && j!=i; ++j)
      {
         distance = sqrt(pow((graph[i].xpos-graph[j].xpos),2.0)
                        +pow((graph[i].ypos-graph[j].ypos),2.0)
                        +pow((graph[i].zpos-graph[j].zpos),2.0));
         if ((not(graph[i].track == true && graph[j].track == true)))
         boost::add_edge(j,i,distance,graph);
      }
   }

}

template <typename E, typename V>
E
examplehcalenergy(Graph<E,V> &graph)
{
   E totalenergy = 0.0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].track == false)
         totalenergy += graph[*viter].energy;
   }
   return totalenergy;
}

template <typename E, typename V>
E
examplehcalenergycomponent(Graph<E,V> &graph, V componentnumber)
{
   E energycomponent = 0.0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].component == componentnumber && graph[*viter].track == false)
         energycomponent += graph[*viter].energy;
   }
   return energycomponent;
}

template <typename E, typename V>
std::vector<E>
examplecountparents(Graph<E,V> graph, V componentnumber)
{
   std::vector<E> parents(2,0.0);
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].component == componentnumber)
      {
         if (graph[*viter].parent == 0)
            ++parents[0];
         if (graph[*viter].parent == 1)
            ++parents[1];
      }
   }
   return parents;
}
template <typename E, typename V>
E
exampleefficiency(Graph<E,V> graph)
{
   E efficiency = examplecountparents<E,V>(graph,0)[0]/(examplecountparents<E,V>(graph,0)[0]+examplecountparents<E,V>(graph,1)[0]);
   return efficiency;
}

template <typename E, typename V>
E
examplepurity(Graph<E,V> graph)
{
   E purity = examplecountparents<E,V>(graph,0)[0]/numberhitscomponent<E,V>(graph,0);
   return purity;
}

template <typename E, typename V>
void
clusterexample(Graph<E,V> &graph, std::vector<std::pair<EDesc<E,V>,E>> &p, E deltaE)
{
   unsigned int i = 0;
   unsigned int numberofedges = p.size();
   E targetenergy             = graph[0].energy;
   E targetenergyplusdeltaE   = graph[0].energy+deltaE;
   E targetenergyminusdeltaE  = graph[0].energy-deltaE;
   E currentenergy            = examplehcalenergycomponent<E,V>(graph,0);
   while (currentenergy > targetenergyplusdeltaE && i != numberofedges)
   {
      if (currentenergy > targetenergy)
      {
         boost::remove_edge(boost::source(p[i].first, graph), boost::target(p[i].first, graph), graph);
         updatecomponents<E,V>(graph);
         currentenergy = examplehcalenergycomponent<E,V>(graph,0);
      }
      if (currentenergy <= targetenergyminusdeltaE)
      {
         boost::add_edge(boost::source(p[i].first, graph), boost::target(p[i].first, graph), p[i].second, graph);
         updatecomponents<E,V>(graph);
         currentenergy = examplehcalenergycomponent<E,V>(graph,0);
      }
      ++i;
   }
}
