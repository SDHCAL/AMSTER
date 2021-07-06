/*
 * functions for graph management
 */

template <typename E, typename V>
struct vertexproperties
{
   bool track;
   V    parent;
   V    component;
   E    energy;
   E    xpos;
   E    ypos;
   E    zpos;
};

template <typename E, typename V>
using Graph = typename boost::adjacency_list
<
   boost::vecS,                            /* edge container    */
   boost::vecS,                            /* vertex container  */
   boost::undirectedS,                     /* undirected graph  */
   vertexproperties<E,V>,                  /* vertex properties */
   boost::property<boost::edge_weight_t,E> /* edge property     */
>;

template <typename E, typename V>
using EIter = typename boost::graph_traits<Graph<E,V>>::edge_iterator;

template <typename E, typename V>
using EDesc = typename boost::graph_traits<Graph<E,V>>::edge_descriptor;

template <typename E, typename V>
using VIter = typename boost::graph_traits<Graph<E,V>>::vertex_iterator;

template <typename E, typename V>
using VDesc = typename boost::graph_traits<Graph<E,V>>::vertex_descriptor;

template <typename E, typename V>
void
slcioinput(lcio::LCEvent* evt, Graph<E,V> &graph)
{
   V numberofhits = 0;
   V component = 0;

   lcio::LCCollection* col1 = evt->getCollection("primaryParticles");
   for (unsigned int i=0, N=col1->getNumberOfElements(); i<N; ++i)
   {
      lcio::MCParticle* particle = (lcio::MCParticle*) col1->getElementAt(i);
      boost::add_vertex(graph);
      graph[i].track     = true;
      graph[i].parent    = particle->getPDG();
      graph[i].component = component;
      graph[i].energy    = particle->getEnergy();
      graph[i].xpos      = particle->getVertex()[0];
      graph[i].ypos      = particle->getVertex()[1];
      graph[i].zpos      = particle->getVertex()[2];
      ++component;
      ++numberofhits;
   }
   int j = col1->getNumberOfElements();
   lcio::LCCollection* col2 = evt->getCollection("HCALEndcap");
   lcio::LCCollection* col3 = evt->getCollection("RelationParticleToHit");
   lcio::LCRelationNavigator navi(col3);
   for (unsigned int i=0, N=col2->getNumberOfElements(); i<N; ++i)
   {
      lcio::CalorimeterHit* particle = (lcio::CalorimeterHit*) col2->getElementAt(i);
      lcio::LCObjectVec MCpartVec = navi.getRelatedFromObjects(particle);
      int parent = 0;
      for (unsigned int k=0; k<MCpartVec.size();++k)
      {
         lcio::MCParticle* particle = (lcio::MCParticle*) MCpartVec[k];
         parent *= 10000;
         parent += particle->getPDG();
      }
      boost::add_vertex(graph);
      graph[j+i].track     = false;
      graph[j+i].parent    = parent;
      graph[j+i].component = 0;
      graph[j+i].energy    = particle->getEnergy();
      graph[j+i].xpos      = particle->getPosition()[0] - 504.788; //- 494.38;
      graph[j+i].ypos      = particle->getPosition()[1] - 504.788; //- 494.38;
      graph[j+i].zpos      = 0.1*particle->getPosition()[2] / 26.131;
      ++numberofhits; 
   }

   E distance;
   for (int i=0; i!=numberofhits; ++i)
   {
      for (int j=0; j<i && j!=i; ++j)
      {
         if
         (
            (not (graph[i].track == true && graph[j].track == true))
            && std::abs(graph[i].zpos-graph[j].zpos) < 0.4
            //(graph[i].parent == graph[j].parent)
            //&& (graph[i].component == graph[j].component)
         )
         {
            distance = sqrt(pow((graph[i].xpos-graph[j].xpos),2.0)
                           +pow((graph[i].ypos-graph[j].ypos),2.0)
                           +pow((graph[i].zpos-graph[j].zpos),2.0));
            //if (distance < )
            boost::add_edge(j,i,distance,graph);
         }
      }
   }
}

template <typename E, typename V>
void
fillmstkruskal(Graph<E,V> &gin, std::vector<EDesc<E,V>> &vecmstkruskal, Graph<E,V> &gout)
{
   for (typename std::vector<EDesc<E,V>>::iterator eiter = vecmstkruskal.begin(); eiter != vecmstkruskal.end(); ++eiter)
   {
      boost::add_edge(boost::source(*eiter, gin),
                      boost::target(*eiter, gin),
                      boost::get(boost::edge_weight, gin, *eiter),
                      gout);

      gout[boost::source(*eiter, gin)].track     = gin[boost::source(*eiter, gin)].track;
      gout[boost::source(*eiter, gin)].parent    = gin[boost::source(*eiter, gin)].parent;
      gout[boost::source(*eiter, gin)].component = gin[boost::source(*eiter, gin)].component;
      gout[boost::source(*eiter, gin)].energy    = gin[boost::source(*eiter, gin)].energy;
      gout[boost::source(*eiter, gin)].xpos      = gin[boost::source(*eiter, gin)].xpos;
      gout[boost::source(*eiter, gin)].ypos      = gin[boost::source(*eiter, gin)].ypos;
      gout[boost::source(*eiter, gin)].zpos      = gin[boost::source(*eiter, gin)].zpos;

      gout[boost::target(*eiter, gin)].track     = gin[boost::target(*eiter, gin)].track;
      gout[boost::target(*eiter, gin)].parent    = gin[boost::target(*eiter, gin)].parent;
      gout[boost::target(*eiter, gin)].component = gin[boost::target(*eiter, gin)].component;
      gout[boost::target(*eiter, gin)].energy    = gin[boost::target(*eiter, gin)].energy;
      gout[boost::target(*eiter, gin)].xpos      = gin[boost::target(*eiter, gin)].xpos;
      gout[boost::target(*eiter, gin)].ypos      = gin[boost::target(*eiter, gin)].ypos;
      gout[boost::target(*eiter, gin)].zpos      = gin[boost::target(*eiter, gin)].zpos;
   }
}

template <typename E, typename V>
void
fillmstprim(Graph<E,V> &graphin, std::vector<VDesc<E,V>> &vecmstprim, Graph<E,V> &graphout)
{
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graphin); viter != viterend; ++viter)
   {
      if (vecmstprim[*viter] != *viter)
      {
         boost::add_edge(vecmstprim[*viter],
                         *viter,
                         boost::get(boost::edge_weight, graphin, boost::edge(vecmstprim[*viter], *viter, graphin).first),
                         graphout);
         graphout[*viter].track     = graphin[*viter].track;
         graphout[*viter].parent    = graphin[*viter].parent;
         graphout[*viter].component = graphin[*viter].component;
         graphout[*viter].energy    = graphin[*viter].energy;
         graphout[*viter].xpos      = graphin[*viter].xpos;
         graphout[*viter].ypos      = graphin[*viter].ypos;
         graphout[*viter].zpos      = graphin[*viter].zpos;
      }
      else //if (vecmstprim[*viter] == *viter)
      {
         boost::add_vertex(graphout);
         graphout[*viter].track     = graphin[*viter].track;
         graphout[*viter].parent    = graphin[*viter].parent;
         graphout[*viter].component = graphin[*viter].component;
         graphout[*viter].energy    = graphin[*viter].energy;
         graphout[*viter].xpos      = graphin[*viter].xpos;
         graphout[*viter].ypos      = graphin[*viter].ypos;
         graphout[*viter].zpos      = graphin[*viter].zpos;
      }
   }
}

template <typename E, typename V>
void
createdirectories()
{
    if (mkdir("graphtxt", 0777) == -1);
    if (mkdir("figures", 0777) == -1);
    if (mkdir("figuresexample", 0777) == -1);
}

template <typename E, typename V>
void
exportdata(std::string graphname, Graph<E,V> graph, V componentnumber, int eventnumber)
{
   std::string output = "graphtxt/" + graphname + std::to_string(eventnumber) + ".txt";
   std::ofstream fout(output);
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      fout << *viter                  << " "
           << graph[*viter].parent    << " "
           << graph[*viter].component << " "
           << graph[*viter].energy    << " "
           << graph[*viter].xpos      << " "
           << graph[*viter].ypos      << " "
           << graph[*viter].zpos      << " "
           << '\n';
   }
   fout.close();
}
