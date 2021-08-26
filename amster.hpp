/*
 * functions for graph management
 */

template <typename E, typename V>
struct
vertexproperties
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
   boost::directedS,                       /* directed graph    */
   //boost::undirectedS,                     /* undirected graph  */
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

      if (parent == -2109689) 
      graph[j+i].parent    = -211;
      if (parent == 3109789) 
      graph[j+i].parent    = 311;

      graph[j+i].component = 0;
      graph[j+i].energy    = particle->getEnergy();
      graph[j+i].xpos      = particle->getPosition()[0] - 504.788; //- 494.38;
      graph[j+i].ypos      = particle->getPosition()[1] - 504.788; //- 494.38;
      graph[j+i].zpos      = 0.1*particle->getPosition()[2] / 26.131;
      ++numberofhits; 
   }

   E distance;
   E x0, y0, z0, x1, y1, z1;
   for (int i=0; i!=numberofhits; ++i)
   {
      for (int j=0; j<i; ++j)
      {
         if
         (
            (not (graph[i].track == true && graph[j].track == true))
            && (std::abs(graph[i].zpos-graph[j].zpos) < 0.4)
            //&& (graph[i].zpos - graph[j].zpos <= 0.)
            //(graph[i].parent == graph[j].parent)
         )
         {
            //distance = sqrt(pow((graph[i].xpos-graph[j].xpos),2.0)
            //               +pow((graph[i].ypos-graph[j].ypos),2.0)
            //               +pow((graph[i].zpos-graph[j].zpos),2.0));
            x0 = graph[i].xpos;
            y0 = graph[i].ypos;
            z0 = graph[i].zpos;
            x1 = graph[j].xpos;
            y1 = graph[j].ypos;
            z1 = graph[j].zpos;
            distance = (x0-x1)*(x0-x1)
                      +(y0-y1)*(y0-y1)
                      +(z0-z1)*(z0-z1);
            //distance = pow((graph[i].xpos-graph[j].xpos),2.0)
            //          +pow((graph[i].ypos-graph[j].ypos),2.0)
            //          +pow((graph[i].zpos-graph[j].zpos),2.0);
            //distance = (graph[i].xpos-graph[j].xpos)*(graph[i].xpos-graph[j].xpos)
             //         +(graph[i].ypos-graph[j].ypos)*(graph[i].ypos-graph[j].ypos)
              //        +(graph[i].zpos-graph[j].zpos)*(graph[i].zpos-graph[j].zpos);
            //if (distance < )
            if (graph[i].zpos < graph[j].zpos)
               boost::add_edge(i,j,distance,graph);
            else 
               boost::add_edge(j,i,distance,graph);
         }
      }
   }
}

template <typename E, typename V>
void
fillmstedmonds(Graph<E,V> &graphin, std::vector<EDesc<E,V>> &vecsamedmonds, Graph<E,V> &graphout)
{
   for (typename std::vector<EDesc<E,V>>::iterator eiter = vecsamedmonds.begin(); eiter != vecsamedmonds.end(); ++eiter)
   {
      boost::add_edge(boost::source(*eiter, graphin),
                      boost::target(*eiter, graphin),
                      boost::get(boost::edge_weight, graphin, *eiter),
                      graphout);

      graphout[boost::source(*eiter, graphin)].track     = graphin[boost::source(*eiter, graphin)].track;
      graphout[boost::source(*eiter, graphin)].parent    = graphin[boost::source(*eiter, graphin)].parent;
      graphout[boost::source(*eiter, graphin)].component = graphin[boost::source(*eiter, graphin)].component;
      graphout[boost::source(*eiter, graphin)].energy    = graphin[boost::source(*eiter, graphin)].energy;
      graphout[boost::source(*eiter, graphin)].xpos      = graphin[boost::source(*eiter, graphin)].xpos;
      graphout[boost::source(*eiter, graphin)].ypos      = graphin[boost::source(*eiter, graphin)].ypos;
      graphout[boost::source(*eiter, graphin)].zpos      = graphin[boost::source(*eiter, graphin)].zpos;

      graphout[boost::target(*eiter, graphin)].track     = graphin[boost::target(*eiter, graphin)].track;
      graphout[boost::target(*eiter, graphin)].parent    = graphin[boost::target(*eiter, graphin)].parent;
      graphout[boost::target(*eiter, graphin)].component = graphin[boost::target(*eiter, graphin)].component;
      graphout[boost::target(*eiter, graphin)].energy    = graphin[boost::target(*eiter, graphin)].energy;
      graphout[boost::target(*eiter, graphin)].xpos      = graphin[boost::target(*eiter, graphin)].xpos;
      graphout[boost::target(*eiter, graphin)].ypos      = graphin[boost::target(*eiter, graphin)].ypos;
      graphout[boost::target(*eiter, graphin)].zpos      = graphin[boost::target(*eiter, graphin)].zpos;
   }
}

template <typename E, typename V>
void
fillmstkruskal(Graph<E,V> &graphin, std::vector<EDesc<E,V>> &vecmstkruskal, Graph<E,V> &graphout)
{
   for (typename std::vector<EDesc<E,V>>::iterator eiter = vecmstkruskal.begin(); eiter != vecmstkruskal.end(); ++eiter)
   {
      boost::add_edge(boost::source(*eiter, graphin),
                      boost::target(*eiter, graphin),
                      boost::get(boost::edge_weight, graphin, *eiter),
                      graphout);

      graphout[boost::source(*eiter, graphin)].track     = graphin[boost::source(*eiter, graphin)].track;
      graphout[boost::source(*eiter, graphin)].parent    = graphin[boost::source(*eiter, graphin)].parent;
      graphout[boost::source(*eiter, graphin)].component = graphin[boost::source(*eiter, graphin)].component;
      graphout[boost::source(*eiter, graphin)].energy    = graphin[boost::source(*eiter, graphin)].energy;
      graphout[boost::source(*eiter, graphin)].xpos      = graphin[boost::source(*eiter, graphin)].xpos;
      graphout[boost::source(*eiter, graphin)].ypos      = graphin[boost::source(*eiter, graphin)].ypos;
      graphout[boost::source(*eiter, graphin)].zpos      = graphin[boost::source(*eiter, graphin)].zpos;

      graphout[boost::target(*eiter, graphin)].track     = graphin[boost::target(*eiter, graphin)].track;
      graphout[boost::target(*eiter, graphin)].parent    = graphin[boost::target(*eiter, graphin)].parent;
      graphout[boost::target(*eiter, graphin)].component = graphin[boost::target(*eiter, graphin)].component;
      graphout[boost::target(*eiter, graphin)].energy    = graphin[boost::target(*eiter, graphin)].energy;
      graphout[boost::target(*eiter, graphin)].xpos      = graphin[boost::target(*eiter, graphin)].xpos;
      graphout[boost::target(*eiter, graphin)].ypos      = graphin[boost::target(*eiter, graphin)].ypos;
      graphout[boost::target(*eiter, graphin)].zpos      = graphin[boost::target(*eiter, graphin)].zpos;
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
         boost::add_edge(*viter,
                         vecmstprim[*viter],
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
exportdata(std::string graphname, Graph<E,V> graph, int eventnumber)
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
