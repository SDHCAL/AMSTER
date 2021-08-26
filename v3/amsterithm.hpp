/* 
 * functions for the clustering algorithm
 */

template <typename E, typename V>
V
numbercomponents(Graph<E,V> graph)
{
   std::vector<V> c(boost::num_vertices(graph));
   return boost::connected_components(graph, boost::make_iterator_property_map(c.begin(), get(boost::vertex_index, graph)));
}

template <typename E, typename V>
V
numbertracks(Graph<E,V> graph)
{
   V numbertracks = 0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].track == true)
         ++numbertracks;
   }
   return numbertracks;
}

template <typename E, typename V>
V
numberhitscomponent(Graph<E,V> graph, V componentnumber)
{
   V nhits = 0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].component == componentnumber)
         ++nhits;
   }
   return nhits;
}

template <typename E, typename V>
E
numberparents(Graph<E,V> graph, V parentnumber)
{
   E numberparents = 0.0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].parent == parentnumber)
         ++numberparents;
   }
   return numberparents;
}
/*
template <typename E, typename V>
std::vector<E>
numberparentscomponent(Graph<E,V> graph, V componentnumber)
{
   std::vector<E> parents(3,0.0);
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].component == componentnumber)
      {
         if (graph[*viter].parent == -211)
            ++parents[0];
         if (graph[*viter].parent == 311)
            ++parents[1];
         else if (graph[*viter].parent != -211 && graph[*viter].parent != 311)
            ++parents[2];
      }
   }
   return parents;
}
*/
template <typename E, typename V>
E
numberparentscomponent(Graph<E,V> graph, V parentnumber, V componentnumber)
{
   E numberparentscomponent = 0.0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].parent == parentnumber && graph[*viter].component == componentnumber)
         ++numberparentscomponent;
   }
   return numberparentscomponent;
}

template <typename E, typename V>
void
sortedges(Graph<E,V> graph, std::vector<std::pair<EDesc<E,V>,E>> &P)
{
   EIter<E,V> eiter, eiterend;
   for (boost::tie(eiter, eiterend) = boost::edges(graph); eiter != eiterend; ++eiter)
   {
        P.push_back(std::make_pair(*eiter, boost::get(boost::edge_weight, graph, *eiter)));
   }
   std::sort
   (
      P.begin(), P.end(), [](auto &left, auto &right)
      {
        return right.second < left.second;
      }
   );
}
/*
template <typename E, typename V>
void
updatecomponents(Graph<E,V> &graph)
{
   typename std::vector<V> c(boost::num_vertices(graph));
   boost::connected_components(graph, boost::make_iterator_property_map(c.begin(), get(boost::vertex_index, graph)));
   typename std::vector<V>::iterator vi;
   for (vi = c.begin(); vi != c.end(); ++vi)
   {
     if (*vi > 1)
         graph[(vi-c.begin())].component = 1;
     else
         graph[(vi-c.begin())].component = *vi;
   }
}
*/
template <typename E, typename V>
std::array<E,3>* 
barycenterComponents(Graph<E,V> &graph, V N_components)
{
   std::array<E,3>* tab = new std::array<E,3>[N_components];
   for (int i=0; i<N_components; ++i)
      tab[i].fill(0);
   V N[N_components] = {};
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      tab[graph[*viter].component][0]+= graph[*viter].xpos;
      tab[graph[*viter].component][1]+= graph[*viter].ypos;
      tab[graph[*viter].component][2]+= graph[*viter].zpos;
      ++N[graph[*viter].component];
   }
   for (int i=0; i<N_components; ++i)
      for (int j=0; j<3; ++j)
         tab[i][j]/=N[i];
   return tab; //beware memory leak
}

template <typename E, typename V>
std::vector<E>
barycentercomponent(Graph<E,V> &graph, V componentnumber)
{
   std::vector<E> barycentercomponent(3,0.0);
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].component == componentnumber)
      {
          barycentercomponent[0] += graph[*viter].xpos;
          barycentercomponent[1] += graph[*viter].ypos;
          barycentercomponent[2] += graph[*viter].zpos;
      }
   }
   E nhitscomponent       = numberhitscomponent<E,V>(graph, componentnumber);
   barycentercomponent[0] = barycentercomponent[0]/nhitscomponent;
   barycentercomponent[1] = barycentercomponent[1]/nhitscomponent;
   barycentercomponent[2] = barycentercomponent[2]/nhitscomponent;
   return barycentercomponent;
}

template <typename E, typename V>
void
updatecomponents(Graph<E,V> &graph)
{
   typename std::vector<V> c(boost::num_vertices(graph));
   boost::connected_components(graph, boost::make_iterator_property_map(c.begin(), get(boost::vertex_index, graph)));
   typename std::vector<V>::iterator vi;
   std::array<E,3>* bar = barycenterComponents(graph, numbercomponents<E,V>(graph));
   std::array<E,3> barycentercomponent0 = bar[0];
   std::array<E,3> barycentercomponent1 = bar[1];
   //std::vector<E> barycentercomponent0 = barycentercomponent<E,V>(graph, 0);
   //std::vector<E> barycentercomponent1 = barycentercomponent<E,V>(graph, 1);
   E distancetocomponent0square;
   E distancetocomponent1square;
   for (vi = c.begin(); vi != c.end(); ++vi)
   {
     if (*vi > 1)
     {
         E x0 = graph[(vi-c.begin())].xpos-barycentercomponent0[0];
         E y0 = graph[(vi-c.begin())].ypos-barycentercomponent0[1];
         E z0 = graph[(vi-c.begin())].zpos-barycentercomponent0[2];
         E x1 = graph[(vi-c.begin())].xpos-barycentercomponent1[0];
         E y1 = graph[(vi-c.begin())].ypos-barycentercomponent1[1];
         E z1 = graph[(vi-c.begin())].zpos-barycentercomponent1[2];
         distancetocomponent0square = x0*x0+y0*y0+z0*z0; 
         distancetocomponent1square = x1*x1+y1*y1+z1*z1;
        if (distancetocomponent0square-distancetocomponent1square <= 0)
           graph[(vi-c.begin())].component = 0;
        else
           graph[(vi-c.begin())].component = 1;
      }
      else
      graph[(vi-c.begin())].component = *vi;
   }
   delete [] bar;
}

template <typename E, typename V>
E
hcalenergy(Graph<E,V> &graph)
{
   V nhit1 = 0;
   V nhit2 = 0;
   V nhit3 = 0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      /* should add track == false if track = 1, 2 or 3 GeV to avoid counting it */
      //if (graph[*viter].track == false)
      {
         if      (graph[*viter].energy == 1)
            ++nhit1;
         else if (graph[*viter].energy == 2)
            ++nhit2;
         else if (graph[*viter].energy == 3)
            ++nhit3;
      }
   }
   V nhit       = nhit1       + nhit2       + nhit3;
   E alpha      = alpha0      + alpha1*nhit + alpha2*nhit*nhit;
   E beta       = beta0       + beta1*nhit  + beta2*nhit*nhit;
   E gamma      = gamma0      + gamma1*nhit + gamma2*nhit*nhit;
   E hcalenergy = alpha*nhit1 + beta*nhit2  + gamma*nhit3;
   return hcalenergy;
}

template <typename E, typename V>
E
hcalenergyparent(Graph<E,V> &graph, V parent)
{
   V nhit1 = 0;
   V nhit2 = 0;
   V nhit3 = 0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].parent == parent)
      {
         if      (graph[*viter].energy == 1)
            ++nhit1;
         else if (graph[*viter].energy == 2)
            ++nhit2;
         else if (graph[*viter].energy == 3)
            ++nhit3;
      } 
   }
   V nhit       = nhit1       + nhit2       + nhit3;
   E alpha      = alpha0      + alpha1*nhit + alpha2*nhit*nhit;
   E beta       = beta0       + beta1*nhit  + beta2*nhit*nhit;
   E gamma      = gamma0      + gamma1*nhit + gamma2*nhit*nhit;
   E hcalenergyparent = alpha*nhit1 + beta*nhit2  + gamma*nhit3;
   return hcalenergyparent;
}

template <typename E, typename V>
E
hcalenergycomponent(Graph<E,V> &graph, V componentnumber)
{
   V nhit1 = 0;
   V nhit2 = 0;
   V nhit3 = 0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      /* should add track == false if track = 1, 2 or 3 GeV to avoid counting it */
      //if (graph[*viter].track == false && graph[*viter].component == componentnumber)
      if (graph[*viter].component == componentnumber)
      {
         if      (graph[*viter].energy == 1)
            ++nhit1;
         else if (graph[*viter].energy == 2)
            ++nhit2;
         else if (graph[*viter].energy == 3)
            ++nhit3;
      } 
   }
   V nhit       = nhit1       + nhit2       + nhit3;
   E alpha      = alpha0      + alpha1*nhit + alpha2*nhit*nhit;
   E beta       = beta0       + beta1*nhit  + beta2*nhit*nhit;
   E gamma      = gamma0      + gamma1*nhit + gamma2*nhit*nhit;
   E hcalenergycomponent = alpha*nhit1 + beta*nhit2  + gamma*nhit3;
   return hcalenergycomponent;
}

template <typename E, typename V>
E
efficiency(Graph<E,V> graph)
{
   E efficiency = numberparentscomponent<E,V>(graph,-211,0)/numberparents<E,V>(graph,-211);
   return efficiency;
}

template <typename E, typename V>
E
purity(Graph<E,V> graph)
{
   E purity = numberparentscomponent<E,V>(graph,-211,0)/numberhitscomponent<E,V>(graph,0);
   return purity;
}

template <typename E, typename V>
E
minimumpurity(Graph<E,V> graph)
{
   E minimumpurity = numberparents<E,V>(graph,-211)/(numberparents<E,V>(graph,-211)+numberparents<E,V>(graph,311));
   return minimumpurity;
}

/* clustering algorithm */
template <typename E, typename V>
void
cluster(Graph<E,V> &graph, std::vector<std::pair<EDesc<E,V>,E>> &p, E deltaE)
{
   int i = 0;
   int numberedges = boost::num_edges(graph);
   E targetenergyplusdeltaE   = graph[0].energy+deltaE;
   E targetenergyminusdeltaE  = graph[0].energy-deltaE;
   E currentenergy            = hcalenergycomponent<E,V>(graph,0);
   while (currentenergy > targetenergyplusdeltaE && i != numberedges)
   {
         boost::remove_edge(boost::source(p[i].first, graph), boost::target(p[i].first, graph), graph);
         updatecomponents<E,V>(graph);
         currentenergy = hcalenergycomponent<E,V>(graph,0);
      if (currentenergy <= targetenergyminusdeltaE)
      {
         boost::add_edge(boost::source(p[i].first, graph), boost::target(p[i].first, graph), p[i].second, graph);
         updatecomponents<E,V>(graph);
         currentenergy = hcalenergycomponent<E,V>(graph,0);
      }
      ++i;
   }
}
/* ordered edges with descending energy removed in edge removed */
