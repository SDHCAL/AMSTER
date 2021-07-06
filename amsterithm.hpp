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
   V numberoftracks = 0;
   VIter<E,V> viter, viterend;
   for (boost::tie(viter, viterend) = boost::vertices(graph); viter != viterend; ++viter)
   {
      if (graph[*viter].track == true)
         ++numberoftracks;
   }
   return numberoftracks;
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
std::vector<E>
countparents(Graph<E,V> graph, V componentnumber)
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
         if (graph[*viter].parent != -211 && graph[*viter].parent != 311)
            ++parents[2];
      }
   }
   return parents;
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
   V nhitscomponent       = numberhitscomponent<E,V>(graph, componentnumber);
   barycentercomponent[0] = barycentercomponent[0]/nhitscomponent;
   barycentercomponent[1] = barycentercomponent[1]/nhitscomponent;
   barycentercomponent[2] = barycentercomponent[2]/nhitscomponent;
   return barycentercomponent;
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
void
updatecomponents(Graph<E,V> &graph)
{
   typename std::vector<V> c(boost::num_vertices(graph));
   boost::connected_components(graph, boost::make_iterator_property_map(c.begin(), get(boost::vertex_index, graph)));
   typename std::vector<V>::iterator vi;
   std::vector<E> barycentercomponent0 = barycentercomponent<E,V>(graph, 0);
   std::vector<E> barycentercomponent1 = barycentercomponent<E,V>(graph, 1);
   E distancetocomponent0;
   E distancetocomponent1;
   for (vi = c.begin(); vi != c.end(); ++vi)
   {
     if (*vi > 1)
     {
         distancetocomponent0 = sqrt(pow((graph[(vi-c.begin())].xpos-barycentercomponent0[0]),2.0)
                                    +pow((graph[(vi-c.begin())].ypos-barycentercomponent0[1]),2.0)
                                    +pow((graph[(vi-c.begin())].zpos-barycentercomponent0[2]),2.0));

         distancetocomponent1 = sqrt(pow((graph[(vi-c.begin())].xpos-barycentercomponent1[0]),2.0)
                                    +pow((graph[(vi-c.begin())].ypos-barycentercomponent1[1]),2.0)
                                    +pow((graph[(vi-c.begin())].zpos-barycentercomponent1[2]),2.0));
        if (distancetocomponent0-distancetocomponent1 <= 0)
           graph[(vi-c.begin())].component = 0;
        else
           graph[(vi-c.begin())].component = 1;
      }
      else
      graph[(vi-c.begin())].component = *vi;
   }
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
   E efficiency = countparents<E,V>(graph,0)[0]/(countparents<E,V>(graph,0)[0]+countparents<E,V>(graph,1)[0]);
   return efficiency;
}

template <typename E, typename V>
E
purity(Graph<E,V> graph)
{
   E purity = countparents<E,V>(graph,0)[0]/numberhitscomponent<E,V>(graph,0);
   return purity;
}

/* clustering algorithm */
template <typename E, typename V>
void
cluster(Graph<E,V> &graph, std::vector<std::pair<EDesc<E,V>,E>> &p, E deltaE)
{
   unsigned int i = 0;
   unsigned int numberofedges = p.size();
   E targetenergy             = graph[0].energy;
   E targetenergyplusdeltaE   = graph[0].energy+deltaE;
   E targetenergyminusdeltaE  = graph[0].energy-deltaE;
   E currentenergy            = hcalenergycomponent<E,V>(graph,0);
   while (currentenergy > targetenergyplusdeltaE && i != numberofedges)
   {
      if (currentenergy > targetenergy)
      {
         boost::remove_edge(boost::source(p[i].first, graph), boost::target(p[i].first, graph), graph);
         updatecomponents<E,V>(graph);
         currentenergy = hcalenergycomponent<E,V>(graph,0);
      }
      if (currentenergy <= targetenergyminusdeltaE)
      {
         boost::add_edge(boost::source(p[i].first, graph), boost::target(p[i].first, graph), p[i].second, graph);
         updatecomponents<E,V>(graph);
         currentenergy = hcalenergycomponent<E,V>(graph,0);
      }
      ++i;
   }
}
