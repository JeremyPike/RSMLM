
#include <iostream>
#include <sstream>
#include <cassert>

#include <algorithm>

#include <Rcpp.h>
#include "Vertex.h"
#include "Distance_ANN.h"
#include "Core.h"
#include "Cluster_Basic.h"
#include "Point.h"
#include "Density.h"	
using namespace Rcpp;

using namespace std;


//rename for brevity
typedef Vertex<ANNPoint,Cluster_Info > Point;

// comparison function object for vector indices
template<class V> class Less_Than {
protected:
  V& v;
public:
  Less_Than (V& v_): v(v_){}
  bool operator()(const int a, const int  b) const 
  {return Point::Less_Than()(v[a], v[b]);}
};



// [[Rcpp::export]]
List tomato(NumericMatrix coords, int num_neighb, double r, double threshold) {
  
  
  int nb_points = coords.nrow();
  int dim = coords.ncol();
  
  vector< Point > point_cloud;
  
  for (int i=0; i<nb_points; i++) {
    ANNPoint p(dim);
    p.coord = new double[dim];
    for (int j=0; j<dim; j++) 
      p.coord[j] = coords(i, j);
    
    Point v(p);
    v.data.boundary_flag=false;
    point_cloud.push_back(v);
  }
  
  
  //create distance structure
  Distance_ANN< vector< Point >::iterator > metric_information;
  metric_information.initialize(point_cloud.begin(),
                                point_cloud.end(),
                                dim);
  
  
  
  //compute density

  
  distance_to_density(point_cloud.begin(),point_cloud.end(),
                      num_neighb, metric_information);
  
  
  // sort point cloud and retrieve permutation (for pretty output)
  vector<int> perm;
  perm.reserve(nb_points);
  for(int i=0; i < nb_points; i++)
    perm.push_back(i);
  std::sort(perm.begin(), perm.end(), Less_Than<vector<Point> >(point_cloud));
  // store inverse permutation as array of iterators on initial point cloud
  vector< vector<Point>::iterator> pperm;
  pperm.reserve(nb_points);
  for (int i=0; i<nb_points; i++)
    pperm.push_back(point_cloud.begin());
  for (int i=0; i<nb_points; i++)
    pperm[perm[i]] = (point_cloud.begin() + i);
  // operate permutation on initial point cloud 
  vector<Point> pc;
  pc.reserve(nb_points);
  for (int i=0; i<nb_points; i++)
    pc.push_back(point_cloud[i]);
  for (int i=0; i<nb_points; i++)
    point_cloud[i] = pc[perm[i]];
  
  //update distance structure --- since it relies on the order of entry
  metric_information.initialize(point_cloud.begin(),point_cloud.end(), dim);
  
  //set rips parameter

  metric_information.mu = r*r;
  
  
  //create cluster data structure
  Cluster< vector< Point >::iterator > output_clusters;
  //set threshold
  output_clusters.tau = threshold;

  
  // perform clustering
  compute_persistence(point_cloud.begin(),point_cloud.end(),
                      metric_information,output_clusters);
  
  
  // compress data structure:
  // attach each data point to its cluster's root directly
  // to speed up output processing
  attach_to_clusterheads(point_cloud.begin(),point_cloud.end());
  
  
  vector<int> clusters;
  output_clusters.output_clusters(clusters, pperm.begin(), pperm.end());
  vector<double> births;
  vector<double> deaths;
  output_clusters.output_intervals(births, deaths);
  
   NumericMatrix diagram(births.size(), 2);
   for(int i=0; i < births.size(); i++) {
     diagram(i, 0) = births.at(i);
     diagram(i, 1) = deaths.at(i);
   }
 


  return List::create(Named("clusters") = clusters,
                      Named("diagram") = diagram);
}

// [[Rcpp::export]]
List tomatoDens(NumericMatrix coords, NumericVector density, double r, double threshold) {
  
  
  int nb_points = coords.nrow();
  int dim = coords.ncol();
  
  vector< Point > point_cloud;
  
  for (int i=0; i<nb_points; i++) {
    ANNPoint p(dim);
    p.coord = new double[dim];
    for (int j=0; j<dim; j++) 
      p.coord[j] = coords(i, j);
    
    Point v(p);
    v.set_func(density[i]);
      v.data.boundary_flag=false;
    point_cloud.push_back(v);
  }
  
  
  //create distance structure
  Distance_ANN< vector< Point >::iterator > metric_information;
  // this one had to be commented out since initialising the metric_information twice causes a memory leak (plus, it's also useless doing it at this point)
  // metric_information.initialize(point_cloud.begin(),
  //                               point_cloud.end(),
  //                               dim);
  
  
  
  
  // sort point cloud and retrieve permutation (for pretty output)
  vector<int> perm;
  perm.reserve(nb_points);
  for(int i=0; i < nb_points; i++)
    perm.push_back(i);
  std::sort(perm.begin(), perm.end(), Less_Than<vector<Point> >(point_cloud));
  // store inverse permutation as array of iterators on initial point cloud
  vector< vector<Point>::iterator> pperm;
  pperm.reserve(nb_points);
  for (int i=0; i<nb_points; i++)
    pperm.push_back(point_cloud.begin());
  for (int i=0; i<nb_points; i++)
    pperm[perm[i]] = (point_cloud.begin() + i);
  // operate permutation on initial point cloud 
  vector<Point> pc;
  pc.reserve(nb_points);
  for (int i=0; i<nb_points; i++)
    pc.push_back(point_cloud[i]);
  for (int i=0; i<nb_points; i++)
    point_cloud[i] = pc[perm[i]];
  
  //update distance structure --- since it relies on the order of entry
  metric_information.initialize(point_cloud.begin(),point_cloud.end(), dim); // <<-- here's the problem
  
  //set rips parameter
  
  metric_information.mu = r*r;
  
  
  //create cluster data structure
  Cluster< vector< Point >::iterator > output_clusters;
  //set threshold
  output_clusters.tau = threshold;
  
  
  // perform clustering
  compute_persistence(point_cloud.begin(),point_cloud.end(),
                      metric_information,output_clusters);
  
  
  // compress data structure:
  // attach each data point to its cluster's root directly
  // to speed up output processing
  attach_to_clusterheads(point_cloud.begin(),point_cloud.end());
  
  
  vector<int> clusters;
  output_clusters.output_clusters(clusters, pperm.begin(), pperm.end());
  vector<double> births;
  vector<double> deaths;
  output_clusters.output_intervals(births, deaths);
  
  NumericMatrix diagram(births.size(), 2);
  for(int i=0; i < births.size(); i++) {
    diagram(i, 0) = births.at(i);
    diagram(i, 1) = deaths.at(i);
  } 
  
  return List::create(Named("clusters") = clusters,
                      Named("diagram") = diagram);
}

    
    
