//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Cluster_Basic.h
// Programmer:		Primoz Skraba
// Description:		Basic Cluster data structure
// Last modified:	Sept 8, 2009 (Version 0.1)
//----------------------------------------------------------------------
//  Copyright (c) 2009 Primoz Skraba.  All Rights Reserved.
//-----------------------------------------------------------------------
//
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
//-----------------------------------------------------------------------
//----------------------------------------------------------------------
// History:
//	Revision 0.1  August 10, 2009
//		Initial release
//----------------------------------------------------------------------
//----------------------------------------------------------------------


#ifndef __CLUSTER__BASIC__H
#define __CLUSTER__BASIC__H

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include "Cluster.h"
#ifndef M_PI
namespace
{
	const double M_PI = std::acos(-1.0);
}
#endif
//---------------------------------------------
// This is heuristic to assign colors 
// for the clusters for outputting to an OFF
// file. Also filters out those clusters whose peak
// is lower than the threshold tau
//---------------------------------------------
template<class Iterator>
int assign_colors(std::map<Iterator, double*>& cluster_colors,
		   std::set<Iterator> &clusters,
		   double tau){

  int nb_clusters = clusters.size();
  int res = 0;

  // reduce nb_clusters to the actual number of clusters 
  for(typename std::set<Iterator>::iterator cit = clusters.begin(); 
      cit!=clusters.end(); cit++)
    // check if cluster does not appear below tau
    if (find_sink(*cit)->func() < tau)
      nb_clusters--;

  //---------------------------------
  // assign distinct colors to clusters (in H,S,V)
  // first, create array of colors;
  //---------------------------------


  std::vector<double*> colors;
  colors.reserve (nb_clusters+1);  // additional entry is for black
  //---------------------------------
  // just choose circularly
  //---------------------------------
  //double ratio =  sqrt(3)/6;
  for (int i = 0; i < nb_clusters; ++i) {
    colors[i] = new double [3];
    double theta = (double)2*i*M_PI/(double)nb_clusters;
    if (theta >= 0 && theta < M_PI/3) {
      colors[i][0] = 1;
      colors[i][1] = 3*theta/M_PI;
      colors[i][2] = 0;
    }
    else if (theta >= M_PI/3 && theta < 2*M_PI/3) {
      colors[i][0] = 1 - 3*(theta-M_PI/3)/M_PI;
      colors[i][1] = 1;
      colors[i][2] = 0;
    }
    else if (theta >= 2*M_PI/3 && theta < M_PI) {
      colors[i][0] = 0;
      colors[i][1] = 1;
      colors[i][2] = 3*(theta-2*M_PI/3)/M_PI;
    }
    else if (theta >= M_PI && theta < 4*M_PI/3) {
      colors[i][0] = 0;
      colors[i][1] = 1 - 3*(theta-M_PI)/M_PI;
      colors[i][2] = 1;
    }
    else if (theta >= 4*M_PI/3 && theta < 5*M_PI/3) {
      colors[i][0] = 3*(theta-4*M_PI/3)/M_PI;
      colors[i][1] = 0;
      colors[i][2] = 1;
    }
    else {  // theta >= 5*M_PI/3 && theta < 2*M_PI
      colors[i][0] = 1;
      colors[i][1] = 0;
      colors[i][2] = 1 - 3*(theta-5*M_PI/3)/M_PI;
    }
  }
  // force black at index nb_clusters
  colors[nb_clusters] = new double [3];
  colors[nb_clusters][0] = 0;
  colors[nb_clusters][1] = 0;
  colors[nb_clusters][2] = 0;

  //---------------------------------
  // compute random permutation
  //---------------------------------
  int *perm = new int[nb_clusters];
  //---------------------------------
  // first, create identity
  //---------------------------------
  for (int i=0; i<nb_clusters; ++i)
    perm [i] = i;
  //---------------------------------
  // then, permute it
  //---------------------------------
  srand((unsigned)time(NULL));
  for (int i=0; i< nb_clusters-1; ++i) {
    int tmp = perm[i];
    int j = i + (int) ((double)rand() * (nb_clusters-i) / RAND_MAX);
    perm[i] = perm[j];
    perm[j] = tmp;
  }

  //---------------------------------
  // now assign color entries to the cluster centers
  //---------------------------------
  int k=0;
  for(typename std::set<Iterator>::iterator cit = clusters.begin(); 
      cit!=clusters.end(); cit++) {
    // check if cluster does not appear below tau
    if (find_sink(*cit)->func() >= tau) {
      cluster_colors[*cit] =  colors[perm[k++]];
      res++;
    }
    else
      cluster_colors[*cit] = colors[nb_clusters];  // black
  }
  
  delete perm;

  return res;
}



//==============================
//basic implementation of 
//cluster class
//==============================
template<class Iterator>
class Cluster : public Cluster_Base<Iterator>{
 private:
  std::vector<Interval> Int_Data;
  //--------------------------------
  // since we are storing it as a 
  // vector, it is equivalent to
  // storing an iterator
  //--------------------------------
  std::map<Iterator, int> Generator;
  //--------------------------------
  //for faster searching
  //--------------------------------


 public:
  double tau;
  
  //----------------------
  //create a new interval
  //----------------------
  void new_cluster(Iterator x){
    Generator[x] = Int_Data.size();
    Int_Data.push_back(Interval(x->func()));

  }
  
  //----------------------
  // Merge two intervals 
  // note this will only output 
  // correctly if persistence 
  // threshold is set to infty
  //----------------------
  bool merge(Iterator x,Iterator y){
    // note: by hypothesis, y->func() >= x->func()
    assert(y->func() >= x->func());

    //---------------------------------
    // test prominences of both clusters
    // assumptions:
    //   - y is its cluster's root
    //   - x is attached to its root directly
    //---------------------------------
    if(std::min(x->get_sink()->func(), y->func()) < x->func() + tau) {
      //---------------------------------
      //kill younger interval
      //---------------------------------
      int i = Generator[x->get_sink()];
      int j = Generator[y];
      if (y->func() <= x->get_sink()->func()) {
	assert(Int_Data[j].inf());
	Int_Data[j].close(x->func());
      }
      else {
	assert(Int_Data[i].inf());
	Int_Data[i].close(x->func());
      }  
      return true;
    }

    return false;
  }

  //----------------------------
  //unwieghted gradient choice
  //----------------------------
  Iterator gradient(Iterator x,std::set<Iterator> &List){
    typename std::set<Iterator>::iterator it;
    Iterator y=x;
    //--------------------
    //find oldest neighbor
    //--------------------
    for(it = List.begin(); it!=List.end();it++){
      if(*it<y)
	y=*it;
    }
    assert(y!=x);
    return y;
  }

  //------------------------------------------------------------
  // weighted gradient choice -- 
  // need access to Distance struct 
  // this simplest way to do this is to 
  // derive from Cluster_Basic and   add a pointer to a 
  // distance structure
  //------------------------------------------------------------
  // Iterator gradient(Iterator x,std::set<Iterator> &List){
  // double dist1 =...;
  // double dist2 =...;
  // return ((x->func()-y1->func())/dist1<(x->func()-y2->func())/dist2) ? y1 : y2; 
  //}


  //-----------------------------------
  // This is a user-defined test
  // if you want to skip over
  // some nodes
  //-----------------------------------
  inline bool test(Iterator){
    return true;
  }
  
  //-----------------------------------
  // Same thing but for valid
  // neighbors
  //-----------------------------------
  bool neighbor_test(Iterator x){
    return !x->data.boundary_flag;
  } 
  
  //-----------------------------------
  // check to make sure this is a true
  // peak
  //----------------------------------- 
  bool valid_neighborhood(std::set<Iterator> &in){
    typename std::set<Iterator>::iterator it;
    if(in.size()==0)
      return true;

    for(it = in.begin(); it!=in.end(); it++){
      if(!((*it)->data.boundary_flag))
	return true;
    }
    return false;
  } 


  //------------------------------
  // Output intervals 
  //------------------------------
  void output_intervals(std::ofstream &out){
    std::vector<Interval>::iterator it;
    for(it = Int_Data.begin();it != Int_Data.end();it++){
      out << it->birth() << " ";
      (it->inf()) ? out<<"-inf" : out<<it->death();
      out << std::endl;
    }
  }
  // added
  void output_intervals(std::vector<double> &outB, std::vector<double> &outD){
    std::vector<Interval>::iterator it;
   
    for(it = Int_Data.begin();it != Int_Data.end();it++){

      outB.push_back(it->birth());
      
      if(it->inf())
        outD.push_back(-1);
      else
        outD.push_back(it->death());
    }
  }

  //------------------------------
  // Ouput points as coordinates
  // with an assignment of clusterheads
  // note: input array is an array
  // of iterators on the real point
  // cloud array (to cope with
  // permutation from initial sort)
  //------------------------------
  template <class IIterator>
  void output_clusters(std::ofstream &out, IIterator start, IIterator finish){
    //------------------------------
    //run through and 
    //create map from prominent peaks 
    //------------------------------
    std::map<Iterator,int> cluster_ids;

    cluster_ids.clear();

    int i=1;
    //------------
    //number nodes
    //------------
    for(IIterator it=start;it!=finish;it++){
      // assert(find_sink(*it) == (*it)->get_sink());
      typename std::map<Iterator,int>::iterator mit = 
	cluster_ids.find(find_sink(*it));
      if(mit == cluster_ids.end()) {
	//--------------
	// check peak height
	//--------------
	if (find_sink(*it)->func() >= tau)
	  cluster_ids[find_sink(*it)]=i++;
      }
    }

    for(IIterator it=start;it!=finish;it++) {
	//--------------
	// check peak height
        // if it is ok, then 
        // output cluster number
        // otherwise output NaN
	//--------------
      // out<<(*it)->geometry<<" ";
      if (find_sink(*it)->func() >= tau)
	out<<cluster_ids[find_sink(*it)];
      else
	out<<"NaN";
      out<<std::endl;
    }
  }  

  //------------------------------
  // Added
  //------------------------------
  template <class IIterator>
  void output_clusters(std::vector<int> &out, IIterator start, IIterator finish) {
	  //------------------------------
	  //run through and 
	  //create map from prominent peaks 
	  //------------------------------
	  std::map<Iterator, int> cluster_ids;

	  cluster_ids.clear();

	  int i = 1;
	  //------------
	  //number nodes
	  //------------
	  for (IIterator it = start; it != finish; it++) {
		  // assert(find_sink(*it) == (*it)->get_sink());
		  typename std::map<Iterator, int>::iterator mit =
			  cluster_ids.find(find_sink(*it));
		  if (mit == cluster_ids.end()) {
			  //--------------
			  // check peak height
			  //--------------
			  if (find_sink(*it)->func() >= tau)
				  cluster_ids[find_sink(*it)] = i++;
		  }
	  }

	  for (IIterator it = start; it != finish; it++) {
		  //--------------
		  // check peak height
		  // if it is ok, then 
		  // output cluster number
		  // otherwise output -1
		  //--------------
		  // out<<(*it)->geometry<<" ";
		  if (find_sink(*it)->func() >= tau) {
	
			  out.push_back(cluster_ids[find_sink(*it)]);
			 
		  } 
		  else
			  out.push_back(-1);
		  
	  }
  }


  //------------------------------
  // Ouput COFF file
  // works only in 2 and 3 dim
  // and outputs height as function
  // value in 2 dim
  //------------------------------
  void output_clusters_coff(std::ofstream &out,Iterator start,Iterator finish){
    //run through and 
    //create map from prominent peaks
    std::map<Iterator,double*> colors;
    std::set<Iterator> peaks;
    
    typename std::set<Iterator>::iterator lit;
  
  
    colors.clear();

    int num_nodes=0;
    //------------
    //number nodes
    //------------
    for(Iterator it=start;it!=finish;it++){
      peaks.insert(find_sink(it));
      num_nodes++;
    }

    //--------------------------------
    //create color scheme
    //--------------------------------
    std::cout << "Number of clusters: "
	      << assign_colors(colors, peaks, tau) << std::endl;
    out<<"COFF"<<std::endl<<num_nodes<<" "<<num_nodes<<" 0"<<std::endl;

    for(Iterator it=start;it!=finish;it++){
      // truncate to first 3 coordinates
      for (int i=0; i<min(3, it->geometry.dim); i++)
	     out << it->geometry.coord[i] << " ";
      // add additional 0 coordinates if necessary
      if(it->geometry.dim<=0)
	out << "0 ";
      if(it->geometry.dim<=1)
	out << "0 ";
      if(it->geometry.dim<=2)
	out << "0 ";
      
      Iterator sink = find_sink(it);
      out << colors[sink][0] << " "
	  << colors[sink][1] << " "
	  << colors[sink][2] << " "
	  << "0" << std::endl;
    }
    for(int i=0;i<num_nodes;i++){
      out << "1 " << i << std::endl;
    }

  }

};


#endif
