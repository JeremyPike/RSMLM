//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Distance_ANN.h
// Programmer:		Primoz Skraba
// Description:		Metric Information using the ANN library 
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

#ifndef __DISTANCE_ANN__H
#define __DISTANCE_ANN__H

#include <ANN/ANN.h>	

#include "Distance.h"


//--------------------------
// when the vertex class
// has a *double coordinate 
// so we can use ANN
// assumes a coord member
//--------------------------
template <class Iterator>
class  Distance_ANN : public Distance<Iterator>{
 public:

  //----------------------
  // persistent data structure
  // for kd tree
  //----------------------
  ANNkd_tree *kd; 
  int num_points;
  int dim;
  Iterator start;

  
 public:
  //----------------------------------
  //
  //rips parameter
  //
  //----------------------------------
  double mu;
  //----------------------------------
  //----------------------------------


  //----------------------------------
  // initialize data structure
  //----------------------------------
  void initialize(Iterator st,Iterator en,int d){
    int i=0;
    start=st;
    num_points = en-st;
    dim = d;
    //------------------
    //memory allocation
    //------------------
    double** data = new double*[num_points];
    for(;st!=en;st++){
      data[i++] = st->geometry.coord;
    }
    // delete kd; // it gives segfault here 
    kd = new ANNkd_tree(data,num_points,dim);    
  }


  Distance_ANN(){
    kd = NULL;
  }
  ~Distance_ANN(){
    delete kd;
  }
  Distance_ANN(Iterator st,Iterator en,int dim){
    initialize(st,en,dim);
  }

  

  //----------------------------------
  //return squared distance
  //----------------------------------
  inline double distance(Iterator x,Iterator y){
    int i;
    double dist=0;
    for(i=0;i<dim;i++)
      dist += (x->geometry.coord[i] - y->geometry.coord[i])* (x->geometry.coord[i] - y->geometry.coord[i]);
    return dist;
  }

  //----------------------------------
  // get rips neighbors, this is what
  // algorithm calls
  //----------------------------------
  void get_neighbors(Iterator q, std::vector<Iterator> &out){

    int nb_neighb,test,i;
    int *neighb;
    double *ndist;

    out.clear();
    nb_neighb =  kd->annkFRSearch(q->geometry.coord, mu, 0);

    //---------------------------   
    // memory efficient version
    //---------------------------
    // allocate mem each time
    //---------------------------

    neighb = new int[nb_neighb];
    ndist = new double[nb_neighb];

    test =  kd->annkFRSearch(q->geometry.coord, mu, nb_neighb, 
					neighb, ndist);
   
    assert(test==nb_neighb);
    for(i=0;i<nb_neighb;i++){
      Iterator nit =  start+neighb[i];
      if(nit!=q){
	out.push_back(nit);
      }
    }
    delete ndist;
    delete neighb;
  }


  //----------------------------------
  // get number of neighbors 
  // for density estimation
  //----------------------------------
  int get_num_neighbors(Iterator q, double radius){
    int nb_neighb;
    int neighb;
    double ndist;
    nb_neighb =  kd->annkFRSearch(q->geometry.coord, radius*radius, 0, 
				neighb, ndist, 0);
    return nb_neighb;
  }

  //----------------------------------
  // get distances to k nearest neighbors
  // for density estimation
  //----------------------------------
  void get_neighbors_dist(Iterator q,int k, double *ndist){
     int *neighb = new int[k];
     
     kd->annkSearch(q->geometry.coord, k,  neighb, ndist, 0); 
     delete neighb;
  }


  //----------------------------------
  // get distances neighbors within r
  // for density estimation
  //----------------------------------  
  int get_neighbors_dist_r(Iterator q,double radius, double *ndist){
    int nb_neighb,test;
    int *neighb;
    nb_neighb =  kd->annkFRSearch(q->geometry.coord, radius*radius, 0, 
				     neighb, ndist, 0);
    
    neighb = new int[nb_neighb];
    test =  kd->annkFRSearch(q->geometry.coord, radius*radius, nb_neighb, 
			     neighb, ndist, 0);


    delete neighb;

    return test;
  }



};

#endif
