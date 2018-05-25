//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Density.h
// Programmer:		Primoz Skraba
// Description:		Templated methods for density estimation
// Last modified:	August 10, 2009 (Version 0.1)
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
#ifndef __DENSITY_H
#define __DENSITY_H

#include <cmath>
#include <cassert>
#include <vector>


//------------------------------
// simple ball density estimator
// iterate over all
//------------------------------
template <class Distance, class Iterator>
  void ball_density(Iterator start, Iterator finish, double radius, int nb_points, Distance &query_struct)
{
  Iterator it;
  std::vector<Iterator> out;
  for(it=start;it!=finish;it++){
    int nb_neighb  = query_struct.get_num_neighbors(it,radius);  
    it->set_func(-((double)nb_neighb+1)/ (double)nb_points); 
  }  

}


//------------------------------
// distance to measure
//------------------------------
template <class Distance, class Iterator>
  void distance_to_density(Iterator start, Iterator finish, int k, Distance &query_struct)
{
  Iterator it;
  double sum;
  double *ndist = new double[k];

  for(it=start;it!=finish;it++){
    query_struct.get_neighbors_dist(it,k,ndist);  
    sum=0;
    for (int i=0; i<k; ++i) {
      sum += ndist[i]*ndist[i];
    }
    it->set_func(sqrt((double)k/(double)sum));
    /* it->set_func(sqrt(sum/(double)k)+1); */
  }  

  delete ndist;
}


//------------------------------
// Gaussian kernel 
// bounded by num of nearest
//  neighbors
//------------------------------
template <class Distance, class Iterator>
  void gaussian_NN(Iterator start, Iterator finish, int k,double h, Distance &query_struct)
{
  Iterator it;
  double sum;
  double *ndist = new double[k];

  for(it=start;it!=finish;it++){
    query_struct.get_neighbors_dist(it,k,ndist);  
    sum=0;

    for (int i=0; i<k; ++i) {
      sum +=  exp(ndist[i]*(-0.5)/h);
    }
    
    it->set_func(-sum);
  }  

  delete ndist;
}


//------------------------------
// Gaussian kernel 
// bounded by cuttoff point
//------------------------------
template <class Distance, class Iterator>
  void gaussian_mu(Iterator start, Iterator finish, double radius, 
		   double h,int nb_points, Distance &query_struct)
{
  Iterator it;
  double sum;
  double *ndist = new double[nb_points];

  for(it=start;it!=finish;it++){
    int nb_neighb  = query_struct.get_neighbors_dist_r(it,radius,ndist);  
    sum=0;
    for (int i=0; i<nb_neighb; ++i) {
      sum +=  exp(ndist[i]*(-0.5)/h);
    }
    
    it->set_func(-sum);
  }  

  delete ndist;
}

#endif
