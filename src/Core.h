//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Core.h
// Programmer:		Primoz Skraba
// Description:		Core method for clustering
// Last modified:       Sept. 8, 2009 (Version 0.1)
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

#ifndef __CORE__H
#define __CORE__H

#include <iostream>
#include <vector>
#include <cassert>
#include <set> 

//--------------------------
// Vertex class should 
// have < defined
// have - defined -> for tau
// singly-linked list
//    sink() ->return sink
//    sink(itself)->set_sink
//--------------------------


//-----------------------------
// Follow linked list to the end
//-----------------------------
template<class Iterator>
Iterator find_sink(Iterator in){
  assert(in != Iterator());
  while(in->get_sink() != in){
    in = in->get_sink();
    assert(in != Iterator());
  }
  return in;
}


//-----------------------------
// collapse union-find data structure
//-----------------------------
template <class Iterator>
void attach_to_clusterheads(Iterator start, Iterator end){
  for(Iterator it=start;it!=end;it++){
    it->set_sink(find_sink(it));
  }
   
}


//-----------------------------
// main algorithm
//-----------------------------
template <class Iterator, class Distance, class Cluster>
  void compute_persistence(Iterator start, 
			   Iterator end, 
			   Distance &Distance_Oracle,
			   Cluster &Cluster_Data_Struct){
  


  Iterator gradient;
  Iterator sink;
  std::vector<Iterator> adjacent_nodes;
  typename std::vector<Iterator>::iterator neighb;

  
  //=================================
  // assume vertex container is 
  // presorted by function value
  //=================================
  for(Iterator vit=start; vit != end; vit++){
        
    //--------------------------------
    // clear adjacency list
    //--------------------------------
    adjacent_nodes.clear();

    //-----------------------------------
    // get neighbors
    //-----------------------------------
    Distance_Oracle.get_neighbors(vit,adjacent_nodes);
   
    //-----------------------------------
    // find gradient, if it exists
    //-----------------------------------
    gradient = vit;

    for(neighb = adjacent_nodes.begin();
	neighb != adjacent_nodes.end();
	neighb++){
      assert(*neighb != vit);
      if (*neighb < gradient)
	gradient = *neighb;
    }
    
    
    //------------------------------
    // if no gradient, then declare
    // vit a peak
    //------------------------------
    if(gradient == vit){
      //-----------------------
      //set the sink to itself
      //-----------------------
      vit->set_sink(vit);

      //-----------------------
      // check to make sure it
      // worked 
      //-----------------------
      assert(find_sink(vit)==vit);

      //-----------------------
      // create new cluster
      //-----------------------
      Cluster_Data_Struct.new_cluster(vit);
    }

    else{
      //----------------------------------
      // gradient has been found
      // attach vit to gradient's cluster
      //----------------------------------
	
      vit->set_sink(find_sink(gradient));

      //----------------------------------
      // check that vit is right
      // below the cluster root
      // (important invariant in the 
      // following)
      //----------------------------------
      assert(vit->get_sink()==find_sink(vit));
      
      //------------------------------
      // Go through the neighbors again
      // to see if their clusters 
      // can be merged with vit's
      //------------------------------
      for(neighb = adjacent_nodes.begin();
	  neighb != adjacent_nodes.end();
	  neighb++){
	
	//----------------------
	// only consider older
	// neighbors
	//----------------------
	if(vit<*neighb)
	  continue;

	//----------------------
	// find sink of neighbor
	//----------------------
	sink = find_sink(*neighb);

	//----------------------------------
	// check that vit is still right
	// below its cluster's root
	// (cf invariant)
	//----------------------------------
	assert(vit->get_sink()==find_sink(vit));

	//----------------------------
	// no need to do anything
	// if neighb's sink is the 
	// same as vit's
	//----------------------------
	if (sink==vit->get_sink())
	  continue;

	//----------------------------
	//check if merge conditions
	// are met
	//----------------------------
	if(Cluster_Data_Struct.merge(vit,sink)){

	  //-------------------------
	  // If so, then merge cluster 
	  // with lower peak into 
	  // cluster with higher peak
	  //-------------------------
	  if (vit->get_sink()->func() >= sink->func())
	    sink->set_sink(vit->get_sink());
	  else {
	    vit->get_sink()->set_sink(sink);

	    //-------------------------
	    // compress path from vit 
	    // to root on the fly
	    // (cf invariant)
	    //-------------------------
	    vit->set_sink(sink);
	  }
	}
      }

    }
  }

}



#endif






