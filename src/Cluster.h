//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Cluster.h
// Programmer:		Primoz Skraba
// Description:		Basic Cluster data structure
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

#ifndef __CLUSTER__H
#define __CLUSTER__H

#include <vector>
#include <map>
#include <set>

#include "Interval.h"


//outline for cluster class
template<class Iterator>
class Cluster_Base{
 public:
  void new_cluster(Iterator);
  bool merge(Iterator,Iterator);
  void update(std::vector<Iterator>,Iterator);
  bool test(Iterator);
  bool neighbor_test(Iterator);
  bool valid_neighborhood(std::set<Iterator>,Iterator);
  Iterator gradient(Iterator x,std::set<Iterator>);
 
};


//==============================
//Auxillary cluster information
//==============================
class Cluster_Info{
public:
  bool boundary_flag;
  
};



#endif
