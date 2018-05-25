//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Distance.h
// Programmer:		Primoz Skraba
// Description:	        Metric information data structure 
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

#ifndef __DISTANCE__H
#define __DISTANCE__H

#include <vector>
#include <cassert>


//----------------------------------------------------
// Distance Oracle Concept
//   provide a list of neighbors
//   and return Distance
//----------------------------------------------------

template <class Iterator>
class Distance{
 public:
  //------------------------------------------------------
  //return vector of iterators of neighbors and distance
  //------------------------------------------------------
  void get_neighbors(Iterator, std::vector<Iterator>& ){}
  double distance(Iterator,Iterator){
    return 1;
  }
};

#endif


