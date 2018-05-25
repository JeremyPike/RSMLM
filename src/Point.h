//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Point.h
// Programmer:		Primoz Skraba
// Description:		Example Point data struture for use with ANN
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

#ifndef __POINT_H
#define __POINT_H
#include <iostream>
#include "Distance.h"

class ANNPoint{

public:

  typedef struct { 
    bool operator()(const ANNPoint a, const ANNPoint  b) const {
      for (int i=0; i<a.dim; i++) {
	if (a.coord[i] < b.coord[i]) return true;
	else if (a.coord[i] > b.coord[i]) return false;
      }
      return false;
    }
  } Less_Than;



  double *coord;
  int dim;
  

  ANNPoint(){}
  
  ANNPoint(int d){
    dim = d;
  }

  
  //------------------------------------------------------------------------
  // I/O functions
  //------------------------------------------------------------------------
  friend std::istream& operator >>(std::istream &in, ANNPoint &out) {
    int i;
    out.coord = new double[out.dim]; 
    for(i=0;i<out.dim;i++){
      in>>out.coord[i];
    }
    return in;
  }
   

  friend std::ostream& operator<<(std::ostream &out, const ANNPoint &in) {
    int i;
    for(i=0;i<in.dim;i++){
      out<<in.coord[i]<<" ";
    }
    return out;
  }

};






#endif
