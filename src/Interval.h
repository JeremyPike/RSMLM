//----------------------------------------------------------------------
//----------------------------------------------------------------------
// File:		Interval.h
// Programmer:		Primoz Skraba
// Description:		Interval data structure
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

#ifndef __INTERVAL__H
#define __INTERVAL__H

#include<cassert>
//======================
// basic data structure
// for holding birth and
// death times
//======================


class Interval{
 private: 
  double bi;
  double de;
  bool infinite;
  
 public:
  Interval(){}

  Interval(double b){
    bi = b;
    infinite = true;
  }
  
  void close(double d){
    de = d;
    infinite = false;
  }
  
  double birth(){
    return bi;
  }
  

  double death(){
    assert(!infinite);
    return de;
  }
  
  bool inf(){
    return infinite;
  }
};




#endif
