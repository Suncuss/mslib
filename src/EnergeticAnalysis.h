/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/

#ifndef ENERGETICANALYSIS
#define ENERGETICANALYSIS

#include "System.h"
#include "MslTools.h"
#include "PyMolVisualization.h"

#include <map>


class EnergeticAnalysis {

	public:
		EnergeticAnalysis();
		~EnergeticAnalysis();

		void analyzePosition(System &_sys, int _position, int _rotamer=0);

		string toString();
		friend ostream & operator<<(ostream &_os, EnergeticAnalysis & _eaObj)  {_os << _eaObj.toString(); return _os;};

	private:


};
#endif