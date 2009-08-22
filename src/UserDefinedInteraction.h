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

#ifndef USERDEFINEDINTERACTION_H
#define USERDEFINEDINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "UserDefinedEnergy.h"

using namespace std;

class UserDefinedInteraction: public TwoBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		UserDefinedInteraction();
		UserDefinedInteraction(Atom & _a1, Atom & _a2, string _type);

		// add an operator= 
		UserDefinedInteraction(const UserDefinedInteraction & _interaction);
		~UserDefinedInteraction();




		double getEnergy();
		double getEnergy(double _distance);

		friend ostream & operator<<(ostream &_os, UserDefinedInteraction & _term) {_os << _term.toString(); return _os;};
		string toString() const;

		string getName() const;
		void setName(string _type);
		
	private:
		void setup(Atom * _a1, Atom * _a2, string _type);
		void copy(const UserDefinedInteraction & _interaction);
		double distance;

		string typeName;
		

};

inline double UserDefinedInteraction::getEnergy() {
	return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
}

inline string UserDefinedInteraction::toString() const { char c [1000]; sprintf(c, "USER DEF %s %s %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], distance, energy); return (string)c; };
inline string UserDefinedInteraction::getName() const {return typeName;}
inline void UserDefinedInteraction::setName(string _type) { typeName = _type; }
#endif

