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

#ifndef HELANAL_H
#define HELANAL_H

#include "CartesianPoint.h"
#include "CartesianGeometry.h"

using namespace std;

/*! \brief  The Helanal object calculates helical axes
 */
class Helanal {



	public:

		Helanal();
		Helanal(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4);
		Helanal(const Helanal & _helanal);
		~Helanal();

		void update();
		void update(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4);
		CartesianPoint & getCA1();
		CartesianPoint & getCA2();
		CartesianPoint & getCA3();
		CartesianPoint & getCA4();

	
		CartesianPoint & getAxis();
		CartesianPoint & getCenter();
		double getTwist() const;
		double getHeight() const;
		double getRadius() const;
		double getResPerTurn() const;
		CartesianPoint & getNpoint();
		CartesianPoint & getCpoint();

		// projection of the CA on the axis
		CartesianPoint & getCA1projection();
		CartesianPoint & getCA2projection();
		CartesianPoint & getCA3projection();
		CartesianPoint & getCA4projection();

		bool fail() const;


	private:
		void setup();
		// input coor
		CartesianPoint * pCA1; 
		CartesianPoint * pCA2; 
		CartesianPoint * pCA3; 
		CartesianPoint * pCA4; 

		// calculated output variables
		CartesianPoint axis;
		CartesianPoint center;
		double twist;
		double height;
		double radius;
		double resPerTurn;
		CartesianPoint Npoint;
		CartesianPoint Cpoint;
		CartesianPoint CA1projection;
		CartesianPoint CA2projection;
		CartesianPoint CA3projection;
		CartesianPoint CA4projection;
		bool errorFlag;


};

#endif
