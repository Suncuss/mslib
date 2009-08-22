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


#include "Tree.h"


int main(){

	Tree<double> t;

	double root = 8.3;
	t.setData(root);

	double d1 = 1.1;
	double d2 = 2.1;
	double d3 = 3.1;
	double d4 = 4.1;

	Tree<double> child1; child1.setData(d1);
	Tree<double> child2; child2.setData(d2);
	Tree<double> child3; child3.setData(d3);
	Tree<double> child4; child4.setData(d4);


	child2.addSubTree(&child3);
	child2.addSubTree(&child4);
	child3.setParent(&child2);
	child4.setParent(&child2);

	child1.setParent(&t);
	child2.setParent(&t);
	t.addSubTree(&child1);
	t.addSubTree(&child2);



	cout << "My Max Tree Depth: "<<t.maxDepth(&t)<<endl;
	cout << "My Tree Looks like: "<<endl;
	t.printTree(50);


}
