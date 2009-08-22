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

#include "PrincipleComponentAnalysis.h"
#include "MslExceptions.h"

PrincipleComponentAnalysis::PrincipleComponentAnalysis(){
	compMat.initialize(4,4);
	eigenvectors.clear();
	eigenvalues.clear();
	
	eigenvectors.resize(4);
	eigenvalues.resize(4);
	for (uint i = 0 ; i < 4; i++){
		eigenvectors[i].resize(4);
	}

	printImag = false;
}

PrincipleComponentAnalysis::PrincipleComponentAnalysis(const PrincipleComponentAnalysis &_pca){
	copy(_pca);
}

void PrincipleComponentAnalysis::operator=(const PrincipleComponentAnalysis &_pca){
	copy(_pca);
}

void PrincipleComponentAnalysis::copy(const PrincipleComponentAnalysis &_pca){

	compMat      = _pca.compMat;
	eigenvectors = _pca.eigenvectors;
	eigenvalues  = _pca.eigenvalues;
	printImag    = _pca.printImag;

	toGeoCenter  = _pca.toGeoCenter;
	toOrigin     = _pca.toOrigin;
	
}

PrincipleComponentAnalysis::~PrincipleComponentAnalysis(){
}


void PrincipleComponentAnalysis::computePrincipleComponents(AtomVector &_av){

	// Initialize the component matrix
	for (uint i = 0 ; i < 4; i++){
		for (uint j = 0;j < 4; j++){
			compMat[i][j] = 0.0;
		}
	}

	// Update geometric center, store the vector from Geocenter to Origin

	_av.updateGeometricCenter();
	toGeoCenter = _av.getGeometricCenter();
	toOrigin    = toGeoCenter * -1;

	Transforms t;
	t.translate(_av, toOrigin);



	for (uint i = 0 ; i < _av.size(); i++){

		compMat[0][0] += ( _av(i)[0] *  _av(i)[0]);
		compMat[0][1] += ( _av(i)[0] *  _av(i)[1]);
		compMat[0][2] += ( _av(i)[0] *  _av(i)[2]);

		
		compMat[1][0] += ( _av(i)[1] *  _av(i)[0]);
		compMat[1][1] += ( _av(i)[1] *  _av(i)[1]);
		compMat[1][2] += ( _av(i)[1] *  _av(i)[2]);


		compMat[2][0] += ( _av(i)[2] *  _av(i)[0]);
		compMat[2][1] += ( _av(i)[2] *  _av(i)[1]);
		compMat[2][2] += ( _av(i)[2] *  _av(i)[2]);

	}

	vector<vector<double> > eigens = compMat.getEigenvectorsEigenValuesGSL();

	for (uint i = 0; i < 4; i++){
		for (uint j = 0; j < 4; j++){
			eigenvectors[i][j] = eigens[i][j];
		}
		eigenvalues[i] = eigens[i][4];
	}

	t.translate(_av, toGeoCenter);

}


void PrincipleComponentAnalysis::printPrincpleComponents(){
	

  for (uint i = 0; i < eigenvectors.size();i++){
    for (uint j = 0; j < eigenvectors[i].size();j++){


      if (printImag){
	      char c[40];
	      //sprintf(c, "%-4.6f + %-4.6f i", eigenvectors[i][j].real(),eigenvectors[i][j].imag());
	      sprintf(c, "%-4.6f + NaN i", eigenvectors[i][j]);
	      cout << string(c); 
      } else {
	      char c[20];
	      sprintf(c, "%-4.6f ", eigenvectors[i][j]);
	      cout << string(c); 
      }

    }
    cout << endl;
  }
}


vector<Line> PrincipleComponentAnalysis::getLines(){

	CartesianPoint xDir(eigenvectors[2][0],eigenvectors[2][1],eigenvectors[2][2]);
	CartesianPoint yDir(eigenvectors[1][0],eigenvectors[1][1],eigenvectors[1][2]);
	CartesianPoint zDir(eigenvectors[0][0],eigenvectors[0][1],eigenvectors[0][2]);

	if (xDir.length() == 0 || yDir.length() == 0 || zDir.length() == 0) {
		throw MslGeneralException("PrincipleComponentAnalysis::getLines() : xDir,yDir,zDir has zero length, can not normalize!");
	}


	xDir /= xDir.length();
	yDir /= yDir.length();
	zDir /= zDir.length();
	Line x(toGeoCenter,xDir);
	Line y(toGeoCenter,yDir);
	Line z(toGeoCenter,zDir);

	vector<Line> coordFrame;

	coordFrame.push_back(z);
	coordFrame.push_back(y);
	coordFrame.push_back(x);

	return coordFrame;
	
}



