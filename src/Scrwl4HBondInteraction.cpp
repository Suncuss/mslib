/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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


#include "Scrwl4HBondInteraction.h"
#include <math.h>

using namespace MSL;
using namespace std;
const string Scrwl4HBondInteraction::typeName = "SCRWL4_HBOND";


// parameters from "G.G.Krivov et al,Improved prediction of protein side-chain conformations with SCRWL4"
const double Scrwl4HBondInteraction::d0 = 2.08;
const double Scrwl4HBondInteraction::sig_d = 0.67;
const double Scrwl4HBondInteraction::B = 35.0;
const double Scrwl4HBondInteraction::cos_alpha_max = cos((37.0 * M_PI)/180.0); //  37 degrees
const double Scrwl4HBondInteraction::cos_beta_max = cos((49.0 *M_PI)/180.0) ; // 49 degrees
const double Scrwl4HBondInteraction::denominator = sig_d * sqrt((1-cos_alpha_max) * (1-cos_beta_max)); 
bool Scrwl4HBondInteraction::debugFlagOn = false;


Scrwl4HBondInteraction::Scrwl4HBondInteraction() {
	setup(NULL, NULL, NULL, NULL, NULL,0.0, 0.0, 0.0, 0.0,1.0);
}

Scrwl4HBondInteraction::Scrwl4HBondInteraction(Atom & _d1, Atom & _d2,Atom & _a1, Atom & _a2,Atom& _a3,double _dist, double _ang,double _e1_dihe,double _e2_dihe,double _scalingFactor) {
	setup (&_d1, &_d2,&_a1, &_a2,&_a3,_dist, _ang,_e1_dihe,_e2_dihe,_scalingFactor);
}

Scrwl4HBondInteraction::Scrwl4HBondInteraction(const Scrwl4HBondInteraction & _interaction) {
	setup(NULL, NULL, NULL, NULL, NULL,0.0, 0.0, 0.0, 0.0,1.0);
	copy(_interaction);
}

Scrwl4HBondInteraction::~Scrwl4HBondInteraction() {
}


void Scrwl4HBondInteraction::setup(Atom * _pA1, Atom * _pA2, Atom * _pA3, Atom * _pA4, Atom * _pA5, double _dist, double _ang, double _e1_dihe, double _e2_dihe,double _scalingFactor) {
	pAtoms.push_back(_pA1); // actual donor 
	pAtoms.push_back(_pA2);  // bonded to donor
	pAtoms.push_back(_pA3);  // actual acceptor
	pAtoms.push_back(_pA4);  // bonded to acceptor
	pAtoms.push_back(_pA5);  // bonded to acceptor_2
	params.push_back(_dist);
	params.push_back(_ang);
	params.push_back(_e1_dihe);
	params.push_back(_e2_dihe);
	scalingFactor = _scalingFactor;
}

void Scrwl4HBondInteraction::copy(const Scrwl4HBondInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	scalingFactor = _interaction.scalingFactor;
}
double Scrwl4HBondInteraction::getEnergy(double _d) {
		cerr << "WARNING 23567: This function is not implemented" << endl;
		return 100.0;
}

void Scrwl4HBondInteraction::printParameters() {
	cout << " d0: " << d0 << endl;
	cout << " sig_d: " << sig_d << endl;
	cout << " B: " << B << endl;
	cout << " denominator: " << denominator << endl;
	cout << " cos_alpha_max " << cos_alpha_max << endl;
	cout << " cos_beta_max " << cos_beta_max << endl;
}

void Scrwl4HBondInteraction::setDebugFlagOn(bool _debugFlagOn) {
	debugFlagOn = _debugFlagOn;
}

double Scrwl4HBondInteraction::getW() {
	CartesianPoint e0(0,0,0); 
	CartesianPoint e1(0,0,0); 
	CartesianPoint e2(0,0,0); 
	CartesianPoint n(0,0,0); 
	
	// atoms[0] actual donor 
	// atoms[1] bonded to donor
	// atoms[2] actual acceptor
	// atoms[3] bonded to acceptor
	// atoms[4] bonded to acceptor_2 

	e0 = (pAtoms[0]->getCoor() - pAtoms[1]->getCoor());
	e1 = (CartesianGeometry::buildRadians(pAtoms[2]->getCoor(),pAtoms[3]->getCoor(),pAtoms[4]->getCoor(),params[0],params[1],params[2])) - pAtoms[2]->getCoor();
	n = (pAtoms[0]->getCoor()-pAtoms[2]->getCoor());


	double d = pAtoms[0]->distance(*pAtoms[2]);

	double t1 = (sig_d *sig_d) - (d - d0) * (d -d0);

	double cos_alpha = cos(CartesianGeometry::angleRadians((n * -1.0),e0));
	double t2 = cos_alpha - cos_alpha_max;
	double cos_beta = cos(CartesianGeometry::angleRadians(n,e1));

	double w = 0;
	if(t1 > 0 && t2 > 0) {
		double t3 = cos_beta - cos_beta_max;
		//cout << "UUU t3 " << t3 << endl;
		//cout << "UUU cos_beta " << cos_beta << endl;
		if(t3 > 0) {
			w = sqrt(t1 * t2 * t3)/denominator; 
		} else {
			e2 = (CartesianGeometry::buildRadians(pAtoms[2]->getCoor(),pAtoms[3]->getCoor(),pAtoms[4]->getCoor(),params[0],params[1],params[3])) - pAtoms[2]->getCoor();
			cos_beta = cos(CartesianGeometry::angleRadians(n,e2));
			t3 = cos_beta - cos_beta_max;
			//cout <<	"e2:" << e2 << endl; 
			//cout << "UUU cos_beta " << cos_beta << endl;
			//cout << "UUU t3 " << t3 << endl;
			if(t3 > 0) {
				w = sqrt(t1 * t2 * t3)/denominator; 
			}
		}
	}        
	
	if(w > 0 && debugFlagOn) {
		cout << " UUU w: " << w << endl; 
		cout << " UUU Donor " << *pAtoms[0] << endl;
		cout << " UUU Bonded to Donor " << *pAtoms[1] << endl;
		cout << " UUU Acceptor " << *pAtoms[2] << endl;
		cout << " UUU Bonded to Acceptor " <<*pAtoms[3] << endl;
		cout << " UUU Angled to Acceptor " <<  *pAtoms[4] << endl;
		cout << " UUU e0: " << e0 + pAtoms[0]->getCoor() << endl; 
		cout << " UUU e1: " << e1 + pAtoms[2]->getCoor() << endl; 
		cout << " UUU e2: " << e2 + pAtoms[2]->getCoor() << endl; 
		cout << " UUU n: " << n + pAtoms[0]->getCoor() << endl; 
		cout << " UUU Energy: " << (scalingFactor * w * pAtoms[0]->getCharge() * pAtoms[2]->getCharge() * B) << endl; 
		cout << " UUU d: " << d << endl;
		cout << " UUU t1 " << t1 << endl;
		cout << " UUU cos_alpha " << cos_alpha << endl;
		cout << " UUU t2 " << t2 << endl;
		
	}
	return w;
}	
 
double Scrwl4HBondInteraction::getEnergy() {
	return(scalingFactor * getW() * pAtoms[0]->getCharge() * pAtoms[2]->getCharge() * B);
}
