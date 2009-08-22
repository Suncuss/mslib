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

#include "PhiPsiStatistics.h"


PhiPsiStatistics::PhiPsiStatistics(){
}
PhiPsiStatistics::PhiPsiStatistics(const PhiPsiStatistics &_phiPsiStat){
	copy(_phiPsiStat);
}
PhiPsiStatistics::~PhiPsiStatistics(){
}

void PhiPsiStatistics::operator=(const PhiPsiStatistics &_phiPsiStat){
	copy(_phiPsiStat);
}
void PhiPsiStatistics::copy(const PhiPsiStatistics &_phiPsiStat){
	map<string,int>::iterator it;

    phiPsiTable = _phiPsiStat.getPhiPsiCounts();
}

void PhiPsiStatistics::addStatisitics(string _residueType, string _phiBin, string _psiBin, int _count){
	stringstream ss;
	ss << _residueType<<":"<<_phiBin<<":"<<_psiBin;
	
	//cout << "Adding: "<<ss.str()<<" "<<_count<<endl;
	phiPsiTable[ss.str()]      = _count;
	phiPsiTable[_residueType] += _count;
}


int PhiPsiStatistics::operator()(string _key){

	map<string,int>::iterator it;
	it = phiPsiTable.find(_key);
	if (it == phiPsiTable.end()){
		return MslTools::intMax;
	}

	return it->second;
}


int PhiPsiStatistics::getCounts(Residue nMinus1, Residue n, Residue nPlus1){	

	double phi = getPhi(nMinus1, n);
	double psi = getPsi(n, nPlus1);


	stringstream ss;
	ss << n.getResidueName() <<":"<<getPhiPsiBin(phi)<<":"<<getPhiPsiBin(psi);

	//cout << "Key: "<<ss.str()<<"."<<endl;
	map<string,int>::iterator it;
	it = phiPsiTable.find(ss.str());
	if (it == phiPsiTable.end()){
		cout << "No Phi/Psi entry found for " << ss.str() << "." <<endl;
		return MslTools::intMax;
	}
	return it->second;

}

double PhiPsiStatistics::getProbability(Residue nMinus1, Residue n, Residue nPlus1){

	/*
	  probability(x,y,z) = 
		  (#AA-Phi-Psi(x,y,z)   / #AA(x)) 
	  
	 */

	
	int AAxyz = getCounts(nMinus1, n, nPlus1);
	if (AAxyz == MslTools::intMax) {
		return MslTools::doubleMax;
	}


	map<string,int>::iterator it;
	it = phiPsiTable.find(n.getResidueName());
	if (it == phiPsiTable.end()) {
		return MslTools::doubleMax;
	}
	int AAx   = it->second;


	double AAxyzDouble = (double)AAxyz;
	double AAxDouble   = (double)AAx;

	//double mlogprob = 0.0;
	if (AAxyzDouble > 0.00001 && AAxDouble > 0.00001){
		//mlogprob = -log(AAxyzDouble / AAxDouble);
		return AAxyzDouble/ AAxDouble;
	}

	//return mlogprob;
	return 0.0;
}


double PhiPsiStatistics::getProbabilityAll(Residue nMinus1, Residue n, Residue nPlus1){
	/*
		  (#AA-Phi-Psi(all,y,z) / #AA(all)) 
	*/
	int AAallyz = 0;
	int AAall   = 0;

	double phiBin = getPhiPsiBin(getPhi(nMinus1,n));
	double psiBin = getPhiPsiBin(getPsi(n,nPlus1));
	stringstream allkey;
	allkey << "ALL"<<":"<<phiBin<<":"<<psiBin;

	map<string,int>::iterator it;
	it = phiPsiTable.find(allkey.str());
	if (it == phiPsiTable.end()){
		return MslTools::doubleMax;
	}
	AAallyz = it->second;

	
	double AAallyzDouble = (double)AAallyz;


	it = phiPsiTable.find("ALL");
	if (it == phiPsiTable.end()){
		return MslTools::doubleMax;
	}
	AAall = it->second;

	double AAallDouble   = (double)AAall;

	return AAallyzDouble / AAallDouble;
	
}
double PhiPsiStatistics::getPropensity(Residue nMinus1, Residue n, Residue nPlus1){

	/*
	  propensity(x,y,z) = 
		  (#AA-Phi-Psi(x,y,z)   / #AA(x)) 
		  -------------------------------
		  (#AA-Phi-Psi(all,y,z) / #AA(all)) 
	  
	 */

	double probRes = getProbability(nMinus1,n,nPlus1);
	double probAll = getProbabilityAll(nMinus1,n,nPlus1);

	if (probRes < 0.0001 || probAll < 0.0001) 
		return 0.0;

        if( (probRes == MslTools::doubleMax) || (probAll == MslTools::doubleMax))
            return MslTools::doubleMax;

	return probRes/probAll;

}


void PhiPsiStatistics::computeTotalCounts(){


	map<string,int>::iterator phiPsiIt;
	map<string,int> runningTotalByRes;
	int runningTotal = 0;
	int runningTotalbySubtotal = 0;
	for (phiPsiIt = phiPsiTable.begin(); phiPsiIt != phiPsiTable.end();phiPsiIt++){
		
		// Get key, get tokens. if token size = 3, then do counts
		string key = phiPsiIt->first;
		vector<string> toks = MslTools::tokenize(key,":");
		if (toks.size() == 3){
			if (toks[0] != "ALL"){
				stringstream newkey;
				newkey << "ALL"<<":"<<toks[1]<<":"<<toks[2];
				phiPsiTable[newkey.str()] += phiPsiIt->second;
				runningTotal += phiPsiIt->second;
				runningTotalByRes[toks[0]] += phiPsiIt->second;
			}

		} else {
			if (toks[0] != "ALL") {
				//cout << "Subtotal-by-read-in for "<<phiPsiIt->first<<" is "<< phiPsiIt->second<<endl;
				runningTotalbySubtotal += phiPsiIt->second;
			}
		}
	}

	//for (phiPsiIt = runningTotalByRes.begin(); phiPsiIt != runningTotalByRes.end();phiPsiIt++){
		//cout << "Subtotal-by-res for "<<phiPsiIt->first<<" is "<< phiPsiIt->second<<endl;
	//}
	if (runningTotalbySubtotal != runningTotal){
		cout << "ERROR 3333 subtotals don't match: "<<runningTotalbySubtotal<<" and "<<runningTotal<<endl;
	}

	phiPsiTable["ALL"] = runningTotal;
}

double PhiPsiStatistics::getPhiPsiBin(double _in){

	double gridSize = 5.0;
	double out;
	int cint = int(_in*100);
	int gint = int(gridSize *100);
	double remainder = gint - (cint % gint);

	out = ( int((_in*100 - (gridSize*100 - remainder))/100));
	if (out < 0) {
		out -= (gridSize / 2);
	} else {
		out += (gridSize / 2);
	}
	
	return out;
}
double PhiPsiStatistics::getPhi(Residue nMinus1, Residue n){
	if (!(nMinus1.exists("C") && n.exists("N") && n.exists("CA") && n.exists("C"))){
		return MslTools::doubleMax;
	}
	return CartesianGeometry::instance()->dihedral(nMinus1("C").getCoor(), n("N").getCoor(), n("CA").getCoor(), n("C").getCoor());
}

double PhiPsiStatistics::getPsi(Residue n, Residue nPlus1){
	if (!(n.exists("N") && n.exists("CA") && n.exists("C") && nPlus1.exists("N"))){
		return MslTools::doubleMax;
	}
	return CartesianGeometry::instance()->dihedral(n("N").getCoor(),n("CA").getCoor(), n("C").getCoor(), nPlus1("N").getCoor());
}

