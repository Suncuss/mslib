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
#include <string>
#include <map>
#include <fstream>
#include <signal.h>
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "Timer.h"
#include "energyOptimizations.h"


using namespace std;

using namespace MSL;

// Global objects
MonteCarloOptimization mc;
Timer t;
double startTime = 0.0;

int main(int argc, char *argv[]){

	// Option Parser
	cout << "Setup Options"<<endl;
	MonteCarloOptions opt = setupMonteCarloOptions(argc, argv);
	cout << "Create System"<<endl;
	// Create a system from the structural input options
	startTime = t.getWallTime();
	System sys;
	createSystem(opt.structOpt, sys);
	cout << "Built system after "<<(t.getWallTime()-startTime)<<" seconds"<<endl;
	startTime = t.getWallTime();

	SelfPairManager scom;
	scom.setSystem(&sys);
	scom.seed(opt.randomSeed);
	scom.setOnTheFly(true);
	EnergySet *e = sys.getEnergySet();
	cout << "SIZEOF: "<<sizeof(*e)<<endl;

	// Calculate Fixed, Self but not Pair Energies (due to on-the-fly flag set)
	if (opt.DEE){
	  scom.setVerbose(true);
	  scom.setEnumerationLimit(10000);
	  scom.setRunDEE(true,true);
	  scom.runOptimizer();
	  scom.setVerbose(false);
	}

	scom.calculateEnergies();
	cout << "Calc Fixed/Self energies after "<<(t.getWallTime()-startTime)<<" seconds"<<endl;
	startTime = t.getWallTime();

	// Setup MCOpt object
	mc.setSelfPairManager(&scom);
	mc.setNumberOfStoredConfigurations(opt.numStoredConfigurations);
	mc.setInitializationState(opt.initAlgo);

	// Random Seed
	mc.seed(opt.randomSeed);
	
	// Run MC
	cout << "Run MC with "<<mc.getNumPositions()<<" positions and "<<mc.getStateEnergy()<<" total energy. random seed: "<<opt.randomSeed<<endl;
	mc.runMC(opt.annealStart,opt.annealEnd,opt.numCycles,opt.annealShape,opt.maxRejections,opt.deltaSteps,opt.minDeltaE);

	fprintf(stdout, "Random Seed Used: %d\n",mc.getSeed());
	cout << "MC Run after "<<(t.getWallTime()-startTime)<<" seconds"<<endl;

	// Either print rotamer selections + energies out, or generate PDBs and print out
	if (opt.structureConfig == ""){

		// Print results..
		cout << "Sampled Energies: "<<endl;
		mc.printSampledConfigurations();
	} else {

		// Get priority queue of resulting conformations
		priority_queue< pair<double,string>, vector< pair<double,string> >, less<pair<double,string> > > &conformations = mc.getSampledConformations();

		int solution = 1;
		while (!conformations.empty()){
		        cout << "Working on solution conformation "<<solution<<" energy "<<conformations.top().first <<" "<<conformations.top().second<<endl;
			
			
			vector<string> toks = MslTools::tokenize(conformations.top().second,":");
			vector<int> rotamerState;
			for (uint i = 0; i < toks.size();i++){
				rotamerState.push_back(MslTools::toInt(toks[i]));
			}
			conformations.pop();

			// Helper function takes structOptions, a System and a rotamer state , putting system into given rotamer state.
			changeRotamerState(opt.structOpt,sys,rotamerState);

			
			string sysString = "";
			for (uint c = 0; c < sys.chainSize();c++){
			  sysString += MslTools::stringf("%1s: %s\n",sys.getChain(c).getChainId().c_str(),PolymerSequence::toOneLetterCode(sys.getChain(c).getAtomPointers(),"CA").c_str());
			}

			fprintf(stdout,"Energy %04d: %8.3f\n%s\n",solution,sys.getEnergySet()->calcEnergy(),sysString.c_str());

			// Write out PDB
			char name[80];
			sprintf(name, "winnerMC-%04d.pdb",solution++);
			sys.writePdb(name);
		
		}
				
		
		
				
	}
	
}






void cleanExit(int sig) {

	cout << "*************SIGNAL CAUGHT*****************\n";
	cout << "\tSIG"<<sig<<endl;

	fprintf(stdout, "Random Seed Used: %d\n",mc.getSeed());
	cout << "Best-Sampled Energies after "<<(t.getWallTime()-startTime)<<" seconds."<<endl;
	mc.printSampledConfigurations();


	cout << "GoodBye."<<endl;
	exit(0);
}







