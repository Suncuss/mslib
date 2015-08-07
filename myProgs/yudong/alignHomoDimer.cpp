#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "System.h"
#include "MslTools.h"
#include "Transforms.h"
#include "OptionParser.h"
#include "PDBReader.h"
#include "AtomSelection.h"
#include "Transforms.h"

using namespace std;
using namespace MSL;


string programName = "alignHomoDimer";
string programDescription = "This program creates  perfect helix dimers of specified geometries and then find the one that has the smallest RMSD value after align with the given imperfect helix.";
string programAuthor = "Yudong Sun, Samantha Anderson and Samson Condon";
string programVersion = "0.0.1";
string programDate = "31 July 2015";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

struct Options {
	string pdbFile;
        string helixFile;
	int startResNum;
        int endResNum;
	string output;
	string outputDir;

	double xShiftStart;
	double zShiftStart;
	double axialRotStart;
	double crossingAngleStart;
	double xShiftEnd;
	double zShiftEnd;
	double axialRotEnd;
	double crossingAngleEnd;
	double xShiftSteps;
	double zShiftSteps;
	double axialRotSteps;
	double crossingAngleSteps;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
};

struct Solution
{
        double rmsd;
        double axialRot;
        double zShift;
        double crossingAngle;
        double xShift;
        bool operator < (const Solution& sol) const
        {
                return(rmsd < sol.rmsd);
        }
};

void printOptions(Options& _opt) {
	cout << "Program            " << programName << " v." << programVersion << "," << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl; 
	cout << "pdbFile            " <<  _opt.pdbFile << endl;
        cout << "helixFile          " <<  _opt.helixFile << endl;
	cout << "output             " <<  _opt.output << endl;
        cout << "outputDir          " <<  _opt.outputDir << endl;
        cout << "startResNum        " <<  _opt.startResNum << endl;
        cout << "endResNum          " <<  _opt.endResNum << endl;
	cout << "xShiftStart        " <<  _opt.xShiftStart << endl;
	cout << "zShiftStart        " <<  _opt.zShiftStart << endl;
	cout << "axialRotStart      " <<  _opt.axialRotStart << endl;
	cout << "crossingAngleStart " <<  _opt.crossingAngleStart << endl;
	cout << "xShiftEnd          " <<  _opt.xShiftEnd << endl;
	cout << "zShiftEnd          " <<  _opt.zShiftEnd << endl;
	cout << "axialRotEnd        " <<  _opt.axialRotEnd << endl;
	cout << "crossingAngleEnd   " <<  _opt.crossingAngleEnd << endl;
	cout << "xShiftSteps        " <<  _opt.xShiftSteps << endl;
	cout << "zShiftSteps        " <<  _opt.zShiftSteps << endl;
	cout << "axialRotSteps      " <<  _opt.axialRotSteps << endl;
	cout << "crossingAngleSteps " <<  _opt.crossingAngleSteps << endl;
}

void usage();
void version();
void help(Options defaults);
string printInfo(Solution sol);
Options parseOptions(int _argc, char * _argv[], Options defaults);
     
int main(int argc, char **argv){
        time(&startTime);
        vector<Solution> solutions;
        Options defaults; 
        Options opt = parseOptions(argc, argv, defaults);
        
        if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;
		usage();
		exit(1);
	}	
	printOptions(opt);

        
        // Setup transfroms
        CartesianPoint origin(0.0,0.0,0.0);
        CartesianPoint xAxis(1.0,0.0,0.0);
        CartesianPoint yAxis(0.0,1.0,1.0);
        CartesianPoint zAxis(0.0,0.0,1.0);
        Transforms trans;

        // Read in the imperfect helix 
        System imperfectHelix;
        if(!imperfectHelix.readPdb(opt.pdbFile)){
                cerr << "Fail to read the file: " << opt.pdbFile << endl;
                exit(1);
        }
        AtomSelection imperfectHelixSelection(imperfectHelix.getAtomPointers());
        char targetSelString[100]; sprintf(targetSelString, "target_bb, name CA and resi %d-%d", opt.startResNum, opt.endResNum);
        AtomPointerVector imperfectHelixCA = imperfectHelixSelection.select(targetSelString);

        // Read in the 69-polyG perfect helix
        System perfectHelix;
        if(!perfectHelix.readPdb(opt.helixFile)){
                cerr << "Fail to read file: " << opt.helixFile << endl;
                exit(1);
        }
        AtomSelection perfectHelixSelection(perfectHelix.getAtomPointers());
        
        // Cut the perfect helix into the same length as the imperfect helix
        unsigned int helixLength = opt.endResNum - opt.startResNum + 1;
        if(helixLength > 69){
                cerr << "Helix length exceeds 69 - Tooooooo long >_< " << endl;
                exit(1);
        }
        unsigned int startRes, endRes;
        if( helixLength % 2 != 0){
                startRes = 35 - helixLength / 2; endRes = 35 + helixLength / 2; 
        }

        // When the desired length is not a odd number, the first half part will have one more resi than the second half of the helix
        else {
                startRes = 35 - helixLength / 2; endRes = 35 + helixLength / 2 - 1;
        }

        char refSelString[100]; 
        sprintf(refSelString, "ref, resi %d-%d", startRes,endRes);
        System perfectHelixDimer(perfectHelixSelection.select(refSelString));
        AtomPointerVector perfectHelixDimerAtoms = perfectHelixDimer.getAtomPointers();
        AtomSelection perfectHelixDimerSelection(perfectHelixDimerAtoms);
        AtomPointerVector chainA = perfectHelixDimerSelection.select("ref_chain_A, chain A");
        AtomPointerVector chainB = perfectHelixDimerSelection.select("ref_chain_B, chain B");
        AtomPointerVector perfectHelixDimerCA = perfectHelixDimerSelection.select("ref_bb, name CA");
        
        
        // Save the initial state of the system for later use
        perfectHelixDimer.saveCoor("initialState");

        
        // Save coor for atoms manually to avoid string compare for better performance
        // Declare containers for atom's coor
        unsigned int chain_size = chainA.size();
        vector<double> atoms_x_init(chain_size); vector<double> atoms_y_init(chain_size); vector<double> atoms_z_init(chain_size);
        vector<double> atoms_x_axialZ(chain_size); vector<double> atoms_y_axialZ(chain_size); vector<double> atoms_z_axialZ(chain_size);
        vector<double> atoms_x_crossing(chain_size); vector<double> atoms_y_crossing(chain_size); vector<double> atoms_z_crossing(chain_size);        
        for(int i = 0; i < chain_size; i++){
                Atom * atom = chainA[i];
                atoms_x_init[i] = atom->getX();
                atoms_y_init[i] = atom->getY();
                atoms_z_init[i] = atom->getZ();
        }
        
        // Check the range of parameter before do the real stuff 
        if(opt.axialRotStart > opt.axialRotEnd || opt.zShiftStart > opt.zShiftEnd || opt.crossingAngleStart > opt.crossingAngleEnd || opt.xShiftStart > opt.xShiftEnd){
                cerr << "Erroooooor: Range of parameters makes zero sense ~!" << endl;
                exit(1);
        }

        // Axial Rotation Loop
        for(double axialRot = opt.axialRotStart; axialRot < opt.axialRotEnd; axialRot += opt.axialRotSteps){
                // Z Shift Loop
                for(double zShift = opt.zShiftStart; zShift < opt.zShiftEnd; zShift += opt.zShiftSteps){
                        // Apply saved coor "initialState"
                        for(int i = 0; i < chain_size; i++){
                                chainA[i]->setCoor(atoms_x_init[i],atoms_y_init[i],atoms_z_init[i]);
                        }
                        trans.rotate(chainA,axialRot,origin,zAxis);
                        trans.Ztranslate(chainA,zShift);
                        // Save coor after transformation
                        for(int i = 0; i < chain_size; i++){
                                Atom * atom = chainA[i];
                                atoms_x_axialZ[i] = atom->getX();
                                atoms_y_axialZ[i] = atom->getY();
                                atoms_z_axialZ[i] = atom->getZ();
                        }
                        // Crossing Angle Loop
                        for(double crossingAngle = opt.crossingAngleStart; crossingAngle < opt.crossingAngleEnd; crossingAngle += opt.crossingAngleSteps){
                                for(int i = 0; i < chain_size; i++){
                                        chainA[i]->setCoor(atoms_x_axialZ[i],atoms_y_axialZ[i],atoms_z_axialZ[i]);
                                }
                                trans.rotate(chainA,crossingAngle/2.0,origin,xAxis);
                                for(int i = 0; i < chain_size; i++){
                                        Atom * atom = chainA[i];
                                        atoms_x_crossing[i] = atom->getX();
                                        atoms_y_crossing[i] = atom->getY();
                                        atoms_z_crossing[i] = atom->getZ();
                                }
                                // X Shift Loop
                                for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps){
                                        for(int i = 0; i < chain_size; i++){
                                               chainA[i]->setCoor(atoms_x_crossing[i],atoms_y_crossing[i],atoms_z_crossing[i]);   
                                        }
                                        trans.Xtranslate(chainA,0.5*xShift);
                                        // Apply the same Transfromations on chian B
                                        for(int i=0; i < chain_size; i++){
                                                Atom * atom = chainA[i];
                                                chainB[i]->setCoor(atom->getX(),atom->getY(),atom->getZ());
                                        }
                                        // Rotate chain B by 180 to set the symmetry
                                        trans.Zrotate180(chainB);
                                        // RMSD Alignment
                                        trans.slientRmsdAlignment(perfectHelixDimerCA,imperfectHelixCA,perfectHelixDimerAtoms);
                                        // Store Possible Solutions
                                        solutions.push_back(Solution {imperfectHelixCA.rmsd(perfectHelixDimerCA),axialRot,zShift,crossingAngle,xShift});
                                }
                        }
                }
        }

        // Sort all the solutions
        sort(solutions.begin(),solutions.end());
       
        // Redirect output
        string outputFileName = opt.outputDir + "/" + opt.output;
        if(!freopen(outputFileName.c_str(), "w", stdout)){
                cerr << "ERROR_freopen: Did not open output file " << opt.output << endl;
		exit(1);
	}
 
        // Print RMSD info
        for(int i = 0; i < solutions.size(); i++){
                cout << "RMSD: " << solutions[i].rmsd << "\tPARAMETERS: " << printInfo(solutions[i]) << endl;
        }

        // Generate pdb files for the best solution
        Solution bestSolution = solutions.front();
        perfectHelixDimer.applySavedCoor("initialState");
        trans.rotate(chainA,bestSolution.axialRot,origin,zAxis);
        trans.rotate(chainB,bestSolution.axialRot,origin,zAxis);
        trans.Ztranslate(chainA,bestSolution.zShift);
        trans.Ztranslate(chainB,bestSolution.zShift);
        trans.rotate(chainA,bestSolution.crossingAngle/2,origin,xAxis);
        trans.rotate(chainB,bestSolution.crossingAngle/2,origin,xAxis);
        trans.Xtranslate(chainA,0.5*bestSolution.xShift);
        trans.Xtranslate(chainB,0.5*bestSolution.xShift);
        trans.Zrotate180(chainB); 
        string bestModelFileName = opt.outputDir + "/best_model.pdb";
        perfectHelixDimer.writePdb(bestModelFileName);
        trans.slientRmsdAlignment(perfectHelixDimerCA,imperfectHelixCA,perfectHelixDimerAtoms);
        string bestAlignmentFileName = opt.outputDir + "/best_alignment.pdb";
        perfectHelixDimer.writePdb(bestAlignmentFileName);

                                
        time(&endTime);
        diffTime = difftime (endTime, startTime);
        cout << "Total Time: " << diffTime << " seconds" << endl;
        return 0;
}

string printInfo(Solution sol){
        char info[100];
        sprintf(info,"model_%03.3f_" "%03.1f_" "%03.3f_" "%03.3f" , sol.axialRot, sol.zShift, sol.crossingAngle, sol.xShift);
        return string(info);
}

Options parseOptions(int _argc, char * _argv[], Options defaults) {

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file 
	 ******************************************/
	
	Options opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuartion file:
	 *  
	 ******************************************/
	vector<string> required;
	vector<string> allowed;


	opt.required.push_back("pdbFile");
        opt.required.push_back("helixFile");
	opt.required.push_back("output");
	opt.required.push_back("outputDir");

	opt.required.push_back("xShiftStart");
	opt.required.push_back("xShiftEnd");
	opt.required.push_back("xShiftSteps");

	opt.required.push_back("zShiftStart");
	opt.required.push_back("zShiftEnd");
	opt.required.push_back("zShiftSteps");

	opt.required.push_back("axialRotStart");
	opt.required.push_back("axialRotEnd");
	opt.required.push_back("axialRotSteps");

	opt.required.push_back("crossingAngleStart");
	opt.required.push_back("crossingAngleEnd");
	opt.required.push_back("crossingAngleSteps");
        opt.required.push_back("startResNum");
        opt.required.push_back("endResNum");

	opt.allowed.push_back("configfile");
        

	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
	}
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version or -v arguments print the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");
	//if (OP.fail()) {
	//	opt.version = OP.getBool("v");
	//}

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");
//	if (OP.fail()) {
//		opt.help = OP.getBool("h");
//	}

	if (opt.help) {
		help(defaults);
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;

	/*****************************************
	 *  OUTPUT DIR AND FILES
	 *****************************************/


	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.warningMessages = "pdb file not specified";
		opt.errorFlag = true;
	}
        opt.helixFile = OP.getString("helixFile");
        if (OP.fail()) {
                opt.warningMessages = "pdb file not specified";
                opt.errorFlag = true;
        }
	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified";
		opt.errorFlag = true;
	}
       	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "ending residue number not specified";
		opt.errorFlag = true;
	}
	opt.output = OP.getString("output");
	if (OP.fail()) {
		opt.errorMessages = "output not specified";
		opt.errorFlag = true;
	}
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.warningMessages = "output directory not specified";
		opt.errorFlag = true;
	}
	opt.xShiftStart = OP.getDouble("xShiftStart");
	if (OP.fail()) {
		opt.errorMessages = "xShiftStart not specified";
		opt.errorFlag = true;
	}
	opt.xShiftEnd = OP.getDouble("xShiftEnd");
	if (OP.fail()) {
		opt.errorMessages = "xShiftEnd not specified";
		opt.errorFlag = true;
	}
	opt.xShiftSteps = OP.getDouble("xShiftSteps");
	if (OP.fail()) {
		opt.errorMessages = "xShiftSteps not specified";
		opt.errorFlag = true;
	}

	opt.zShiftStart = OP.getDouble("zShiftStart");
	if (OP.fail()) {
		opt.errorMessages = "zShiftStart not specified";
		opt.errorFlag = true;
	}
	opt.zShiftEnd = OP.getDouble("zShiftEnd");
	if (OP.fail()) {
		opt.errorMessages = "zShiftEnd not specified";
		opt.errorFlag = true;
	}
	opt.zShiftSteps = OP.getDouble("zShiftSteps");
	if (OP.fail()) {
		opt.errorMessages = "zShiftSteps not specified";
		opt.errorFlag = true;
	}

	opt.axialRotStart = OP.getDouble("axialRotStart");
	if (OP.fail()) {
		opt.errorMessages = "axialRotStart not specified";
		opt.errorFlag = true;
	}
	opt.axialRotEnd = OP.getDouble("axialRotEnd");
	if (OP.fail()) {
		opt.errorMessages = "axialRotEnd not specified";
		opt.errorFlag = true;
	}
	opt.axialRotSteps = OP.getDouble("axialRotSteps");
	if (OP.fail()) {
		opt.errorMessages = "axialRotSteps not specified";
		opt.errorFlag = true;
	}
	opt.crossingAngleStart = OP.getDouble("crossingAngleStart");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleStart not specified";
		opt.errorFlag = true;
	}
	opt.crossingAngleEnd = OP.getDouble("crossingAngleEnd");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleEnd not specified";
		opt.errorFlag = true;
	}
	opt.crossingAngleSteps = OP.getDouble("crossingAngleSteps");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleSteps not specified";
		opt.errorFlag = true;
	}

		return opt;

}

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help(Options defaults) {
	cout << "This program runs as:" << endl;
	cout << " % alignHomoDimer --pdbFile <pdbfile> --helixFile <filename> --output <filename>" << endl;
        cout << " --startResNum <int> --endResNum <int>" << endl;
	cout << " --xShiftStart <double> --xShiftEnd <double> --xShiftSteps <double> " << endl;
	cout << " --zShiftStart <double> --zShiftEnd <double> --zShiftSteps <double> " << endl;
	cout << " --axialRotStart <double> --axialRotEnd <double> --axialRotSteps <double> " << endl;
	cout << " --crossingAngleStart <double> --crossingAngleEnd <double> --crossingAngleSteps <double> " << endl;
	cout << endl;
}
