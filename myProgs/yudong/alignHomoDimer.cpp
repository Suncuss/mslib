#include <iostream>
#include <fstream>
#include <cstdlib>
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
string programDescription = "This program creates helix dimers of specified geometries and find the one that has the smallest RMSD value after align with the given target helix.";
string programAuthor = "Yudong Sun and Samantha Anderson";
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

        
        // Declare system

        CartesianPoint origin(0.0,0.0,0.0);
        CartesianPoint xAxis(1.0,0.0,0.0);
        CartesianPoint yAxis(0.0,1.0,1.0);
        CartesianPoint zAxis(0.0,0.0,1.0);
        Transforms trans;

        // Read in the target helix 
        System ip_helix;
        if(!ip_helix.readPdb(opt.pdbFile)){
                cerr << "Fail to read the file: " << opt.pdbFile << endl;
        }
        AtomSelection ip_helix_sl(ip_helix.getAtomPointers());
        AtomPointerVector ip_helix_ca = ip_helix_sl.select("target_bb, name CA");

        // Read in the polyG helix
        System pf_helix;
        if(!pf_helix.readPdb(opt.helixFile)){
                cerr << "Fail to read file: " << "XXX" << endl;
        }
        AtomPointerVector chainA = pf_helix.getAtomPointers();
        pf_helix.saveCoor("initialState");


        // Axial Rotation Loop
        for(double axialRot = opt.axialRotStart; axialRot < opt.axialRotEnd; axialRot += opt.axialRotSteps){
                // Z Shift Loop
                for(double zShift = opt.zShiftStart; zShift < opt.zShiftEnd; zShift += opt.zShiftSteps){
                        pf_helix.applySavedCoor("initialState");
                        trans.rotate(chainA,axialRot,origin,zAxis);
                        trans.Ztranslate(chainA,zShift);
                        pf_helix.saveCoor("AxialZState");
                        // Crossing Angle Loop
                        for(double crossingAngle = opt.crossingAngleStart/2; crossingAngle < opt.crossingAngleEnd/2; crossingAngle += opt.crossingAngleSteps){
                                pf_helix.applySavedCoor("AxialZState");
                                trans.rotate(chainA, crossingAngle,origin,xAxis);
                                pf_helix.saveCoor("crossingAngleState");
                                // X Shift Loop
                                for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps){
                                        pf_helix.applySavedCoor("crossingAngleState");
                                        trans.Xtranslate(chainA,0.5*xShift);
                                        // Creat Chain B
                                        System yetAnother_pf_helix(chainA);
                                        AtomPointerVector chainB = yetAnother_pf_helix.getAtomPointers();
                                        trans.Zrotate180(chainB);
                                        for(int i = 0; i < chainB.size(); i++){
                                                chainB[i]->setChainId("B");
                                        }
                                        System model(chainA+chainB);
                                        AtomSelection model_sl(model.getAtomPointers());
                                        AtomPointerVector model_ca = model_sl.select("reference_bb, name CA");
                                        // RMSD Alignment
                                        trans.rmsdAlignment(model_ca,ip_helix_ca,model.getAtomPointers());
                                        // Store Possible Solutions
                                        solutions.push_back(Solution {ip_helix_ca.rmsd(model_ca),axialRot,zShift,crossingAngle,xShift});

                                }
                        }
                }
        }

        // Sort Solutions
        sort(solutions.begin(),solutions.end());
        

        // Redirect Output
        string outputFileName = opt.outputDir + "/" + opt.output;
        if(!freopen(outputFileName.c_str(), "w", stdout)){
                cerr << "ERROR_freopen: Did not open output file " << opt.output << endl;
		exit(0);
	}
 
        // Print RMSD info
        for(int i = 0; i < solutions.size(); i++){
                cout << "RMSD: " << solutions[i].rmsd << "\tPARAMETERS: " << printInfo(solutions[i]) << endl;
        }

        // Generate the best solution
        Solution bestSolution = solutions.front();
        pf_helix.applySavedCoor("initialState");
        trans.rotate(chainA,bestSolution.axialRot,origin,zAxis);
        trans.Ztranslate(chainA,bestSolution.zShift);
        trans.rotate(chainA, bestSolution.crossingAngle,origin,xAxis);
        trans.Xtranslate(chainA,0.5*bestSolution.xShift);
        System yetAnother_pf_helix(chainA);
        AtomPointerVector chainB = yetAnother_pf_helix.getAtomPointers();
        trans.Zrotate180(chainB);
        for(int i = 0; i < chainB.size(); i++){
                chainB[i]->setChainId("B");
        }
        System bestModel(chainA+chainB);
        string bestModelFileName = opt.outputDir + "/best_model.pdb";
        bestModel.writePdb(bestModelFileName);
        AtomSelection bestModel_sl(bestModel.getAtomPointers());
        AtomPointerVector bestModel_ca = bestModel_sl.select("best_model_bb, name CA");

        trans.rmsdAlignment(bestModel_ca,ip_helix_ca,bestModel.getAtomPointers());

        string bestAlignmentFileName = opt.outputDir + "/best_alignment.pdb";
        bestModel.writePdb(bestAlignmentFileName);

                                
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
		opt.startResNum = 1;
		opt.warningMessages += "starting residue number not specified; defaulting to 1";
		opt.warningFlag = true;
	}
       	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "ending residue number not specified";
		opt.warningFlag = true;
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
	cout << " --xShiftStart <double> --xShiftEnd <double> --xShiftSteps <double> " << endl;
	cout << " --zShiftStart <double> --zShiftEnd <double> --zShiftSteps <double> " << endl;
	cout << " --axialRotStart <double> --axialRotEnd <double> --axialRotSteps <double> " << endl;
	cout << " --crossingAngleStart <double> --crossingAngleEnd <double> --crossingAngleSteps <double> " << endl;
	cout << endl;
}