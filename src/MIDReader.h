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

#ifndef MIDREADER_H
#define MIDREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"


// Storage Formats
#include "MoleculeInterfaceDatabase.h"


// STL Includes
#include <vector>
using namespace std;
using namespace MslTools;

/**
 * This class will be used to read a MoleculeInterfaceDatabase from
 * a binary file.
 */
class MIDReader : public Reader {

	public:
		// Constructors/Destructors
		MIDReader();
		MIDReader(const string &_filename);
		MIDReader(const MIDReader & _reader);
		MIDReader(stringstream &_stream);
		MIDReader(string &_string);
		virtual ~MIDReader();

		bool read(string _archiveType="text");

		// Get/Set
        /**
         * This method returns the MoleculeInterfaceDatabase read in.
         *
         * @return The MoleculeInterfaceDatabase.
         */
		MoleculeInterfaceDatabase & getMoleculeInterfaceDatabase() { return mid;}
		void reset();

		// Operators

	protected:		
	private:
		void deletePointers();

		MoleculeInterfaceDatabase mid;

};

//Inlines go HERE
/**
 * The default constructor.
 */
inline MIDReader::MIDReader() : Reader() {}
/**
 * This constructor takes the filename of the binary
 * file which holds the MoleculeInterfaceDatabase.
 *
 * @param _filename The name of the binary file holding the MoleculeInterfaceDatabase.
 */
inline MIDReader::MIDReader(const string &_filename) : Reader(_filename) {}
/**
 * A copy constructor.
 *
 * @param _reader The MIDReader to be copied.
 * @todo Does the MIDReader copy constructor actually copy anything currently?
 */
inline MIDReader::MIDReader(const MIDReader & _reader) { }
/**
 * A constructor that takes in a stringstream as input.
 *
 * @input _ss The stringstream to use to read in the MoleculeInterfaceDatabase.
 */
inline MIDReader::MIDReader(stringstream &_ss) : Reader(_ss)     {read();}
inline MIDReader::MIDReader(string &_string)   : Reader(_string) {read();}
/**
 * The deconstructor.
 */
inline MIDReader::~MIDReader() {deletePointers(); close();}
/**
 * This method will reset the reader.  Currently, reset doesn't really mean
 * much.  This object doesn't hold any pointers, so nothing to delete really.
 */
inline void MIDReader::reset() {deletePointers();}
/**
 * This method would delete all the pointers held by
 * this object if there were any pointers to delete.
 * Basically a place-holder in case this object
 * holds more information.
 */
inline void MIDReader::deletePointers(){}
#endif
