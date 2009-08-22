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

#ifndef FILE_H
#define FILE_H

// STL Includes
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


using namespace std;

class File {
	public:
		enum fileHandlerType { cstyle=0, cppstyle=1, stringstyle=2 };
		enum fileModes { read=0, write=1, append=2 };

		// Constructors/Destructors
		File(const string &_filename, int _mode);
		File(stringstream &_ss);
		File(const File &_anotherFile);
		void operator=(const File &_anotherFile);
	
		virtual ~File();
		void init(const string &_filename, int _mode);  // all constructors call this

		// Set/Get functions
		void setFileName(const string &_filename);
		string getFileName() const;
	
		void setOpenMode(int _mode);
		fileModes getOpenMode() const;

		void addRemark(const string &_remark);
		vector<string> getRemarks() const;
		void clearRemarks();

		void setFileHandlerType(int _fileHandler);
		fileHandlerType getFileHandlerType() const;

	
		// Member functions
		bool open();               // There is a default implementation
		bool open(const string &_filename); // There is a default implementation
		bool open(const string &_filename, int mode); // There is a default implementation
		bool open(stringstream &_ss);
		bool is_open();

		void close();

		bool doesFileExist();

		bool endOfFileTest();



	protected:
		void copy(const File &_anotherFile);

		// Different types of file handlers
		fstream fileStream;
		FILE *filePtr;
		stringstream *stringStreamPtr;

		// Meta data variables
		string  fileName;
		vector<string>  remarks; // Put in top of files
		fileHandlerType fileHandler;
		fileModes fileOpenMode;

	private:
		File(); // Default constructor is private , don't know if reader or writer.
       
};	

// INLINES GO HERE
inline File::File(const File &_anotherFile) { copy(_anotherFile); }
inline void File::operator=(const File &_anotherFile) { copy(_anotherFile); }
inline void               File::setFileName(const string &_filename) { fileName = _filename; }
inline string             File::getFileName() const { return fileName; }
inline void               File::addRemark(const string &_remark) { remarks.push_back(_remark); }
inline vector<string>     File::getRemarks() const     { return remarks; }

#endif
