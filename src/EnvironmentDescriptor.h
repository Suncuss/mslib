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

#ifndef ENVIRONMENTDESCRIPTOR_H
#define ENVIRONMENTDESCRIPTOR_H


// MSL Includes
#include "AtomPointerVector.h"
#include "Frame.h"
#include "System.h"

// STL Includes
#include <map>

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#endif


using namespace std;

class EnvironmentDescriptor {

	public:

		EnvironmentDescriptor();
		EnvironmentDescriptor(EnvironmentDescriptor & _ed);
		~EnvironmentDescriptor();

		void operator=(EnvironmentDescriptor & _ed); // assignment

		// Complete setup 
		bool setupDescriptor(Residue  &_res, System &_sys, string type);
		
		// Get Set
		AtomPointerVector & getCore();
		void setCore(AtomPointerVector &_atoms);

		Frame & getReferenceFrame();
		void setReferenceFrame(Frame &_frame);
		
		Frame & getEnvironmentFrame(string _environmentType);
		AtomPointerVector & getEnvironment(string _environmentType);
		void setEnvironment(string _environmentType, AtomPointerVector &_atoms);

		map<string,AtomPointerVector*> & getEnvironmentMap();
		map<string,Frame*> & getFrameMap();

		string generateLookupKey(string _envType);

		void setName(string _name);
		string getName();



	private:

		void copy(EnvironmentDescriptor & _ed);


		map<string, AtomPointerVector*> environmentMap;
		map<string, Frame*>      frameMap;

		AtomPointerVector *core;
		Frame *frame;

		string name;

		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__

	public:

		void save_checkpoint(string filename) const{
			std::ofstream fout(filename.c_str());
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

		void load_checkpoint(string filename){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}
	private:

		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			ar & environmentMap;
			ar & frameMap;
			ar & core;
			ar & frame;
			ar & name;
		}

#endif
		
};
#endif
