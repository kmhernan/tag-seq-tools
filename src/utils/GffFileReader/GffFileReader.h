/*******************************************************************************
  GffFileReader.h 
  Aug 2013

  Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
*******************************************************************************/
#ifndef __GFF_FILE_READER_H__
#define __GFF_FILE_READER_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <stdint.h>
#include <cstdio>
#include <vector>
#include <map>
#include "IntervalTree.h"

// Data typedefs
typedef uint32_t POS;

// Structure for Gff records
struct Gff {
    std::string scaffold;
    std::string type;
    std::string gene;
    char ori;
    POS start;
    POS end;
};

// Type defs for data structures 
typedef std::vector<Interval<Gff> > GffInterval;
typedef std::map<std::string, std::map<char, IntervalTree<Gff> > > IntervalTreeMap;

// Class for reading Gff Files
class GffFileReader {

public:
    // Constructor/Destructor
    GffFileReader(std::string&);
    ~GffFileReader(void);
 
    // the Gff file 
    std::string gffFile;
    
    // File management
    // Opens GFF file for reading
    void Open(void);
    // Closes an open GFF file
    void Close(void);
    // Converts the current line string into the Gff struct 
    Gff& GetRecord(std::string&, Gff&);
    // Loads the GFF file into a map of IntervalTrees<Gff> keyed by scaffold and strand.
    // Intervals depend on the user input query
    void LoadIntervals(int&, IntervalTreeMap&);

private:
    // Input file and streams
    std::istream *_gffStream;
    std::string _file;
};

#endif // GFF_FILE_READER_H
