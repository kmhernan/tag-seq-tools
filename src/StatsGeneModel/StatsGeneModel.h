/*******************************************************************************
  StatsGeneModel.h
  Sep 2013

  Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
*******************************************************************************/

#ifndef __STATS_GENE_MODEL_H__
#define __STATS_GENE_MODEL_H__

#include "api/BamReader.h"
#include "IntervalTree.h"
#include "GffFileReader.h"

#include <string>
#include <iostream>
#include <fstream>
#include <map>

struct Counts {
    // Read-level stats
    int total;
    int reverse;    
    int forward;
    int unique;
    int duplicate;
    int overlap;
    int rev_overlap;
    int frd_overlap;
    int nonoverlap;
    int rev_nonover;
    int frd_nonover;
    int rev_one_gm;
    int frd_one_gm;
    int rev_multi_gm;
    int frd_multi_gm;

    // Gene-level stats
    int total_pac; 
    int rev_pac_hit;
    int frd_pac_hit;
    int rev_dup_pac;
    int rev_unique_pac;
    int frd_dup_pac;
    int frd_unique_pac;

    Counts() {total = 0; reverse = 0; forward = 0; unique = 0; duplicate = 0;
	      overlap = 0; rev_overlap = 0; frd_overlap = 0; nonoverlap = 0;
	      rev_nonover = 0; frd_nonover = 0; rev_one_gm = 0; frd_one_gm = 0;
              rev_multi_gm = 0; frd_multi_gm = 0; total_pac = 0; rev_dup_pac = 0;
              rev_unique_pac = 0; frd_dup_pac = 0; frd_unique_pac = 0;} 
};

class TagSeqGMStats {

public:
    // constructor
    TagSeqGMStats(std::string &infile, std::string &out, std::string &ingff);

    // destructor
    ~TagSeqGMStats(void);

private:

    // IO variables
    std::string _infile;
    std::string _out;
    std::string _ingff;

    // Initialize GffFileReader pointer
    GffFileReader * _gff;

    // Initialize Counts pointer
    Counts * _counts;

    // Containers
    std::map<std::string, std::map<char, IntervalTree<Gff> > > _tree_list;
    
    // Functions
    void Run();
    void GetStatistics();
    void PrintSummary();
};

#endif // __STATS_GENE_MODEL_H__
