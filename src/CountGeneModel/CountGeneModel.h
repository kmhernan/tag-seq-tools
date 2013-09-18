/*******************************************************************************
  CountGeneModel.h
  Aug 2013

  Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
*******************************************************************************/

#ifndef __COUNT_GENE_MODEL_H__
#define __COUNT_GENE_MODEL_H__

#include "api/BamReader.h"
#include "IntervalTree.h"
#include "GffFileReader.h"

#include <string>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <random>
#include <set>
#include <math.h>

// Structure for summary counts
struct GeneralCounts
{
  int total;
  int overlapping;
  int nonOverlapping;
  int dups;
  int multiHitModel;
  int notPrimaryAln;
  int totalPass; 
  GeneralCounts(){total=0; overlapping = 0; nonOverlapping=0; 
		  dups=0; multiHitModel=0; notPrimaryAln=0;
		  totalPass= 0;}
};

// Class methods and elements
class TagSeqGMCounts {

public:

    // constructor
    TagSeqGMCounts(std::string &infile, std::string &ingff, 
		   std::string &out, std::string &nonoverlapping,
		   int& query, bool noDups, bool primaryAln,
                   bool randomOne, bool splitCounts, bool annotateDups);

    // destructor
    ~TagSeqGMCounts(void);

private:

    // IO files.
    std::string _infile;
    std::string _ingff;
    std::string _out;
    std::string _nonoverlapping;
    
    // query
    int _query;

    // Instance of GFF file class
    GffFileReader *_gff;

    // Duplicate read handling
    bool _noDups;
    bool _primaryAln;
    bool _randomOne;
    bool _splitCounts;
    bool _annotateDups;

    // Containers
    std::map<std::string, std::map<char, IntervalTree<Gff> > > _tree_list;
    std::map<std::string, int> _gene_counts;
    std::set<std::string> _dup_pac;
    GeneralCounts _counts;

    // Process the BAM
    void Run();
    void ProcessBam_noDups();
    void ProcessBam_primaryAln();
    void ProcessBam_randomOne();
    void ProcessBam_splitCounts();

    // Write results to file
    void WriteCounts();

    // Print reports
    void PrintSummary();
};
#endif // __COUNT_GENE_MODEL_H__
