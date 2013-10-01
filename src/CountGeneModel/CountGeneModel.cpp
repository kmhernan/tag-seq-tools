/*******************************************************************************
  CountGeneModel.cpp
  Aug 2013

  Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
*******************************************************************************/
#include "CountGeneModel.h"
using namespace BamTools;
using namespace std;

// build
TagSeqGMCounts::TagSeqGMCounts(string &infile, string &ingff,
			       string &out, string &nonoverlapping,
			       int &query, bool noDups, bool primaryAln,
 			       bool randomOne, bool splitCounts, bool annotateDups)
{

    _infile	    = infile;
    _ingff	    = ingff;
    _out	    = out;
    _nonoverlapping = nonoverlapping;

    _query	    = query;
    _gff	    = new GffFileReader(ingff);

    _noDups         = noDups;
    _primaryAln     = primaryAln;
    _randomOne      = randomOne;
    _splitCounts    = splitCounts;
    _annotateDups   = annotateDups;
 
    // Process
    Run();
}

// destroy
TagSeqGMCounts::~TagSeqGMCounts(void) {
    delete _gff;
}

/**
 * Wrapper function for running the application
 */
void TagSeqGMCounts::Run() {
    // load the GFF file
    _gff->LoadIntervals(_query, _tree_list);

    // Process counts
    if (_noDups)
        ProcessBam_noDups(); 
    else if (_primaryAln)
        ProcessBam_primaryAln();
    else if (_randomOne)
        ProcessBam_randomOne();
    else
        ProcessBam_splitCounts();

    // Write the results to a file.
    WriteCounts();

    // Print summaries to console.
    PrintSummary();
}

/**
 * Run GMCounts but remove all duplicate reads
 */
void TagSeqGMCounts::ProcessBam_noDups() {
    // Open the BAM file
    BamReader reader;
    if ( !reader.Open(_infile)) {
        cerr << "Failed to open BAM file " << _infile << endl;
        exit(1);
    }

    // See BamTools Api 
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    // Variables for reading the Bam file
    string scaffold;
    BamAlignment al;
    GffInterval results;
    string readName;
    string pacid;

    // Write the nonoverlapping sites to file.
    ofstream o(_nonoverlapping.c_str());

    // Cache for duplicate reads.
    map<string, string> saw_read;
    bool readDup;

    while ( reader.GetNextAlignment(al) ) {

        ++_counts.total;

        readName = al.Name;
 
        // If the read already has been processed, it is multi mapped.
        if (saw_read.count(readName) > 0)
            readDup = true;
        else {
            readDup = false;
            saw_read[readName] = "";
        }

        // Clear the results vector
        results.clear(); 

        // BWA Sometimes does weird things with the reference dictionary, so
        // this check is needed to skip a possible segmentation fault.
        if (al.RefID >= 0) {
            scaffold = references[al.RefID].RefName;
            // Check if the scaffold is even in the map 
            if (_tree_list.count(scaffold) > 0) {
                // Check the strandedness and if the map contains intervals for that 
                // strand on this scaffold
                if (al.IsReverseStrand() && _tree_list[scaffold].count('-') > 0)
                    _tree_list[scaffold]['-'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else if (!al.IsReverseStrand() && _tree_list[scaffold].count('+') > 0)
                    _tree_list[scaffold]['+'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                    continue;
                }
                // Deal with dups 
                if (results.size() == 1) {
                    ++_counts.overlapping;
                    pacid = results[0].value.gene;
                    if (!readDup) { 
                        saw_read[readName] = pacid;
                        _gene_counts[pacid]++;
                    }
                    else if (readDup && saw_read[readName] != "") {
                        _gene_counts[saw_read[readName]]--;
                        _counts.dups++;
                        saw_read[readName] = "";
                    }
                }

                else if (results.size() > 1) {
                    ++_counts.multiHitModel;
                }

                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                }
            }

            else {
                ++_counts.nonOverlapping;
                o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
            }
        }
    }
   
    o.close();
    reader.Close();
}

void TagSeqGMCounts::ProcessBam_primaryAln() {
    // Open the BAM file
    BamReader reader;
    if ( !reader.Open(_infile)) {
        cerr << "Failed to open BAM file " << _infile << endl;
        exit(1);
    }

    // See BamTools Api 
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    // Variables for reading the Bam file
    string scaffold;
    BamAlignment al;
    GffInterval results;
    string readName;
    string pacid;

    // Write the nonoverlapping sites to file.
    ofstream o(_nonoverlapping.c_str());

    while ( reader.GetNextAlignment(al) )
    {
        ++_counts.total;

        if (!al.IsPrimaryAlignment()) {
            ++_counts.notPrimaryAln;
            continue;
        }

        results.clear(); 

        // BWA Sometimes does weird things with the reference dictionary, so
        // this check is needed to skip a possible segmentation fault.
        if (al.RefID >= 0) {
            scaffold = references[al.RefID].RefName;
            // Check if the scaffold is even in the map; else it's non-overlapping 
            if (_tree_list.count(scaffold) > 0) {
                // Check the standedness and if the map contains intervals for that 
                // strand on this scaffold. It it's not found on either interval, it's
                // non-overlapping.
                if (al.IsReverseStrand() && _tree_list[scaffold].count('-') > 0)
                    _tree_list[scaffold]['-'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else if (!al.IsReverseStrand() && _tree_list[scaffold].count('+') > 0)
                    _tree_list[scaffold]['+'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                    continue;
                }
                // Now we count overlapping sites, if any, and deal with multiHitModel sites
                if (results.size() == 1) {
                    ++_counts.overlapping;
                    pacid = results[0].value.gene;
                    _gene_counts[pacid]++;
                }
                else if (results.size() > 1)
                    ++_counts.multiHitModel;
                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                }
            } 

            else {
                ++_counts.nonOverlapping;
                o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
            }
        }
    }
    
    o.close();
    reader.Close();
}

/**
 * Run GMCounts but choose a random one for multi-mapped reads.
 * Here we actually have to store PACIDs in a vector 
 */
void TagSeqGMCounts::ProcessBam_randomOne() {
    // Open the BAM file
    BamReader reader;
    if ( !reader.Open(_infile)) {
        cerr << "Failed to open BAM file " << _infile << endl;
        exit(1);
    }

    // See BamTools Api 
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    // Variables for reading the Bam file
    string scaffold;
    BamAlignment al;
    GffInterval results;
    string readName;
    string pacid;

    // Write the nonoverlapping sites to file.
    ofstream o(_nonoverlapping.c_str());

    // Cache for duplicate reads.
    map<string, map<string, int> > saw_read;

    while ( reader.GetNextAlignment(al) ) {

        ++_counts.total;

        readName = al.Name;
 
        // Clear the results vector
        results.clear(); 

        // BWA Sometimes does weird things with the reference dictionary, so
        // this check is needed to skip a possible segmentation fault.
        if (al.RefID >= 0) {
            scaffold = references[al.RefID].RefName;
            // Check if the scaffold is even in the map 
            if (_tree_list.count(scaffold) > 0) {
                // Check the strandedness and if the map contains intervals for that 
                // strand on this scaffold
                if (al.IsReverseStrand() && _tree_list[scaffold].count('-') > 0)
                    _tree_list[scaffold]['-'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else if (!al.IsReverseStrand() && _tree_list[scaffold].count('+') > 0)
                    _tree_list[scaffold]['+'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                    continue;
                }
                // Deal with dups 
                if (results.size() == 1) {
                    ++_counts.overlapping;
                    pacid = results[0].value.gene;
                    saw_read[readName][pacid]++;
                }

                else if (results.size() > 1) {
                    ++_counts.multiHitModel;
                }

                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                }
            }

            else {
                ++_counts.nonOverlapping;
                o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
            }
        }
    }
   
    default_random_engine generator;
    vector<string> keys;
    string choice;
    // Now iterate through the map, count duplicates, iterate through vector and add counts to _gene_counts
    // I don't know if the conversion to and from a set is efficient, but it still runs in about 105 seconds
    // with 6260336 reads on my computer.
    for (map<string, map<string, int> >::iterator i = saw_read.begin(); i != saw_read.end(); ++i) {
        if (i->second.size() == 1) {
            for (map<string, int>::iterator j = i->second.begin(); j != i->second.end(); ++j)
<<<<<<< HEAD
                _gene_counts[j->first] += j->second;
=======
                _gene_counts[j->first] += 1; 
>>>>>>> develop
        }
        else {
            ++_counts.dups;
            keys.clear(); 
            for (map<string, int>::iterator j = i->second.begin(); j != i->second.end(); ++j)
	        keys.push_back(j->first);

            uniform_int_distribution<int> distribution(0, keys.size() - 1);
            choice = keys[distribution(generator)];
            if (_annotateDups)
	        _dup_pac.insert(choice);
<<<<<<< HEAD
            _gene_counts[choice] += i->second[choice];
=======
            _gene_counts[choice] += 1; 
>>>>>>> develop
        }
    }
  
    o.close();
    reader.Close();
}

/**
 * Run GMCounts but split counts evenly across genes. 
 * Here we actually have to store PACIDs in a vector 
 */
void TagSeqGMCounts::ProcessBam_splitCounts() {
    // Open the BAM file
    BamReader reader;
    if ( !reader.Open(_infile)) {
        cerr << "Failed to open BAM file " << _infile << endl;
        exit(1);
    }

    // See BamTools Api 
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    // Variables for reading the Bam file
    string scaffold;
    BamAlignment al;
    GffInterval results;
    string readName;
    string pacid;

    // Write the nonoverlapping sites to file.
    ofstream o(_nonoverlapping.c_str());

    // Cache for duplicate reads.
    map<string, map<string, int> > saw_read;
    
    while ( reader.GetNextAlignment(al) ) {

        ++_counts.total;

        readName = al.Name;
 
        // Clear the results vector
        results.clear(); 

        // BWA Sometimes does weird things with the reference dictionary, so
        // this check is needed to skip a possible segmentation fault.
        if (al.RefID >= 0) {
            scaffold = references[al.RefID].RefName;
            // Check if the scaffold is even in the map 
            if (_tree_list.count(scaffold) > 0) {
                // Check the strandedness and if the map contains intervals for that 
                // strand on this scaffold
                if (al.IsReverseStrand() && _tree_list[scaffold].count('-') > 0)
                    _tree_list[scaffold]['-'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else if (!al.IsReverseStrand() && _tree_list[scaffold].count('+') > 0)
                    _tree_list[scaffold]['+'].findOverlapping(al.Position, al.GetEndPosition(), results);
                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                    continue;
                }
                // Deal with dups 
                if (results.size() == 1) {
                    ++_counts.overlapping;
                    pacid = results[0].value.gene;
                    saw_read[readName][pacid]++;
                }

                else if (results.size() > 1) {
                    ++_counts.multiHitModel;
                }

                else {
                    ++_counts.nonOverlapping;
                    o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
                }
            }

            else {
                ++_counts.nonOverlapping;
                o << scaffold << "\t" << al.Position << "\t" << al.GetEndPosition() << endl;
            }
        }
    }
  

    // Now iterate through the map, count duplicates, iterate through vector and split counts
    // I have to iterate through the map twice in the else statment... I'm sure there is a better way.
    for (map<string, map<string, int> >::iterator i = saw_read.begin(); i != saw_read.end(); ++i) {
        if (i->second.size() == 1) {
            for (map<string, int>::iterator j = i->second.begin(); j != i->second.end(); ++j)
                _gene_counts[j->first] += j->second;
        }
        else {
            ++_counts.dups;
            double sum_counts = 0.0;

            for (map<string, int>::iterator j = i->second.begin(); j != i->second.end(); ++j) 
                sum_counts += j->second;

	    double split = sum_counts / i->second.size();
	    long int rnd = lrint( split );
            for (map<string, int>::iterator j = i->second.begin(); j != i->second.end(); ++j) {
                if (_annotateDups)
                    _dup_pac.insert(j->first);
                _gene_counts[j->first] += rnd;
            }
        }
    }
  
    o.close();
    reader.Close();
}

/**
 * Writes the count data to the -o file. Sums the total number of expression counts.
 */
void TagSeqGMCounts::WriteCounts() {
    ofstream o_cts(_out.c_str());

    for (map<string, int>::const_iterator i = _gene_counts.begin(); i != _gene_counts.end(); ++i) {
        if (i->second > 0){
            if (_annotateDups && _dup_pac.count(i->first) > 0) {
                o_cts << i->first << "_DUP\t" << i->second << endl;
                _counts.totalPass += i->second;
            } else {
                o_cts << i->first << "\t" << i->second << endl;
                _counts.totalPass += i->second;
            }
        }
    }
 
    o_cts.close();
}

/**
 * Prints summary counts to stdout.
 */
void TagSeqGMCounts::PrintSummary() {
    cout     << "# Reads                = " << _counts.total          << endl;
    if (_primaryAln) 
        cout << "# Secondary alignments = " << _counts.notPrimaryAln  << endl;
    if (!_primaryAln) 
        cout << "# Duplicate reads      = " << _counts.dups           << endl;
    cout     << "# Gene-model hits      = " << _counts.overlapping    << endl;
    cout     << "# Gene-model multi     = " << _counts.multiHitModel  << endl;
    cout     << "# Gene-model misses    = " << _counts.nonOverlapping << endl;
    cout     << "# Final transcripts    = " << _counts.totalPass      << endl;
}
