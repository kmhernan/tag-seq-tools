/***************************************
  StatsGeneModel.cpp

  2013 - Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
***************************************/

#include "StatsGeneModel.h"
using namespace BamTools;
using namespace std;

// build
TagSeqGMStats::TagSeqGMStats(string &infile, string &out, string &ingff) {

    _infile = infile;
    _out    = out;
    _ingff  = ingff;
    _counts = new Counts(); 
    _gff    = new GffFileReader(ingff);

    // Process
    Run();
}

// destroy
TagSeqGMStats::~TagSeqGMStats(void) {
    delete _gff;
    delete _counts;
}

/**
 * Wrapper function for running the application.
 */
void TagSeqGMStats::Run() {
    // load the GFF file
    cout << "Loading intervals..." << endl;
    int query = 0;
    _gff->LoadIntervals(query, _tree_list);

    // Run statistics
    cout << "Gathering statistics..." << endl;
    GetStatistics();

    // Print summary
    cout << "Writing summary..." << endl;
    PrintSummary();
}

void TagSeqGMStats::GetStatistics() {
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
    string pacid;

    // Duplicate read cache
    //map<string, map<string, int> > saw_read;
    map<string, int> saw_read;

    // PACID cache
    map<string, map<string, int> > plus_gene;
    map<string, map<string, int> > minus_gene;

    // Process the BAM records
    while ( reader.GetNextAlignment(al) ) {

        // Need statistics at both the READ level and the
        // PACID gene level.
        _counts->total++;

        results.clear();
      
        saw_read[al.Name]++;

        if (al.RefID >= 0) {
            scaffold = references[al.RefID].RefName;
            if (al.IsReverseStrand()) {
	        _counts->reverse++;
                if (_tree_list.count(scaffold) > 0 && _tree_list[scaffold].count('-') > 0) {
		    _counts->overlap++;
		    _counts->rev_overlap++;
                    _tree_list[scaffold]['-'].findOverlapping(al.Position, al.GetEndPosition(), results); 

                    // Check if there is only one gene in the interval
                    if (results.size() == 1) {
                        _counts->rev_one_gm++;
                        pacid = results[0].value.gene;
                        minus_gene[pacid][al.Name]++;
                    }
                    else
                        _counts->rev_multi_gm++;
                }
		else {
		    _counts->nonoverlap++;
		    _counts->rev_nonover++;
        	}
            }
            else {
                _counts->forward++;
                if (_tree_list.count(scaffold) > 0 && _tree_list[scaffold].count('+') > 0) {
		    _counts->overlap++;
		    _counts->frd_overlap++;
                    _tree_list[scaffold]['+'].findOverlapping(al.Position, al.GetEndPosition(), results); 

                    // Check if there is only one gene in the interval
                    if (results.size() == 1) {
                        _counts->frd_one_gm++;
                        pacid = results[0].value.gene;
                        plus_gene[pacid][al.Name]++;
                    }
                    else
                        _counts->frd_multi_gm++;
                }
	 	else {
		    _counts->nonoverlap++;
		    _counts->frd_nonover++;
 		}
            }
        }
    }
 
    reader.Close();

    // Iterate through map to count unique and duplicate reads
    for (map<string, int>::iterator i = saw_read.begin(); i != saw_read.end(); ++i) {
        if (i->second == 1)
            _counts->unique++;
        else
            _counts->duplicate++;
    }

    // Size of this map == number of PACIDS hit
    _counts-> rev_pac_hit = minus_gene.size();
    _counts-> frd_pac_hit = plus_gene.size();
    _counts-> total_pac   = _counts->rev_pac_hit + _counts->frd_pac_hit;

    // Iterate through pac map to count the number of unique and duplicate reads
    for (map<string, map<string, int> >::iterator j = minus_gene.begin(); j != minus_gene.end(); j++){
	bool is_dup = false;
        for(map<string, int>::iterator k = j->second.begin(); k != j->second.end(); k++) {
            if (saw_read[k->first] > 1) {
                is_dup = true;
                break;
            }
        }
        if (is_dup)
            _counts->rev_dup_pac++;
        else
            _counts->rev_unique_pac++;
    }

    for (map<string, map<string, int> >::iterator j = plus_gene.begin(); j != plus_gene.end(); j++){
	bool is_dup = false;
        for(map<string, int>::iterator k = j->second.begin(); k != j->second.end(); k++) {
            if (saw_read[k->first] > 1) {
                is_dup = true;
                break;
            }
        }
        if (is_dup)
            _counts->frd_dup_pac++;
        else
            _counts->frd_unique_pac++;
    }
}

/**
 * Write summary statistics to file
 */
void TagSeqGMStats::PrintSummary() {
    ofstream o(_out.c_str());
    double c = _counts->total;

    o << "StatLevel\tStat\tCount\tPercentage"       << endl;
    o << "Read\t" << "Total\t"                      << _counts->total << "\t1.000000" << endl; 
    o << "Read\t" << "Reverse\t"                    << _counts->reverse          << "\t" 
						       << _counts->reverse / c      << endl; 
    o << "Read\t" << "Forward\t"                    << _counts->forward          << "\t" 
						       << _counts->forward / c      << endl; 
    o << "Read\t" << "Unique\t"                     << _counts->unique           << "\t" 
						       << _counts->unique / c       << endl; 
    o << "Read\t" << "Duplicate\t"                  << _counts->duplicate        << "\t" 
						       << _counts->duplicate / c    << endl;
    o << "Read\t" << "GM Overlaps\t"                << _counts->overlap          << "\t" 
						       << _counts->overlap / c      << endl;
    o << "Read\t" << "GM Reverse Overlaps\t"        << _counts->rev_overlap      << "\t" 
                                                       << _counts->rev_overlap / c  << endl;
    o << "Read\t" << "GM Forward Overlaps\t"        << _counts->frd_overlap      << "\t" 
                   			               << _counts->frd_overlap / c  << endl;
    o << "Read\t" << "GM Nonoverlaps\t"             << _counts->nonoverlap       << "\t" 
				                       << _counts->nonoverlap / c   << endl;
    o << "Read\t" << "GM Reverse Nonoverlap\t"      << _counts->rev_nonover      << "\t" 
                                                       << _counts->rev_nonover / c  << endl;
    o << "Read\t" << "GM Forward Nonoverlap\t"      << _counts->frd_nonover      << "\t" 
					               << _counts->frd_nonover / c  << endl;
    o << "Read\t" << "Reverse Single GM Overlaps\t" << _counts->rev_one_gm       << "\t" 
						       << _counts->rev_one_gm / c   << endl;
    o << "Read\t" << "Reverse Multi GM Overlaps\t"  << _counts->rev_multi_gm     << "\t" 
						       << _counts->rev_multi_gm / c << endl;
    o << "Read\t" << "Forward Single GM Overlaps\t" << _counts->frd_one_gm       << "\t" 
						       << _counts->frd_one_gm / c   << endl;
    o << "Read\t" << "Forward Multi GM Overlaps\t"  << _counts->frd_multi_gm     << "\t" 
					               << _counts->frd_multi_gm / c << endl;
    // Now we print out the PACID level stats
    c = _counts->total_pac;
    o << "PACID\t" << "Total PACID Hits\t"          << _counts->total_pac << "\t1.000000" << endl;
    o << "PACID\t" << "Reverse PACID Hits\t"        << _counts->rev_pac_hit         << "\t"
						       << _counts->rev_pac_hit / c     << endl;
    o << "PACID\t" << "Reverse PACID Duplicates\t"  << _counts->rev_dup_pac         << "\t"
						       << _counts->rev_dup_pac / c     << endl;
    o << "PACID\t" << "Reverse PACID Unique\t"      << _counts->rev_unique_pac      << "\t"
						       << _counts->rev_unique_pac / c  << endl;
    o << "PACID\t" << "Forward PACID Hits\t"        << _counts->frd_pac_hit         << "\t"
						       << _counts->frd_pac_hit / c     << endl;
    o << "PACID\t" << "Forward PACID Duplicates\t"  << _counts->frd_dup_pac         << "\t"
						       << _counts->frd_dup_pac / c     << endl;
    o << "PACID\t" << "Forward PACID Unique\t"      << _counts->frd_unique_pac      << "\t"
						       << _counts->frd_unique_pac / c  << endl;

    // Close stream
    o.close();
}
