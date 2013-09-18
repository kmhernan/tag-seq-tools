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
#include "GffFileReader.h"
using namespace std;

// build
GffFileReader::GffFileReader(string& inGff) { 
  gffFile = inGff; 
}

// destroy
GffFileReader::~GffFileReader(void) { }

// Load the interval tree
void GffFileReader::LoadIntervals(int& query, IntervalTreeMap& tree_list) {
    Open();
    string line;
    string currScaff;
    Gff record;
    GffInterval forward_intervals;
    GffInterval reversed_intervals;
    IntervalTree<Gff> tree_f;
    IntervalTree<Gff> tree_r;
    string read_id;
    string queryString = ((query == 0) ? "mRNA" : "gene");

    while (getline( *_gffStream, line )) {

        if (line[0] == '#') continue;

        else {
            GetRecord(line, record);

            if (record.type == queryString){

                if (currScaff.empty()){
                    currScaff = record.scaffold;
                    if (record.ori == '+') 
                        forward_intervals.push_back(Interval<Gff>(record.start, record.end, record));
                    else
                        reversed_intervals.push_back(Interval<Gff>(record.start, record.end, record));
                }

                else if (record.scaffold != currScaff){
                    tree_f = IntervalTree<Gff>(forward_intervals);
                    tree_r = IntervalTree<Gff>(reversed_intervals);
                    tree_list[currScaff]['+'] = tree_f;
                    tree_list[currScaff]['-'] = tree_r;
                    currScaff = record.scaffold;
                    forward_intervals.clear();
                    reversed_intervals.clear();
                    if (record.ori == '+') 
                        forward_intervals.push_back(Interval<Gff>(record.start, record.end, record));
                    else
                        reversed_intervals.push_back(Interval<Gff>(record.start, record.end, record));
                }

                else {
                    if (record.ori == '+') 
                        forward_intervals.push_back(Interval<Gff>(record.start, record.end, record));
                    else
                        reversed_intervals.push_back(Interval<Gff>(record.start, record.end, record));
                }
            }
        }
    }

    tree_f = IntervalTree<Gff>(forward_intervals);
    tree_r = IntervalTree<Gff>(reversed_intervals);
    tree_list[currScaff]['+'] = tree_f;
    tree_list[currScaff]['-'] = tree_r;
    Close();
}

Gff& GffFileReader::GetRecord(string& line, Gff& record){
  stringstream ss( line );
  int i = 0;
  string field;

  while (getline( ss, field, '\t'))
  {
    if (i == 0)
      record.scaffold = field;
    else if (i == 2)
      record.type = field;
    else if (i == 3)
      record.start = atoi(field.c_str());
    else if (i == 4)
      record.end = atoi(field.c_str());
    else if (i == 6)
      record.ori = field[0];
    else if (i == 8)
    {
      size_t found = field.find(';');
      if (found != string::npos)
        record.gene = field.substr(3, found-3);
    }
    i++;
  }

  return record;
}

// Open the Gff file
void GffFileReader::Open(void) {
  _gffStream = new ifstream(gffFile.c_str(), ios::in);
  if ( _gffStream->fail() ) {
      cerr << "*****Error: The requested file ("
           << gffFile
           << ") "
           << "could not be opened. "
           << "Error message: ("
           << strerror(errno)
           << "). Exiting!" << endl;
      exit(1);
  }
}

// Close the Gff file
void GffFileReader::Close() {
  delete _gffStream;
}
