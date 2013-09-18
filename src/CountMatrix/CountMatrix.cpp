/***************************************
  CountMatrix.cpp

  2013 - Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
***************************************/
#include "CountMatrix.h" 
using namespace std;

//build
TagSeqCountMatrix::TagSeqCountMatrix(string& indir, string& out) {
    _indir = indir; 
    _out   = out;
    
    // Process
    Run();
}

// destroy
TagSeqCountMatrix::~TagSeqCountMatrix(void) { }

/**
 * Wrapper function for running the application
 */
void TagSeqCountMatrix::Run() {
    LoadFiles();

    for (vector<string>::iterator i = _files.begin(); i != _files.end(); ++i) {
        ProcessFile(*i);
    }

    WriteCounts();
}

/**
 * Read file names into vector
 */
void TagSeqCountMatrix::LoadFiles() {
    DIR * curr;
    struct dirent *dp;
    curr = opendir(_indir.c_str());
    while ((dp = readdir(curr)) != NULL) {
        if (strcmp(dp->d_name, ".")) {
            if(strcmp(dp->d_name, "..")) 
                _files.push_back(dp->d_name);
        }
    }
    (void)closedir(curr);
}

/**
 * Build map
 */
void TagSeqCountMatrix::ProcessFile(string &fil) {
    string full   = _indir + fil;
    size_t ext    = fil.find_last_of('.');
    string sample = fil.substr(0, ext);
    ifstream f(full.c_str());
    string line;
    size_t found;
    size_t dup;
    string pacid;
    int cts;
    
    _header += "\t" + sample;

    while ( getline(f, line) ) {
        found = line.find('\t');
        pacid = line.substr(0, found);
        dup = pacid.find('_');
        if (dup != string::npos) {
            pacid = line.substr(0, dup);
            _dup_pac.insert(pacid);
        }
        cts = atoi(line.substr(found, line.size()).c_str());
        _ct_map[pacid][sample] = cts;
    }

    f.close();
}

/**
 * Write to file.
 */
void TagSeqCountMatrix::WriteCounts() {
    ofstream o(_out.c_str());
    size_t ext;
    string sample;
    string curr;

    o << _header << endl;
    for (ct_map::const_iterator j = _ct_map.begin(); j != _ct_map.end(); ++j) {
        if (_dup_pac.count(j->first) > 0)
            o << j->first << "_DUP";
        else 
            o << j->first;

        for (vector<string>::iterator i = _files.begin(); i != _files.end(); ++i) {
            curr   = *i;
            ext    = curr.find_last_of('.');
            sample = curr.substr(0, ext);
            o << "\t" << _ct_map[j->first][sample]; 
        }
        o << "\n";
    }
    o.close(); 
}
