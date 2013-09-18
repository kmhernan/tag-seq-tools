/***************************************
  CountMatrix.h

  2013 - Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
***************************************/
#ifndef __COUNT_MATRIX_H__
#define __COUNT_MATRIX_H__

#include <string>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <vector>
#include <cstring>
#include <sys/types.h>
#include <dirent.h>
#include <map>
#include <set>

typedef std::map<std::string, std::map<std::string, int> > ct_map;
// Class methods and elements
class TagSeqCountMatrix {

public:

    // constructor
    TagSeqCountMatrix(std::string &indir, std::string &out);

    // destructor
    ~TagSeqCountMatrix(void);

private:
    
    // IO files
    std::string _indir;
    std::string _out;
    std::vector<std::string> _files;
     
    // Count container
    ct_map _ct_map;
    std::set<std::string> _dup_pac; 
    // Header string
    std::string _header = "PACID";
    // Processing
    void Run();
    void LoadFiles();
    void ProcessFile(std::string&);
    void WriteCounts(); 
};
#endif // __COUNT_MATRIX_H__
