/*******************************************************************************
  TagSeqTools.cpp 

  Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
*******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <time.h>
using namespace std;

// Define parameter checker as used in bedtools
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

int countGeneModel_main(int argc, char* argv[]);
int countMatrix_main(int argc, char* argv[]);
<<<<<<< HEAD
=======
int statsGeneModel_main(int argc, char* argv[]);
>>>>>>> develop

int tagseq_help(void);

int main(int argc, char *argv[]) {

    // Status of exit 
    int STATUS = 0; 

    // Timer
    clock_t t;
    t = clock();
    
    // Make sure ther user has entered a sub_command
    if (argc < 2) { return tagseq_help(); }

    string sub_cmd = argv[1];

    // help
    if (sub_cmd == "-h" || sub_cmd == "--help" || sub_cmd == "-help")
        return tagseq_help();

    // tools
    else if (sub_cmd == "GMCounts") { STATUS = countGeneModel_main(argc-1, argv+1); }
    else if (sub_cmd == "CountMatrix") { STATUS = countMatrix_main(argc-1, argv+1); }
<<<<<<< HEAD
=======
    else if (sub_cmd == "GMStats") { STATUS = statsGeneModel_main(argc-1, argv+1); }
>>>>>>> develop

    // Unknown command
    else {
        cerr << "error: unrecognized command: " << argv[1] << endl << endl;
        return 1;
    }

    t = clock() - t;
    printf ("Finished; Took %f seconds.\n", ((float)t)/CLOCKS_PER_SEC);

    return STATUS;
}

int tagseq_help(void) {
    cout << "TagSeqTools: multiple tools for 3' tagSeq data.\n";
    cout << "usage:     TagSeqTools <subcommand> [options]" << endl << endl;
    cout << "The sub-commands include:" << endl;

    cout << endl;
    cout << "    GMCounts    " << "Create expression counts from reads mapped to genome [GFF3 required].\n";
    cout << "    CountMatrix " << "Merge GMCounts files to create an expression count matrix.\n";
<<<<<<< HEAD
=======
    cout << "    GMStats     " << "Run summary statistics on reads mapped to genome and overlap \n";
    cout << "                " << "mRNA intervals in a GFF file.\n";
>>>>>>> develop
    cout << endl;
    return 0;
}
