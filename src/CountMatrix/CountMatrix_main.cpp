/***************************************
  CountMatrix_main.cpp

  2013 - Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
***************************************/
#include "CountMatrix.h"

using namespace std;

// Define parameter checking macro as used in bedtools
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function prototypes
void countMatrix_help(void);

int countMatrix_main(int argc, char* argv[]) {

    // configuration variables
    bool showHelp = false;

    // I/O files
    string indir; 
    string out;

    // parm flags
    bool haveInput = false;
    bool haveOut   = false;

    if (argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 6, parameterLength))) {
          showHelp = true;
        }
    }

    if(showHelp) countMatrix_help(); 

    // parse command line
    for (int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);
         
        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveInput = true;
                indir = argv[i + 1];
                i++; 
            }
      }

      else if (PARAMETER_CHECK("-o", 2, parameterLength)) {
        if ((i+1) < argc) {
          haveOut = true;
          out = argv[i + 1];
          i++;
        }
      }

      else {
        cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << endl;
        showHelp = true;
      }
    }

    // Check required elements 
    if (!haveInput) {
      cerr << "*****ERROR: Need -i input directory of GMCounts files." << endl;
      showHelp = true; 
    }
    else if (!haveOut) {
      cerr << "*****ERROR: Need -o output expression matrix file." << endl;
      showHelp = true; 
    }

    if (!showHelp) {
      TagSeqCountMatrix *ts  = new TagSeqCountMatrix(indir, out);
      delete ts;
    }

    else {
      countMatrix_help();
    }

    return 0;
}

void countMatrix_help(void)
{
    cout << endl;
    cout << "TagSeqTools CountMatrix: Create expression count matrix by merging "                 << endl;
    cout << "                         count files created with GMCounts into one file."           << endl;
    cout << "Kyle Hernandez, 2013. UNLICENSED <www.unlicense.org>"                                << endl;
    cout << endl; 
    cout << "USAGE:TagSeqTools CountMatrix -i <in-directory> -o <output.tab>"			  << endl; 
    cout << endl; 

    cout << "Required Arguments:"                                                                 << endl;
    cout << "Input either a directory with files:\n"                                              << endl;
    cout << "    -i                   " << "Input parent directory of all the GMCount files you " << endl;
    cout << "                         " << "want to combine.\n"                                         << endl;
    cout << "    -o                   " << "Output expression matrix file"                        << endl;
    cout << endl;

    cout << "Optional Arguments:"                                                                 << endl;
    cout << "    -h                   " << "Print this message and exit"                          << endl;
    cout << endl;

    exit(1);
}
