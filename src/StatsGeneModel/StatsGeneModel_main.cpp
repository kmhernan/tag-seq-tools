/***************************************
  StatsGeneModel_main.cpp

  2013 - Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
***************************************/
#include "StatsGeneModel.h"

using namespace std;

// Define parameter checking macro as used in bedtools
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function prototypes
void statsGeneModel_help(void);

int statsGeneModel_main(int argc, char* argv[]) {

    // configuration variables
    bool showHelp = false;

    // I/O files
    string inbam; 
    string out;
    string gff;

    // parm flags
    bool haveInput = false;
    bool haveOut   = false;
    bool haveGff   = false;

    if (argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 6, parameterLength))) {
          showHelp = true;
        }
    }

    if(showHelp) statsGeneModel_help(); 

    // parse command line
    for (int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);
         
        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveInput = true;
                inbam = argv[i + 1];
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

      else if (PARAMETER_CHECK("-g", 2, parameterLength)) {
        if ((i+1) < argc) {
          haveGff = true;
          gff = argv[i + 1];
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
      cerr << "*****ERROR: Need -i input BAM file." << endl;
      showHelp = true; 
    }
    else if (!haveOut) {
      cerr << "*****ERROR: Need -o output stats file." << endl;
      showHelp = true; 
    }
    else if (!haveGff) {
      cerr << "*****ERROR: Need -g gff file." << endl;
      showHelp = true; 
    }

    if (!showHelp) {
      TagSeqGMStats *ts  = new TagSeqGMStats(inbam, out, gff);
      delete ts;
    }

    else {
      statsGeneModel_help();
    }

    return 0;
}

void statsGeneModel_help(void)
{
    cout << endl;
    cout << "TagSeqTools GMStats: Summary statistics for transcripts mapped to a reference "      << endl;
    cout << "                     genome and overlapped with gene-model files."                   << endl;
    cout << "Kyle Hernandez, 2013. UNLICENSED <www.unlicense.org>"                                << endl;
    cout << endl; 
    cout << "USAGE:TagSeqTools GMStats -i <file.bam> -g <file.gff> -o <output.tab>"      	  << endl; 
    cout << endl; 

    cout << "Required Arguments:"                                                                 << endl;
    cout << "    -i                   " << "Input BAM file. "					  << endl;
    cout << "    -g                   " << "Input GFF file. "                                     << endl;
    cout << "    -o                   " << "Output tab-delimited statistics file. "               << endl;
    cout << endl;

    cout << "Optional Arguments:"                                                                 << endl;
    cout << "    -h                   " << "Print this message and exit"                          << endl;
    cout << endl;

    exit(1);
}
