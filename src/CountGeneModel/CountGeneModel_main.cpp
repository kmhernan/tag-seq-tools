/***************************************
  CountGeneModel_main.cpp

  2013 - Kyle Hernandez
  Juenger Laboratory
  Section of Integrative Biology
  The University of Texas at Austin
  kmhernan84@gmail.com

  Unlicensed, <http://unlicense.org/>
***************************************/
#include "CountGeneModel.h"

using namespace std;

// Define parameter checking macro as used in bedtools
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function prototypes
void countGeneModel_help(void);

int countGeneModel_main(int argc, char* argv[]) {

    // configuratin variables
    bool showHelp = false;

    // I/O files
    string infile;
    string ingff;
    string out;
    string nonoverlapping;

    // Gff query choice
    int query = 0;

    // parm flags
    bool haveInput = false;
    bool haveGff = false;
    bool haveOut = false;
    bool haveNon = false;
    // Duplicate reads options
    bool noDups      = true;
    bool primaryAln  = false;
    bool randomOne   = false; 
    bool splitCounts = false;
    bool annotateDups = false;

    if (argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 6, parameterLength))) {
          showHelp = true;
        }
    }

    if(showHelp) countGeneModel_help(); 

    // parse command line
    for (int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);
         
        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveInput = true;
                infile = argv[i + 1];
                i++; 
            }
      }

      else if (PARAMETER_CHECK("-g", 2, parameterLength)) {
        if ((i+1) < argc) {
          haveGff = true;
          ingff = argv[i + 1];
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

      else if (PARAMETER_CHECK("-n", 2, parameterLength)) {
        if ((i+1) < argc) {
          haveNon = true;
          nonoverlapping = argv[i + 1];
          i++;
        }
      }

      else if (PARAMETER_CHECK("--gene", 6, parameterLength)) {
        if ((i+1) < argc) {
          query = 1;
          i++;
        }
      }

      else if (PARAMETER_CHECK("--mRNA", 6, parameterLength)) {
        if ((i+1) < argc) {
          query = 0;
          i++;
        }
      }

      // Options for duplicately mapped reads
      else if (PARAMETER_CHECK("--primary-alignment", 19, parameterLength)) {
          primaryAln = true;
          noDups = false; 
      }

      else if (PARAMETER_CHECK("--random-one", 12, parameterLength)) {
          randomOne = true; 
          noDups = false; 
      }

      else if (PARAMETER_CHECK("--split-counts", 14, parameterLength)) {
          splitCounts = true; 
          noDups = false; 
      }

      else if (PARAMETER_CHECK("--flag-dups", 11, parameterLength)) {
          annotateDups = true; 
      }

      else {
        cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << endl;
        showHelp = true;
      }
    }

    // Check required elements 
    if (!haveInput) {
      cerr << "*****ERROR: Need -i input BAM." << endl;
      showHelp = true; 
    }
    else if (!haveGff) {
      cerr << "*****ERROR: Need -g input GFF3 file." << endl;
      showHelp = true; 
    }
    else if (!haveOut) {
      cerr << "*****ERROR: Need -o output expression file." << endl;
      showHelp = true; 
    }
    else if (!haveNon) {
      cerr << "*****ERROR: Need -n non-overlapping transcripts file." << endl;
      showHelp = true; 
    }
    // Check mutally exclusive dups options
    else if (primaryAln && (randomOne || splitCounts)) {
      cerr << "*****ERROR: You can only choose one duplicate read function." << endl;
      showHelp = true; 
    }
    else if (randomOne && (primaryAln || splitCounts)) {
      cerr << "*****ERROR: You can only choose one duplicate read function." << endl;
      showHelp = true; 
    }
    else if (splitCounts && (primaryAln || randomOne)){ 
      cerr << "*****ERROR: You can only choose one duplicate read function." << endl;
      showHelp = true; 
    }
    else if (noDups && annotateDups) {
      cerr << "*****ERROR: No duplicate reads will be processed, so they can't be annotated." << endl;
      showHelp = true;
    }
    else if (primaryAln && annotateDups) {
      cerr << "*****ERROR: --flag-dups is not compatible with --primary-alignment" << endl;
      showHelp = true;
    }

    // Run the application
    if (!showHelp) {
      TagSeqGMCounts *ts = new TagSeqGMCounts(infile, ingff, out, nonoverlapping, query,
					      noDups, primaryAln, randomOne, splitCounts,
					      annotateDups);
      delete ts;
    }

    else {
      countGeneModel_help();
    }

    return 0;
}

void countGeneModel_help(void)
{
    cout << endl;
    cout << "TagSeqTools GMCounts: Create transcript counts from reads mapped to a genome "       << endl;
    cout << "                      reference using a GFF annotation file for appropriate  "       << endl;
    cout << "                      intervals based on user input arguments."                      << endl;
    cout << "Kyle Hernandez, 2013. UNLICENSED <www.unlicense.org>"                                << endl;
    cout << endl;
 
    cout << "USAGE: TagSeqTools GMCounts -i <file.bam> -g <file.gff> -o <file.tab> "              << endl;
    cout << "                            -n <nonoverlapping.tab> [--mRNA | --gene] "              << endl;
    cout << "                            [--primary-alignment | --random-one | --split-counts]"   << endl;
    cout << "                            [--flag-dups]"                 			  << endl;
    cout << endl; 

    cout << "Required Arguments:"                                                                 << endl;
    cout << "    -i                   " << "Input bam file. Requires strand flags"                << endl;
    cout << "    -g                   " << "<String>\tInput Gff file"                             << endl;
    cout << "    -o                   " << "Output expression matrix file"                        << endl;
    cout << "    -n                   " << "Output nonoverlapping reads file"                     << endl;
    cout << endl;

    cout << "Duplicate read handling (Removing duplicate reads is the default):"                  << endl;
    cout << "    --primary-alignment  " << "Use only primary-alignments. Requires that the "      << endl;
    cout << "                               BAM files have such information (e.g., bwa mem -M)"   << endl;
    cout << "    --random-one         " << "Pick a random read."                                  << endl;
    cout << "    --split-counts       " << "Distribute counts across all targets evenly."         << endl;
    cout << "    --flag-dups          " << "Append '_DUP' to PACIDs with duplicate reads."        << endl;
    cout << endl;

    cout << "Optional Arguments:"                                                                 << endl;
    cout << "    --gene               " << "Query intervals on 'gene' in GFF [mRNA is default]"   << endl;
    cout << "    -h                   " << "Print this message and exit"                          << endl;
    cout << endl;

    exit(1);
}
