#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <vector>
#include <cstdlib>
#include <climits>
#include <algorithm>

#include <getopt.h>
#include <errno.h>
#include <err.h>

#include "config.h"

using namespace std;

size_t reads_count=0;
size_t sequences_count=0;
enum INPUT_TYPE {
	INPUT_AUTO_DETECT=1,
	INPUT_FASTA,
	INPUT_FASTQ,
};
INPUT_TYPE input_type=INPUT_AUTO_DETECT;
enum COLLAPSED_MODE {
	COLLAPSED_AUTO_DETECT=1,
	COLLAPSED_READS,
	NO_COLLAPSED_READS
};
COLLAPSED_MODE collapse_mode=COLLAPSED_AUTO_DETECT;
bool verbose;
vector<string> input_files;

size_t current_line;
string current_file_name;
size_t lines_in_current_file;
size_t reads_in_current_file;
size_t sequences_in_current_file;


void show_help()
{
	cout << "HELP!!!!!" << endl;

}

enum OPTS {
	OPT_INPUT_FASTA = CHAR_MAX+1,
	OPT_INPUT_FASTQ,
	OPT_COLLAPSED,
	OPT_NO_COLLAPSED,
};

struct option long_options[] = {
	{"fasta",       no_argument, 0, OPT_INPUT_FASTA},
	{"fastq",       no_argument, 0, OPT_INPUT_FASTQ},
	{"collapsed",   no_argument, 0, OPT_COLLAPSED},
	{"nocollapsed", no_argument, 0, OPT_NO_COLLAPSED},
	{"help",        no_argument, 0, 'h'},
	{"verbose",     no_argument, 0, 'v'},
	{NULL,0,0,0}
};

void parse_commandline_options(int argc, char* argv[])
{
	int c;
	int long_opt_index;
	while ( (c=getopt_long(argc,argv,"h",long_options,&long_opt_index)) != -1 ) {
		switch(c)
		{
		case OPT_INPUT_FASTA:
			input_type = INPUT_FASTA;
			break;
		case OPT_INPUT_FASTQ:
			input_type = INPUT_FASTQ;
			break;
		case OPT_COLLAPSED:
			collapse_mode = COLLAPSED_READS;
			break;
		case OPT_NO_COLLAPSED:
			collapse_mode = NO_COLLAPSED_READS;
			break;
		case 'h':
			show_help();
			exit(0);
			break;
		case 'v':
			verbose=true;
			break;
		default:
			errx(1,"Error: invalid command-line option");
		}
	}

	for (int i=optind; i<argc; ++i)
		input_files.push_back(argv[i]);

	if (input_files.empty())
		errx(1,"Error: missing input file names. Use 'stdin' to read from STDIN. See --help for more details");
}

// returns TRUE
// if the input string looks like a valid collapsed read ID:
//  NNN-NNNN
inline bool is_collapsed_id(const std::string &str, size_t& /*output*/ read_count)
{
	errno = 0 ;
	char* endptr=NULL;

	//Try the first number
	const char* nptr = str.c_str();
	size_t result = strtoul(nptr, &endptr, 10);
	if (endptr == nptr)
		return false; //no digits found

	if (errno != 0)
		return false; //another error occured (ERANGE, usually)

	if (*endptr !='-')
		return false; //we expact a minus/dash after the first number

	nptr = endptr+1;
	result = strtoul(nptr, &endptr, 10);

	if (endptr == nptr)
		return false; //no digits found

	if (errno != 0)
		return false; //another error occured (ERANGE, usually)

	if (*endptr != 0)
		return false; //we expect the end of the string here

	read_count = result;

	return true;
}

void count_reads_from_id(const std::string& id)
{
	if (id.empty())
		errx(1,"Input error: file '%s' line '%zu' is empty", current_file_name.c_str(), current_line);

	if (input_type==INPUT_FASTA && id[0] != '>')
		errx(1,"Input error: file '%s' line '%zu': expecting FASTA file with '>', got '%s'",
				current_file_name.c_str(), current_line, id.c_str());

	if (input_type==INPUT_FASTQ && id[0] != '@')
		errx(1,"Input error: file '%s' line '%zu': expecting FASTQ file with '@', got '%s'",
				current_file_name.c_str(), current_line, id.c_str());

	// skip the  '>' / '@' prefix
	const std::string &id2 ( id.substr(1) ) ;

	size_t count=0;

	if (collapse_mode == COLLAPSED_AUTO_DETECT) {
		// check if this looks like a collapsed read ID
		// must match /^\d+-\d+$/
		if (is_collapsed_id(id2,count)) {
			collapse_mode = COLLAPSED_READS;
			if (verbose)
				cerr << "Detected collapsed READ-IDs in '" << current_file_name << "'" << endl;
		} else {
			collapse_mode = NO_COLLAPSED_READS;
			if (verbose)
				cerr << "Detected non-collapsed READ-IDs in '" << current_file_name << "'" << endl;
		}
	}


	switch(collapse_mode)
	{
	case COLLAPSED_READS:
		if (!is_collapsed_id(id2,count))
			errx(1,"Input error: file '%s' line '%zu': expecting collapsed read-id, got '%s'",
				current_file_name.c_str(), current_line, id.c_str());
		reads_in_current_file+=count;
		reads_count+=count;
		sequences_in_current_file++;
		sequences_count++;
		break;

	case NO_COLLAPSED_READS:
		reads_in_current_file++;
		sequences_in_current_file++;
		reads_count++;
		sequences_count++;
		break;

	case COLLAPSED_AUTO_DETECT:
	default:
		errx(1,"Internal error: collapsed mdoe not detected (file = '%s' line %zu)",
				current_file_name.c_str(), current_line);
	}
}

inline bool read_one_line(istream& strm, std::string& /*output*/ line)
{
	getline(strm, line);

	if (strm.eof())
		return false;

	if (strm.bad() || strm.fail()) {
		err(1,"Error: failed to read from file '%s'", current_file_name.c_str());
	}

	++current_line;
	++lines_in_current_file;

	return true;
}


// Consume and ignore linse from the input stream
// (but verify the lines were read successfully)
void consume_and_discard_lines(istream &strm, size_t count)
{
	string s;
	for (size_t i=0;i<count;++i) {
		if (!read_one_line(strm,s))
			errx(1,"Error: file '%s' line '%zu': got empty line in the middle of a read (file truncated?)",
				current_file_name.c_str(), current_line);
	}
}

bool detect_file_type(const std::string& first_line)
{
	// Detect type based on first line
	if (first_line.length()==0)
		errx(1,"Error: file '%s' has empty first line (can't detect FASTA or FASQT)",
				current_file_name.c_str());

	switch(first_line[0])
	{
	case '>':
		if (verbose)
			cerr << "type auto-detection: FASTA" << endl;
		input_type = INPUT_FASTA;
		break;

	case '@':
		if (verbose)
			cerr << "type auto-detection: FASTQ" << endl;
		input_type = INPUT_FASTQ;
		break;

	default:
		errx(1,"Error: file '%s' doesn't look like FASTA or FASTQ, first line = \n%s\n",
				current_file_name.c_str(), first_line.c_str());
	}


	return true;
}

void count_reads_in_stream(istream &strm)
{
	while (1) {
		string s;
		if (!read_one_line(strm,s)) {
			//can happen if the file is empty
			if (current_line==0 && verbose)
				cerr << "file '" << current_file_name << "' is empty" << endl;
			return;
		}

		if (input_type==INPUT_AUTO_DETECT)
			detect_file_type(s);

		//Count the reads of this first line
		count_reads_from_id(s);

		//Continue to consume the reads of the lines for this read
		consume_and_discard_lines(strm, (input_type==INPUT_FASTA)?1:3);
	}
}

void count_reads_in_file(const std::string& filename)
{
	//TODO:
	// If the file is a regular file (not a pipe),
	// try to detect compressed files (gzip/bzip2) and automatically decompress
	// using "gzstream" or "bz2stream" classes
	current_line = 0 ;
	lines_in_current_file = 0 ;
	reads_in_current_file = 0 ;
	sequences_in_current_file = 0 ;
	current_file_name = filename;

	if (verbose)
		cerr << "Reading from '" << filename << "'..." << endl;

	if (filename == "stdin" || filename == "-") {
		count_reads_in_stream(cin);
	} else {
		ifstream f(filename.c_str());
		if (!f)
			err(1,"Error: failed to open file '%s'", filename.c_str());
		count_reads_in_stream(f);
	}

	if (verbose)
		cerr << "File '" << filename << "': " << lines_in_current_file << " lines, "
			<< reads_in_current_file << " reads" << endl;
}


int main(int argc, char* argv [])
{
	//Chapter 13, Section 13.14.1 "Synchronization with C's standard streams" page 682
	std::ios_base::sync_with_stdio(false);

	//Chapter 13, Section 13.10.1 "Loose Coupling using tie()" page 638
	std::cin.tie(static_cast<std::ostream*>(0));

	parse_commandline_options(argc,argv);
	for_each(input_files.begin(),input_files.end(),count_reads_in_file);

	//reporting
	cout << reads_count;
	return 0;
}
