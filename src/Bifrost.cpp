#include "CompactedDBG.hpp"
#include "ColoredCDBG.hpp"
#include <iomanip>

using namespace std;

void PrintVersion() {

    cout << BFG_VERSION << endl;
}

void PrintUsage() {

    cout << "Bifrost " << BFG_VERSION << endl << endl;

    cout << "Highly parallel construction, indexing and querying of colored and compacted de Bruijn graphs" << endl << endl;

    cout << "Usage: Bifrost [COMMAND] [PARAMETERS]" << endl << endl;

    cout << "[COMMAND]:" << endl << endl;

    cout << "   build                   Build a compacted de Bruijn graph, with or without colors" << endl;
    cout << "   update                  Update a compacted (possible colored) de Bruijn graph with new sequences" << endl;
    cout << "   query                   Query a compacted (possible colored) de Bruijn graph" << endl << endl;

    cout << "[PARAMETERS]: build" << endl << endl;

    cout << "   > Mandatory with required argument:" << endl << endl;

    cout << "   -s, --input-seq-file     Input sequence file (FASTA/FASTQ possibly gzipped)" << endl;
    cout << "                            Multiple files can be provided as a list in a TXT file (one file per line)" << endl;
    cout << "                            K-mers with exactly 1 occurrence in the input sequence files will be discarded" << endl;
    cout << "   -r, --input-ref-file     Input reference file (FASTA/FASTQ possibly gzipped and GFA)" << endl;
    cout << "                            Multiple files can be provided as a list in a TXT file (one file per line)" << endl;
    cout << "                            All k-mers of the input reference files are used" << endl;
    cout << "   -o, --output-file        Prefix for output file(s)" << endl << endl;

    cout << "   > Optional with required argument:" << endl << endl;

    cout << "   -t, --threads            Number of threads (default is 1)" << endl;
    cout << "   -k, --kmer-length        Length of k-mers (default is 31)" << endl;
    cout << "   -m, --min-length         Length of minimizers (default is automatically chosen)" << endl;
    cout << "   -b, --bloom-bits         Number of Bloom filter bits per k-mer with 1+ occurrences in the input files (default is 14)" << endl;
    cout << "   -B, --bloom-bits2        Number of Bloom filter bits per k-mer with 2+ occurrences in the input files (default is 14)" << endl;
    cout << "   -l, --load-mbbf          Input Blocked Bloom Filter file, skips filtering step (default is no input)" << endl;
    cout << "   -w, --write-mbbf         Output Blocked Bloom Filter file (default is no output)" << endl;
    cout << "   -u, --chunk-size         Read chunk size per thread (default is 64)" << endl;
    cout << "   -C, --compress-tigs      Compress the tigs with the given algorithm (default is no compression)" << endl << endl;

    cout << "   > Optional with no argument:" << endl << endl;

    cout << "   -c, --colors             Color the compacted de Bruijn graph (default is no coloring)" << endl;
    cout << "   -y, --keep-mercy         Keep low coverage k-mers connecting tips" << endl;
    cout << "   -i, --clip-tips          Clip tips shorter than k k-mers in length" << endl;
    cout << "   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length" << endl;
    cout << "   -a, --fasta              Output file is in FASTA format (only sequences) instead of GFA" << endl;
    cout << "   -v, --verbose            Print information messages during execution" << endl << endl;

    cout << "[PARAMETERS]: update" << endl << endl;

    cout << "  > Mandatory with required argument:" << endl << endl;

    cout << "   -g, --input-graph-file   Input graph file to update (GFA format)" << endl;
    cout << "   -s, --input-seq-file     Input sequence file (FASTA/FASTQ possibly gzipped)" << endl;
    cout << "                            Multiple files can be provided as a list in a TXT file (one file per line)" << endl;
    cout << "                            K-mers with exactly 1 occurrence in the input sequence files will be discarded" << endl;
    cout << "   -r, --input-ref-file     Input reference file (FASTA/FASTQ possibly gzipped and GFA)" << endl;
    cout << "                            Multiple files can be provided as a list in a TXT file (one file per line)" << endl;
    cout << "                            All k-mers of the input reference files are used" << endl;
    cout << "   -o, --output-file        Prefix for output file(s)" << endl << endl;

    cout << "   > Optional with required argument:" << endl << endl;

    cout << "   -f, --input-color-file   Input color file associated with the input graph file to update" << endl;
    cout << "   -t, --threads            Number of threads (default is 1)" << endl;
    cout << "   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)" << endl;
    cout << "   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or automatically chosen)" << endl << endl;

    cout << "   > Optional with no argument:" << endl << endl;

    cout << "   -i, --clip-tips          Clip tips shorter than k k-mers in length" << endl;
    cout << "   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length" << endl;
    cout << "   -v, --verbose            Print information messages during execution" << endl << endl;

    cout << "[PARAMETERS]: query" << endl << endl;

    cout << "  > Mandatory with required argument:" << endl << endl;

    cout << "   -g, --input-graph-file   Input graph file to query (GFA format)" << endl;
    cout << "   -q, --input-query-file   Input query file (FASTA/FASTQ possibly gzipped)" << endl;
    cout << "                            Multiple files can be provided as a list in a TXT file (one file per line)" << endl;
    cout << "   -o, --output-file        Prefix for output file" << endl;
    cout << "   -e, --ratio-kmers        Ratio of k-mers from queries that must occur in the graph (default is 0.8)" << endl << endl;

    cout << "   > Optional with required argument:" << endl << endl;

    cout << "   -f, --input-color-file   Input color file associated with the input graph file to query" << endl;
    cout << "                            Presence/absence of queries will be output for each color" << endl;
    cout << "   -t, --threads            Number of threads (default is 1)" << endl;
    cout << "   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)" << endl;
    cout << "   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or or automatically chosen)" << endl << endl;

    cout << "   > Optional with no argument:" << endl << endl;

    cout << "   -n, --inexact            Graph is searched with exact and inexact k-mers (1 substitution or indel) from queries" << endl;
    cout << "   -v, --verbose            Print information messages during execution" << endl << endl;
}

int parse_ProgramOptions(int argc, char **argv, CCDBG_Build_opt& opt) {

    int option_index = 0, c;

    const char* opt_string = "s:r:q:g:f:o:t:k:m:e:b:B:l:w:u:nidvcya";

    static struct option long_options[] = {

        {"input-seq-file",      required_argument,  0, 's'},
        {"input-ref-file",      required_argument,  0, 'r'},
        {"input-query-file",    required_argument,  0, 'q'},
        {"input-graph-file",    required_argument,  0, 'g'},
        {"input-color-file",    required_argument,  0, 'f'},
        {"output-file",         required_argument,  0, 'o'},
        {"threads",             required_argument,  0, 't'},
        {"kmer-length",         required_argument,  0, 'k'},
        {"min-length",          required_argument,  0, 'm'},
        {"ratio-kmers",         required_argument,  0, 'e'},
        {"bloom-bits",          required_argument,  0, 'b'},
        {"bloom-bits2",         required_argument,  0, 'B'},
        {"load-mbbf",           required_argument,  0, 'l'},
        {"write-mbbf",          required_argument,  0, 'w'},
        {"chunk-size",          required_argument,  0, 'u'},
        {"inexact_search",      no_argument,        0, 'n'},
        {"clip-tips",           no_argument,        0, 'i'},
        {"del-isolated",        no_argument,        0, 'd'},
        {"verbose",             no_argument,        0, 'v'},
        {"colors",              no_argument,        0, 'c'},
        {"keep-mercy",          no_argument,        0, 'y'},
        {"fasta",               no_argument,        0, 'a'},
        {"compress-tigs",       required_argument,  0, 'C'},
        {0,                     0,                  0,  0 }
    };

    if (strcmp(argv[1], "--version") == 0) return 1; // Print version
    else if (strcmp(argv[1], "--help") == 0) return 2; // print help

    if (strcmp(argv[1], "build") == 0) opt.build = true;
    else if (strcmp(argv[1], "update") == 0) opt.update = true;
    else if (strcmp(argv[1], "query") == 0) opt.query = true;
    else if (strcmp(argv[1], "helsitests") == 0) opt.helsitests = true;

    if (opt.build || opt.update || opt.query){

        while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

            switch (c) {

                case 's':
                    opt.filename_seq_in.push_back(optarg);
                    break;
                case 'r':
                    opt.filename_ref_in.push_back(optarg);
                    break;
                case 'q':
                    opt.filename_query_in.push_back(optarg);
                    break;
                case 'g':
                    opt.filename_graph_in = optarg;
                    break;
                case 'f':
                    opt.filename_colors_in = optarg;
                    break;
                case 'o':
                    opt.prefixFilenameOut = optarg;
                    break;
                case 't':
                    opt.nb_threads = atoi(optarg);
                    break;
                case 'k':
                    opt.k = atoi(optarg);
                    break;
                case 'm':
                    opt.g = atoi(optarg);
                    break;
                case 'e':
                    opt.ratio_kmers = atof(optarg);
                    break;
                case 'b':
                    opt.nb_bits_unique_kmers_bf = atoi(optarg);
                    break;
                case 'B':
                    opt.nb_bits_non_unique_kmers_bf = atoi(optarg);
                    break;
                case 'w':
                    opt.outFilenameBBF = optarg;
                    break;
                case 'l':
                    opt.inFilenameBBF = optarg;
                    break;
                case 'u':
                    opt.read_chunksize = atoi(optarg);
                    break;
                case 'n':
                    opt.inexact_search = true;
                    break;
                case 'i':
                    opt.clipTips = true;
                    break;
                case 'd':
                    opt.deleteIsolated = true;
                    break;
                case 'v':
                    opt.verbose = true;
                    break;
                case 'c':
                    opt.outputColors = true;
                    break;
                case 'y':
                    opt.useMercyKmers = true;
                    break;
                case 'a':
                    opt.outputGFA = false;
                    break;
                case 'C':
                    if (strcmp(optarg, "unitigs") == 0) {
                        opt.compress_tigs = Tigs::UNITIGS;
                    } else if (strcmp(optarg, "pathtigs") == 0) {
                        opt.compress_tigs = Tigs::PATHTIGS;
                    } else if (strcmp(optarg, "helsitigs") == 0) {
                        opt.compress_tigs = Tigs::HELSITIGS;
                    } else {
                        cerr << "-C or --compress-tigs was given, but argument is unknown: " << optarg << endl;
                        return 2;
                    }
                    break;
                default: break;
            }
        }
    }

    return 0;
}

bool check_ProgramOptions(CCDBG_Build_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [&](vector<string>& v_files) {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "Error: File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const int format = FileParser::getFileFormat(file.c_str());

                if (format >= 0) files_tmp.push_back(file); // File is FASTA/FASTQ/GFA
                else {

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != NULL){

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        size_t i = 0;

                        while (i_file_txt.getline(buffer, 4096).good()){

                            fp = fopen(buffer, "r");

                            if (fp == NULL) {

                                cerr << "Error: Could not open file at line " << i << " in file " << file << " for reading." << endl;
                                ret = false;
                                break;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }

                            ++i;
                        }

                        if (i_file_txt.fail() && (i == 0)) {

                            cerr << "Error: File " << file << " is neither FASTA, FASTQ nor GFA." << endl;
                            cerr << "If it is a list of files, it is either empty or has a line with >4096 characters." << endl;
                            ret = false;
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "Error: Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
            }
        }

        v_files = move(files_tmp);

        delete[] buffer;
    };

    // Check general parameters

    if (!opt.build && !opt.update && !opt.query && !opt.helsitests){

        cerr << "Error: No command selected (can be 'build' or 'update' or 'query' or 'helsitests')." << endl;
        ret = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if (opt.k <= 2){

        cerr << "Error: Length k of k-mers cannot be less than 3." << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed " << (MAX_KMER_SIZE - 1) << "." << endl;
        cerr << "To enable a larger k, recompile Bifrost with the appropriate MAX_KMER_SIZE variable." << endl;
        ret = false;
    }

    if (opt.g == 0){

        cerr << "Error: Length m of minimizers cannot be equal to 0." << endl;
        ret = false;
    }

    if ((opt.g >= 0) && (opt.g > opt.k - 2)) {

        cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
        ret = false;
    }

    if (opt.query){  // Check param. command build

        if (opt.prefixFilenameOut.length() == 0) {

            cerr << "Error: No output filename prefix given." << endl;
            ret = false;
        }
        else {

            const string out = opt.prefixFilenameOut + ".tsv";

            FILE* fp = fopen(out.c_str(), "w");

            if (fp == NULL) {

                cerr << "Error: Could not open file for writing output of query in TSV format: " << out << "." << endl;
                ret = false;
            }
            else {

                fclose(fp);
                if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
            }
        }

        if (opt.filename_query_in.size() == 0) {

            cerr << "Error: Missing input query files." << endl;
            ret = false;
        }
        else check_files(opt.filename_query_in);

        if ((opt.ratio_kmers < 0.0) || (opt.ratio_kmers > 1.0)) {

            cerr << "Error: Ratio of k-mers from queries that must occur in the graph cannot be less than 0.0 or more than 1.0 (" << opt.ratio_kmers << ")." << endl;
            ret = false;
        }

        if (opt.g > opt.k - 2){

            cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
            ret = false;
        }
    }
    else if (!opt.helsitests) {

        if (opt.prefixFilenameOut.length() == 0) {

            cerr << "Error: No output filename prefix given." << endl;
            ret = false;
        }
        else {

            const string out = opt.prefixFilenameOut + (opt.outputGFA ? ".gfa" : ".fasta");

            FILE* fp = fopen(out.c_str(), "w");

            if (fp == NULL) {

                cerr << "Error: Could not open file for writing output graph in GFA format: " << out << "." << endl;
                ret = false;
            }
            else {

                fclose(fp);
                if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
            }
        }

        if ((opt.filename_seq_in.size() + opt.filename_ref_in.size()) == 0) {

            cerr << "Error: Missing input files." << endl;
            ret = false;
        }
        else {

            check_files(opt.filename_seq_in);
            check_files(opt.filename_ref_in);
        }
    }

    if (opt.build){ // Check param. command build

        if (opt.read_chunksize <= 0) {

            cerr << "Error: Chunk size of reads to share among threads cannot be less than or equal to 0." << endl;
            ret = false;
        }

        if (opt.outFilenameBBF.length() != 0){

            FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

            if (fp == NULL) {

                cerr << "Error: Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                ret = false;
            }
            else {

                fclose(fp);

                if (remove(opt.outFilenameBBF.c_str()) != 0){

                    cerr << "Error: Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
                }
            }
        }

        if (opt.inFilenameBBF.length() != 0){

            if (check_file_exists(opt.inFilenameBBF)){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "Error: Input Blocked Bloom filter " << opt.inFilenameBBF << " file does not exist." << endl;
                ret = false;
            }
        }
    }

    if (opt.update || opt.query){

        if (opt.filename_graph_in.length() == 0){

            cerr << "Error: No graph file was provided in input." << endl;
            ret = false;
        }
        else if (!check_file_exists(opt.filename_graph_in)){

            cerr << "Error: The graph file does not exist." << endl;
            ret = false;
        }
        else {

            FILE* fp = fopen(opt.filename_graph_in.c_str(), "r");

            if (fp == NULL) {

                cerr << "Error: Could not read input graph file " << opt.filename_graph_in << "." << endl;
                ret = false;
            }
            else fclose(fp);
        }

        if (opt.filename_colors_in.length() != 0){

            if (!check_file_exists(opt.filename_colors_in)){

                cerr << "Error: The input color file does not exist." << endl;
                ret = false;
            }
            else {

                FILE* fp = fopen(opt.filename_colors_in.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input color file " << opt.filename_colors_in << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
        }
    }

    return ret;
}

std::basic_ostream<char>& operator<<(std::basic_ostream<char>& stream, const std::vector<std::string>& strings) {
    stream << "[";
    bool first = true;
    for (const auto& string : strings) {
        if (first) {
            first = false;
        } else {
            stream << ", ";
        }
        stream << string;
    }
    stream << "]";
    return stream;
}

int main(int argc, char **argv){

    if (argc < 2) PrintUsage();
    else {

        CCDBG_Build_opt opt;

        opt.outputColors = false; // We dont know yet if we want colors or not

        const int print = parse_ProgramOptions(argc, argv, opt); // Parse input parameters

        if (print == 1) PrintVersion();
        else if (print == 2) PrintUsage();
        else if (check_ProgramOptions(opt)) {

            if (opt.build){ // Build the graph

                if (opt.outputColors) {
                    auto start_build = std::chrono::high_resolution_clock::now();

                    ColoredCDBG<> ccdbg(opt.k, opt.g);

                    ccdbg.buildGraph(opt);
                    ccdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                    ccdbg.buildColors(opt);

                    auto stop_build = std::chrono::high_resolution_clock::now();
                    double duration_build = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_build - start_build).count()) / 1e6;
                    cout << "Took " << std::fixed << std::setprecision(3) << duration_build << "s for building dbg (input_filename = " << opt.filename_graph_in << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;

                    bool ok = true;
                    if (opt.compress_tigs != Tigs::NONE) {
                        auto start_tig = std::chrono::high_resolution_clock::now();

                        ColoredCDBG<> tigged_ccdbg(opt.k, opt.g);
                        //ccdbg.checkColors(ccdbg.getColorNames());
                        ok = tigged_ccdbg.convert_tigs(&ccdbg, opt.compress_tigs, opt.nb_threads,
                                                                 opt.prefixFilenameOut);
                        //ccdbg.checkColors(ccdbg.getColorNames());
                        //tigged_ccdbg.checkColors(ccdbg.getColorNames());
                        ccdbg = move(tigged_ccdbg);

                        auto stop_tig = std::chrono::high_resolution_clock::now();
                        double duration_tig = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_tig - start_tig).count()) / 1e6;
                        cout << "Took " << std::fixed << std::setprecision(3) << duration_tig << "s for compressing tigs (input_filename = " << opt.filename_graph_in << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;
                    }

                    if (ok) {
                        ccdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
                    }
                }
                else {
                    auto start_build = std::chrono::high_resolution_clock::now();
                    CompactedDBG<> cdbg(opt.k, opt.g);

                    cdbg.build(opt);
                    cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);

                    auto stop_build = std::chrono::high_resolution_clock::now();
                    double duration_build = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_build - start_build).count()) / 1e6;
                    cout << "Took " << std::fixed << std::setprecision(3) << duration_build << "s for building dbg (input_filename = " << opt.filename_graph_in << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;

                    bool ok = true;
                    if (opt.compress_tigs != Tigs::NONE) {
                        auto start_tig = std::chrono::high_resolution_clock::now();

                        CompactedDBG<> tigged_cdbg(opt.k, opt.g);
                        ok = tigged_cdbg.convert_tigs(&cdbg, opt.compress_tigs, opt.nb_threads,
                                                      opt.prefixFilenameOut);
                        cdbg = move(tigged_cdbg);

                        auto stop_tig = std::chrono::high_resolution_clock::now();
                        double duration_tig = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_tig - start_tig).count()) / 1e6;
                        cout << "Took " << std::fixed << std::setprecision(3) << duration_tig << "s for compressing tigs (input_filename = " << opt.filename_graph_in << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;
                    }

                    if (ok) {
                        cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.outputGFA, opt.verbose);
                    }
                }
            }
            else if (opt.update){

                if (opt.filename_colors_in.size() != 0){ // If colors in or out

                    ColoredCDBG<> ccdbg1(opt.k, opt.g);
                    ColoredCDBG<> ccdbg2(opt.k, opt.g);

                    ccdbg1.read(opt.filename_graph_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);
                    ccdbg2.buildGraph(opt);
                    ccdbg2.buildColors(opt);

                    const size_t ccdbg1_len = ccdbg1.length();
                    const size_t ccdbg2_len = ccdbg2.length();

                    ColoredCDBG<>& ccdbg_a = (ccdbg1_len > ccdbg2_len) ? ccdbg1 : ccdbg2;
                    ColoredCDBG<>& ccdbg_b = (ccdbg1_len > ccdbg2_len) ? ccdbg2 : ccdbg1;

                    ccdbg_a.merge(move(ccdbg_b), opt.nb_threads, opt.verbose);

                    ccdbg_a.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                    ccdbg_a.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
                }
                else {

                    CompactedDBG<> cdbg1(opt.k, opt.g);
                    CompactedDBG<> cdbg2(opt.k, opt.g);

                    cdbg1.read(opt.filename_graph_in, opt.nb_threads, opt.verbose);
                    cdbg2.build(opt);

                    const size_t cdbg1_len = cdbg1.length();
                    const size_t cdbg2_len = cdbg2.length();

                    CompactedDBG<>& cdbg_a = (cdbg1_len > cdbg2_len) ? cdbg1 : cdbg2;
                    CompactedDBG<>& cdbg_b = (cdbg1_len > cdbg2_len) ? cdbg2 : cdbg1;

                    cdbg_a.merge(cdbg_b, opt.nb_threads, opt.verbose);
                    cdbg_b.clear();

                    cdbg_a.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                    cdbg_a.write(opt.prefixFilenameOut, opt.nb_threads, opt.outputGFA, opt.verbose);
                }
            }
            else if (opt.query){

                if (opt.filename_colors_in.size() != 0){
                    auto start_read = std::chrono::high_resolution_clock::now();
                    ColoredCDBG<> ccdbg(opt.k, opt.g);
                    ccdbg.read(opt.filename_graph_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);
                    auto stop_read = std::chrono::high_resolution_clock::now();
                    double duration_read = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_read - start_read).count()) / 1e6;
                    cout << "Took " << std::fixed << std::setprecision(3) << duration_read << "s for reading dbg (input_filename = " << opt.filename_graph_in << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;

                    auto start_search = std::chrono::high_resolution_clock::now();
                    ccdbg.search(opt.filename_query_in, opt.prefixFilenameOut, opt.ratio_kmers, opt.inexact_search, opt.nb_threads, opt.verbose);
                    auto stop_search = std::chrono::high_resolution_clock::now();
                    double duration_search = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_search - start_search).count()) / 1e6;
                    cout << "Took " << std::fixed << std::setprecision(3) << duration_search << "s for querying dbg (query_filenames = " << opt.filename_query_in << ", out_filename_prefix = " << opt.prefixFilenameOut << ", ratio_kmers = " << opt.ratio_kmers << ", inexact_search = " << opt.inexact_search << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;
                }
                else {
                    auto start_read = std::chrono::high_resolution_clock::now();
                    CompactedDBG<> cdbg(opt.k, opt.g);
                    cdbg.read(opt.filename_graph_in, opt.nb_threads, opt.verbose);
                    auto stop_read = std::chrono::high_resolution_clock::now();
                    double duration_read = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_read - start_read).count()) / 1e6;
                    cout << "Took " << std::fixed << std::setprecision(3) << duration_read << "s for reading dbg (input_filename = " << opt.filename_graph_in << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;

                    auto start_search = std::chrono::high_resolution_clock::now();
                    cdbg.search(opt.filename_query_in, opt.prefixFilenameOut, opt.ratio_kmers, opt.inexact_search, opt.nb_threads, opt.verbose);
                    auto stop_search = std::chrono::high_resolution_clock::now();
                    double duration_search = ((double) std::chrono::duration_cast<std::chrono::microseconds>(stop_search - start_search).count()) / 1e6;
                    cout << "Took " << std::fixed << std::setprecision(3) << duration_search << "s for querying dbg (query_filenames = " << opt.filename_query_in << ", out_filename_prefix = " << opt.prefixFilenameOut << ", ratio_kmers = " << opt.ratio_kmers << ", inexact_search = " << opt.inexact_search << ", nb_threads = " << opt.nb_threads << ", verbose = " << opt.verbose << ")" << std::endl;
                }
            } else if (opt.helsitests) {
                cout << "Executing helsitests" << endl;

                ColoredCDBG<> ccdbg(4, 2);
                ccdbg.initialiseHashtable(20);
                ccdbg.addUnitig("ACACACA", 0);
                ccdbg.addUnitig("ACACACT", 1);

                for (const auto& unitig : ccdbg) {
                    const string seq(unitig.referenceUnitigToString());
                    cout << seq << endl;
                }

                {
                    const auto unitig = ccdbg.findUnitig("ACAC", 0, 4);
                    const string seq(unitig.referenceUnitigToString());
                    cout << "Kmer ACAC is in " << seq << endl;
                }

                {
                    const auto unitigs = ccdbg.searchSequence("ACACAG", false, false, false, false, false);
                    cout << "Search results of ACACAG" << endl;
                    for (const auto& entry : unitigs) {
                        const auto& id = entry.first;
                        const string seq(entry.second.referenceUnitigToString());
                        cout << "Hit " << id << " in unitig " << seq << endl;
                    }
                }

                {
                    const auto unitigs = ccdbg.searchSequence("ACACA", false, false, false, false, false);
                    cout << "Search results of ACACAG" << endl;
                    for (const auto& entry : unitigs) {
                        const auto& id = entry.first;
                        const string seq(entry.second.referenceUnitigToString());
                        cout << "Hit " << id << " in unitig " << seq << endl;
                    }
                }

                {
                    const auto unitigs = ccdbg.findAll(Kmer("ACAC"), false);cout << "Search results of ACACAG" << endl;
                    cout << "findAll results of ACAC" << endl;
                    for (const auto& entry : unitigs) {
                        const string seq(entry.referenceUnitigToString());
                        cout << "Hit in unitig " << seq << endl;
                    }

                    ccdbg.buildColors(opt);
                }

                cout << "Finished helsitests" << endl;
            }
        }
    }
}
