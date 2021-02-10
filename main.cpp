// C/C++/C++11/C++17 headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include <algorithm>
#include <set>
#include <cinttypes>
#include <cstdlib>
#include <iterator>
#include <functional>
#include <cassert>
#include <experimental/filesystem>

// bifrost header
#include <bifrost/ColoredCDBG.hpp>

// robinhood hashing header
#include "robin_hood.h"

int main(int argc, char *argv[]) {
    const std::string graph_infile = argv[1];
    const std::string kmer_count = argv[2];
    const std::string output_prefix = argv[3];

    bool b;
    istringstream(argv[4]) >> std::boolalpha >> b;

    bool build = b;
    int num_threads = 4;

    int kmer = 31;

    bool verb = true;

    // build graph
    CompactedDBG<> cdbg;

    if (build) {
        CCDBG_Build_opt opt;
        opt.k = kmer;
        opt.nb_threads = num_threads;
        opt.verbose = verb;
        opt.prefixFilenameOut = output_prefix;

        std::string filename;
        std::ifstream in_graph(graph_infile);
        while (std::getline(in_graph, filename)) {
            opt.filename_ref_in.push_back(filename);
        }

        cdbg.build(opt);
        cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
        cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
    } else {
        cdbg.read(graph_infile, num_threads);
        // set local variables
        kmer = cdbg.getK();
    }

    // parse fasta
    robin_hood::unordered_map<std::string, int> kmer_count_map;
    std::string line, id, DNA_sequence;
    std::ifstream in_count(kmer_count);
    while (std::getline(in_count, line)) {

        // line may be empty so you *must* ignore blank lines
        // or you have a crash waiting to happen with line[0]
        if (line.empty())
            continue;

        if (line[0] == '>') {
            // output previous line before overwriting id
            id = line.substr(1);
            DNA_sequence.clear();
        } else {//  if (line[0] != '>'){ // not needed because implicit
            DNA_sequence += line;
        }

        // output final entry
        // but ONLY if id actually contains something
        if (!id.empty())
            kmer_count_map[DNA_sequence] = std::stoi(id);
    }

    // create map to store each unitigs average coverage
    robin_hood::unordered_map<std::string, size_t> unitig_coverage;
    // get all head kmers for parrellelisation
    for (const auto um : cdbg)
    {
        size_t num_kmers = um.size - kmer + 1;
        int total_coverage = 0;

        std::string unitig_seq = um.referenceUnitigToString();

        for (size_t i = 0; i < num_kmers; i++)
        {
            size_t kmer_cov = 0;
            std::string kmer_seq = um.getUnitigKmer(i).toString();
            if (kmer_count_map.find(kmer_seq) != kmer_count_map.end())
            {
                kmer_cov = kmer_count_map[kmer_seq];
            } else
            {
                kmer_seq = reverse_complement(kmer_seq);
                kmer_cov = kmer_count_map[kmer_seq];
            }

            total_coverage += kmer_cov;
        }

        double average_coverage = total_coverage / num_kmers;
        unitig_coverage[unitig_seq] = average_coverage;
    }

    // write to file
    ofstream outfile;
    outfile.open(output_prefix + "_unitig_coverage.txt");

    for (const auto& unitig : unitig_coverage)
    {
        outfile << unitig.first << "\t" << unitig.second << "\n";
    }
    outfile.close();

    return 0;
}
