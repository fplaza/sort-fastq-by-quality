#include "Parameters.hh"
#include <iostream>
#include <cstdlib>
#include <parallel/algorithm>
#include <boost/unordered_map.hpp>
#include "FastqMultiReader.hh"
#include "Phred.hh"
#include "FastqPrinter.hh"

struct DecreasingHighQualityBasesReadsSorter
{
    DecreasingHighQualityBasesReadsSorter(const boost::unordered_map<std::string,size_t>& reads_num_high_quality_bases)
        : reads_num_high_quality_bases_(reads_num_high_quality_bases)
    {};

    bool operator()(const FastqEntry& fastq_entry1, const FastqEntry& fastq_entry2)
    {
        return reads_num_high_quality_bases_.find(fastq_entry1.id)->second > reads_num_high_quality_bases_.find(fastq_entry2.id)->second;
    }

    const boost::unordered_map<std::string,size_t>& reads_num_high_quality_bases_;
};
int main(int argc, char *argv[])
{
    try
    {
        const Parameters& parameters = Parameters::parse(argc, argv);
        std::cout << parameters << std::endl;

        FastqMultiReader fastq_multi_reader(parameters.fastq_files);
        std::vector<FastqEntry> all_fastq_entries;
        std::vector<FastqEntry> chunk_fastq_entries;

        std::cout << "Loading reads..." << std::endl;
        while (fastq_multi_reader.next_chunk(chunk_fastq_entries))
        {
            all_fastq_entries.insert(all_fastq_entries.end(), chunk_fastq_entries.begin(), chunk_fastq_entries.end());
        }
        std::cout << all_fastq_entries.size() << " reads loaded.\n" << std::endl;

        std::cout << "Counting high quality bases..." << std::endl;
        boost::unordered_map<std::string,size_t> reads_num_high_quality_bases;
        reads_num_high_quality_bases.rehash(all_fastq_entries.size());
        for (std::vector<FastqEntry>::const_iterator fastq_entry = all_fastq_entries.begin(); fastq_entry < all_fastq_entries.end(); ++fastq_entry)
        {

            reads_num_high_quality_bases[fastq_entry->id] = 0;
        }

        #pragma omp parallel for
        for (std::vector<FastqEntry>::const_iterator fastq_entry = all_fastq_entries.begin(); fastq_entry < all_fastq_entries.end(); ++fastq_entry)
        {
            size_t num_high_quality_bases = 0;
            for (size_t curr_qual = 0; curr_qual < fastq_entry->quals.size(); ++curr_qual)
            {
                const size_t curr_qual_value = Phred::decode_phred_33(fastq_entry->quals[curr_qual]);
                if (curr_qual_value >= parameters.high_qual_thr)
                    ++num_high_quality_bases;
            }
            reads_num_high_quality_bases[fastq_entry->id] = num_high_quality_bases;
        }

        std::cout << "Done.\n" << std::endl;

        
        std::cout << "Sorting reads by number of high quality bases..." << std::endl;
        __gnu_parallel::sort(all_fastq_entries.begin(), all_fastq_entries.end(), DecreasingHighQualityBasesReadsSorter(reads_num_high_quality_bases)); 
        std::cout << "Done.\n" << std::endl;

        std::cout << "Writing sorted reads..." << std::endl;
        FastqPrinter fastq_printer(parameters.output_file);
        fastq_printer.print_fastq_entries(all_fastq_entries);
        std::cout << "Done." << std::endl;

        std::exit(0);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::exit(-1);
    }
}
