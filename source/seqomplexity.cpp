#include <sharg/all.hpp> // includes all necessary headers
#include <filesystem>
 
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
 

size_t kmer_to_hash(std::vector<seqan3::dna5> & sequence, size_t it_start, size_t it_end, size_t base=3)
{
    size_t hashvalue = 0;
    for (size_t i = it_start; i < it_end; i++)
    {
        hashvalue = (hashvalue << base) + seqan3::to_rank(sequence[i]);
    }
    return hashvalue;
}

std::vector<size_t> sequence_to_kmer_hashes(
    std::vector<seqan3::dna5> & sequence,
    size_t k,
    size_t it_start,
    size_t it_end)
{
    //body
    // return [kmer_to_hash(seq[i:i+k]) for i in range(len(seq)-k+1)]
    std::vector<size_t> hashvalues(it_end - it_start -k+1);
    for (size_t i = it_start; i < it_end-k+1; i++)
    {
        hashvalues[i-it_start] = kmer_to_hash(sequence,i,i+k);
    }
    return hashvalues;
}

size_t extend_hash(size_t kmer_hash, uint8_t letter_dna5 , size_t k, size_t base=3)
{
    return ((kmer_hash << base) + letter_dna5) & ((1 << (k*base))-1);
}

float product(std::vector<size_t> & a, std::vector<size_t> & b)
{
    float result(1.0);
    for (size_t i = 0; i < std::min(a.size(),b.size()); i++)
    {
        result *= (float(a[i]) / float(b[i]));
    }
    return result;
}

void sequence_complexity(
    std::vector<seqan3::dna5> & sequence,
    size_t W,
    std::vector<uint8_t> kmers)
{
    size_t N(sequence.size());
    // --- setup --- //
    //std::vector<float> results(N,float(0.0));
    std::vector<size_t> MAX_UNIQUE_HASHES(kmers.size());
    for (size_t i = 0; i < kmers.size(); i++){
        MAX_UNIQUE_HASHES[i] = std::min(W-kmers[i]+1,size_t(pow(4,kmers[i])));
    }
    // --- init --- //
    // starting k-mers
    std::vector<std::vector<size_t>> kmer_hashes(kmers.size());
    for (size_t i = 0; i < kmers.size(); i++)
    {
        kmer_hashes[i] = sequence_to_kmer_hashes(sequence, kmers[i], 0, W);
    }
    // starting unique k-mer counts
    std::vector<size_t> current_unique_n_hashes(kmers.size());
    for (size_t i = 0; i < kmers.size(); i++)
    {
        current_unique_n_hashes[i] = std::set<size_t>(kmer_hashes[i].begin(), kmer_hashes[i].end() ).size();
    }
    // save resulting scores
    //results[size_t(W/2)] = product(current_unique_n_hashes,MAX_UNIQUE_HASHES);

    // cout padding
    for (size_t i = 0; i < size_t(W/2); i++)
    {
        std::cout << product(current_unique_n_hashes,MAX_UNIQUE_HASHES) << '\n';
    }
    // main loop
    for (size_t i = size_t(W/2)+1; i < N-size_t(W/2); i++)
    {
        size_t h = seqan3::to_rank(sequence[i+size_t(W/2)]);
        for (size_t j = 0; j < kmers.size(); j++)
        {
            size_t k = kmers[j];
            size_t position = size_t((i-size_t(W/2)-2)%(W-k+1));
            size_t last_position=size_t((i-size_t(W/2)-2)%(W-k+1));
            size_t old_hash = kmer_hashes[j][last_position];
            bool last_hash_unique = std::count(kmer_hashes[j].begin(), kmer_hashes[j].end(), kmer_hashes[j][position]) == 1;
            current_unique_n_hashes[j] -= size_t(last_hash_unique);
            size_t new_hash = extend_hash(old_hash,h,k);
            kmer_hashes[j][position]=new_hash;
            // update count of unique elements in k-mers
            bool new_hash_unique = std::count(kmer_hashes[j].begin(), kmer_hashes[j].end(), new_hash) == 1;
            current_unique_n_hashes[j] += int(new_hash_unique);
        }
        std::cout << product(current_unique_n_hashes,MAX_UNIQUE_HASHES) << '\n';
        //results[i] = product(current_unique_n_hashes,MAX_UNIQUE_HASHES);
    }
    for (size_t i = N-size_t(W/2)+1; i < N; i++)
    {
        std::cout << product(current_unique_n_hashes,MAX_UNIQUE_HASHES) << '\n';
    }
    
    return ;

}

void run_program(
        std::filesystem::path & input,
        //std::filesystem::path & output,
        size_t wsize,
        std::vector<uint8_t> kmers)
{
    seqan3::sequence_file_input fin{input};
    for (auto & record : fin)
    {
        //seqan3::debug_stream << "ID:  " << record.id() << '\n'; // prints first ID in batch
        sequence_complexity(record.sequence(), 15, kmers);
    }
}
// -----------------------------------------------------------------------------
 
struct cmd_arguments
{
    std::filesystem::path input{};
    //std::filesystem::path output{};
    size_t wsize{};
    std::vector<uint8_t> kmers{};
};
 
void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "VM";
    parser.info.short_description = "computes the sequence complexity for different k within a window w.";
    parser.info.version = "0.0.1";

    // ... add information, options, flags and positional options
    parser.add_option(args.input, sharg::config{
        .short_id = 'i',
        .long_id = "input",
        .description = "reference file in fasta format.",
        .required = true,
        .validator = sharg::input_file_validator{{"fa", "fasta"}}});
    // parser.add_option(args.output, sharg::config{
    //     .short_id = 'o',
    //     .long_id = "output",
    //     .description = "path to output memmap.",
    //     .required = true,
    //     .validator = sharg::output_file_validator{"fa", "fasta"}});
    parser.add_option(args.wsize, sharg::config{
        .short_id = 'w',
        .long_id = "wsize",
        .description = "window size w must be 2 < w < 21."});
    parser.add_option(args.kmers, sharg::config{
        .short_id = 'k',
        .long_id = "kmers",
        .description = "any k must be 1 < k < w+1."});
}
 
int main(int argc, char ** argv)
{
    sharg::parser parser{"sequence-complexity", argc, argv}; // initialise myparser
    cmd_arguments args{};
    initialise_parser(parser, args);
 
    try
    {
        parser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << "[Winter has come] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
 
    // parsing was successful !
    // we can start running our program
    run_program(args.input, args.wsize, args.kmers);
 
    return 0;
}

