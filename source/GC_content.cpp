#include <sharg/all.hpp> // includes all necessary headers
#include <filesystem>
 
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
 

size_t is_base_GT(seqan3::dna5 base)
{
    switch (seqan3::to_rank(base))
    {
        case 1:
            return 1;
        case 2:
            return 1;
        default:
            return 0;
    }
}

size_t GC_of_interval(std::vector<seqan3::dna5> & sequence, size_t it_start, size_t w)
{
    size_t result(0);
    for (size_t i = it_start; i != it_start+w; ++i)
    {
        result = result << 1;
        result += is_base_GT(sequence[i]);
        result &= ((1 << (w))-1);
    }
    return result;
}

void sequence_GC(
    std::vector<seqan3::dna5> & sequence,
    size_t W)
{
    // init: GC for window pat of sequence for position W/2
    // GC of 
    size_t N(sequence.size());
    // --- setup --- //
    //std::vector<float> results(N,float(0.0));
    size_t GCs = GC_of_interval(sequence,0,W);
    float result = float(std::popcount(GCs)) / float(W);
    // padding
    for (size_t i = 0; i < size_t(W/2)+1; i++)
    {
        std::cout << result << ' ' << GCs << ' ' << seqan3::to_char(sequence[i]) << '\n';
    }
    
    // main loop
    for (size_t i = size_t(W/2)+1; i < N-size_t(W/2); i++)
    {
        GCs = GCs << 1;
        GCs += size_t(is_base_GT(sequence[i+size_t(W/2)]));
        GCs &= (1 << W)-1;
        result = float(std::popcount(GCs)) / float(W);
        std::cout << result << ' ' << seqan3::to_char(sequence[i+size_t(W/2)]) << ' ' << GCs << ' ' << float(std::popcount(GCs)) << ' ' << float(W) << '\n';
    }
    // padding
    for (size_t i = N-size_t(W/2); i < N; i++)
    {
        std::cout << result << ' ' << GCs << ' ' << seqan3::to_char(sequence[i]) << '\n';
    }
    return ;

}

void run_program(
        std::filesystem::path & input,
        size_t wsize)
{
    seqan3::sequence_file_input fin{input};
    for (auto & record : fin)
    {
        //seqan3::debug_stream << "ID:  " << record.id() << '\n'; // prints first ID in batch
        sequence_GC(record.sequence(), wsize);
    }
    return ;
}
// -----------------------------------------------------------------------------
 
struct cmd_arguments
{
    std::filesystem::path input{};
    size_t wsize{};
};
 
void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "VM";
    parser.info.short_description = "computes the sequence GC content for a given sliding window wsize.";
    parser.info.version = "0.0.1";

    // ... add information, options, flags and positional options
    parser.add_option(args.input, sharg::config{
        .short_id = 'i',
        .long_id = "input",
        .description = "reference file in fasta format.",
        .required = true,
        .validator = sharg::input_file_validator{{"fa", "fasta"}}});
    parser.add_option(args.wsize, sharg::config{
        .short_id = 'w',
        .long_id = "wsize",
        .description = "window size w must be 1 < w < 64."});
}
 
int main(int argc, char ** argv)
{
    sharg::parser parser{"sequence-GC-content", argc, argv}; // initialise myparser
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
    run_program(args.input, args.wsize);
 
    return 0;
}

