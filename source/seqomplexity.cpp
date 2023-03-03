#include <sharg/all.hpp> // includes all necessary headers
#include <filesystem>
 
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
 
// This is the program!
// Take a look at it if you are interested in an example of parsing a data file.
// -----------------------------------------------------------------------------
 
void run_program(
        std::filesystem::path & input,
        std::filesystem::path & output,
        size_t wsize,
        std::vector<uint8_t> kmers)
{
    seqan3::sequence_file_input fin{input};
 
    // `&&` is important because seqan3::views::chunk returns temporaries!
    for (auto & record : fin)
    {
        // `records` contains 10 elements (or less at the end)
        seqan3::debug_stream << "ID:  " << record.id() << '\n'; // prints first ID in batch
    }
}
// -----------------------------------------------------------------------------
 
struct cmd_arguments
{
    std::filesystem::path input{};
    std::filesystem::path output{};
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
    parser.add_option(args.output, sharg::config{
        .short_id = 'o',
        .long_id = "output",
        .description = "path to output memmap.",
        .required = true,
        .validator = sharg::output_file_validator{"fa", "fasta"}});
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
    run_program(args.input, args.output, args.wsize, args.kmers);
 
    return 0;
}

