#include <sharg/all.hpp> // includes all necessary headers
#include <filesystem>
 
#include <iostream>
#include <string>
#include <vector>
#include <cmath>


bool is_base_GC(char c)
{
    switch (c)
    {
        case 'G':
        case 'C':
        case 'g':
        case 'c':
            return true;
        default:
            return false;
    }
}

size_t GC_of_interval(std::string & sequence, size_t it_start, size_t w)
{
    size_t result(0);
    for (size_t i = it_start; i != it_start+w; ++i)
    {
        result = result << 1;
        result += is_base_GC(sequence[i]);
        result &= ((1 << (w))-1);
    }
    return result;
}

void run_program(size_t wsize)
{
    std::string line;
    std::string buffer;
    
    size_t GCs(0);
    float result(0.0f);
    while (std::getline(std::cin, line))
    {
        if (line[0] != '>')
        {
            // append the line to the buffer
            buffer += line;
            // if the buffer is long enough, compute the hash values
            if (buffer.size() >= wsize)
            {
                // compute GC content of first wsize letters
                GCs = GC_of_interval(buffer,0,wsize);
                result = float(std::popcount(GCs)) / float(wsize);
                // padding of first wsize/2 result values
                for (size_t i = 0; i < size_t(wsize/2)+1; i++)
                {
                    // cout result
                    std::cout << result << '\n';
                }
                // remove the first w letters from the buffer
                buffer.erase(0,wsize);
                break;
            }
        }
    }
    // main loop
    // while loop through the rest of standard input
    while (std::getline(std::cin, line))
    {
        if (line[0] != '>')
        {
            // if the buffer is not empty, then add the buffer to the front of the line
            if (buffer.size() > 0) { line = buffer + line;}
            // iterate over the line
            for (char c : line)
            {
                // shift the GCs value to the left
                GCs = GCs << 1;
                // add the new base to the GCs value
                GCs += is_base_GC(c);
                // mask the GCs value to the correct size
                GCs &= ((1 << (wsize))-1);
                // compute the result
                result = float(std::popcount(GCs)) / float(wsize);
                // cout result
                std::cout << result << '\n';
            }
        }
    }
    // padding of last wsize/2 result values
    for (size_t i = 0; i < size_t(wsize/2); i++)
    {
        // cout result
        std::cout << result << '\n';
    }
}


struct cmd_arguments
{
    size_t wsize{};
};
 
void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "VM";
    parser.info.short_description = "computes the GC content within a window w.";
    parser.info.version = "0.0.1";

    parser.add_option(args.wsize, sharg::config{
        .short_id = 'w',
        .long_id = "wsize",
        .description = "window size w must be 2 < w < 21."});
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
    run_program(args.wsize);
 
    return 0;
}