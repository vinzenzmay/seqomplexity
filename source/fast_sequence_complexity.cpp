#include <sharg/all.hpp> // includes all necessary headers
#include <filesystem>
 
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

size_t hash_dna5(char c)
{
    switch (c)
    {
        case 'A':
            return 0;
        case 'a':
            return 0;
        case 'C':
            return 1;
        case 'c':
            return 1;
        case 'G':
            return 2;
        case 'g':
            return 2;
        case 'T':
            return 3;
        case 't':
            return 3;
        default:
            return 4;
    }
}

// this function recieves a start iterator and an end iterator
// and returns a hash value for the sequence between the two iterators
size_t kmer_to_hash(std::string::iterator it_start, std::string::iterator it_end, size_t base=3)
{
    size_t hashvalue = 0;
    for (auto it = it_start; it != it_end; it++)
    {
        hashvalue = (hashvalue << base) + hash_dna5(*it);
    }
    return hashvalue;
}

// this function is ran once at initialization
// it computes the maximum number of unique hashes for each kmer
// in a given sequence of characters
std::vector<size_t> sequence_to_kmer_hashes(
    std::string::iterator it_start,
    std::string::iterator it_end,
    size_t k)
{
    // loop over the sequence and calculate the hash value for each kmer
    std::vector<size_t> hashvalues{};
    for (auto it = it_start; it != it_end-k+1; it++)
    {
        hashvalues.push_back(kmer_to_hash(it,it+k));
    }
    return hashvalues;
}

size_t extend_hash(size_t kmer_hash, size_t letter_dna5 , size_t k, size_t base=3)
{
    return ((kmer_hash << base) + letter_dna5) & ((1 << (k*base))-1);
}


// this function computes the product of the element wise division of two arrays.
// the arrays must have the same size. The result is a float.
float product(size_t *a, size_t *b, size_t size)
{
    float result(1.0);
    for (size_t i = 0; i < size; i++)
    {
        result *= (float(a[i]) / float(b[i]));
    }
    return result;
}

size_t count_equal_kmer_hashes(size_t * kmer_hashes, size_t start, size_t end, size_t value)
{
    size_t count(0);
    for (size_t i = start; i < end; i++)
    {
        if (kmer_hashes[i] == value)
        {
            count++;
        }
    }
    return count;
}

size_t calc_position(size_t w, size_t k, size_t i, size_t index)
{
    size_t alpha = w-k+1;
    return index+((alpha+i-1)%alpha);
}
size_t calc_last_position(size_t w, size_t k, size_t i, size_t index)
{
    size_t alpha = w-k+1;
    return index+((alpha+i-2)%alpha);
}

int run_program(
        size_t wsize,
        std::vector<uint8_t> kmers)
{
    // check if the window size is larger than the maximum kmer size
    if (wsize < *std::max_element(kmers.begin(), kmers.end()))
    {
        std::cerr << "The window size must be larger than the maximum kmer size." << std::endl;
        exit(1);
    }
    // check if wsize is odd and 2 < wsize < 21.
    if (wsize % 2 == 0 || wsize < 2 || wsize > 21)
    {
        std::cerr << "The window size must be an odd number between 2 and 21." << std::endl;
        exit(1);
    }
    size_t nk(kmers.size());
    size_t kmers_array[nk] = {0};
    std::copy(kmers.begin(), kmers.end(), kmers_array);
    size_t MAX_UNIQUE_HASHES[nk];
    // compute the maximum number of unique hashes for each k
    for (size_t i = 0; i < nk; i++){
        MAX_UNIQUE_HASHES[i] = std::min(wsize-kmers_array[i]+1,size_t(pow(4,kmers_array[i]))+1);
    }
    std::string line;
    std::string buffer;
    // replace the vector of vectors with an array. every array according to a kmer is saved in
    // a certain interval on the array. this way, we can use the same array for all kmers.
    size_t kmer_hashes_index[nk] = {0};
    for (size_t i = 1; i < nk; i++)
    {
        kmer_hashes_index[i] = kmer_hashes_index[i-1]
            + wsize - kmers_array[i-1]+1;
    }
    size_t kmer_hashes_size = kmer_hashes_index[nk-1]+wsize-kmers_array[nk-1]+1;
    size_t kmer_hashes[kmer_hashes_size]= {0};
    // std::vector<std::vector<size_t>> kmer_hashes(kmers.size());
    size_t current_unique_n_hashes[nk] = {0};
    // first, fill the buffer and compute the initial hash values
    while (std::getline(std::cin, line))
    {
        if (line[0] != '>')
        {
            // append the line to the buffer
            buffer += line;
            // if the buffer is long enough, compute the hash values
            if (buffer.size() >= wsize)
            {
                for (size_t ki = 0; ki < nk; ki++)
                {
                    std::vector<size_t> kmer_hashes_init = sequence_to_kmer_hashes(buffer.begin(),buffer.begin()+wsize,kmers_array[ki]);
                    // i is the index of the kmer in the window - offset
                    // the offset is given in the index
                    // std::cout << "calc first hashes with k = " << kmers_array[ki] << std::endl;
                    for (size_t i = 0; i < wsize-kmers_array[ki]+1; i++)
                    {
                        kmer_hashes[kmer_hashes_index[ki] + i] = kmer_hashes_init[i];
                        // DEBUG START
                        // print kmer_hashes from 0 to kmer_hashes_index[nk-1]+kmers_array[nk-1]
                        // for (size_t i = 0; i < kmer_hashes_size; i++)
                        // {
                        //     std::cout << kmer_hashes[i] << " ";
                        // }
                        // std::cout << std::endl;
                        // DEBUG END
                    }
                }
                // starting unique k-mer counts
                for (size_t ki = 0; ki < nk; ki++)
                {
                    // std::cout << "count unique k-mers" << std::endl;
                    size_t interval_start = kmer_hashes_index[ki];
                    size_t interval_end = interval_start + wsize - kmers_array[ki]+1;
                    // copy the interval of kmer_hashes to a set and count the elements
                    current_unique_n_hashes[ki] = std::set<size_t>(
                        kmer_hashes+interval_start,
                        kmer_hashes+interval_end).size();
                }
                // padding to fill the first half of the window
                for (size_t i = 0; i < size_t(wsize/2)+1; i++)
                {
                    // std::cout << i << '\n';
                    std::cout << product(current_unique_n_hashes,MAX_UNIQUE_HASHES,nk) << '\n';
                }
                // remove the first w letters from the buffer
                buffer.erase(0,wsize);
                // print buffer with cout
                break;
            }
        }
    }
    // main loop
    // while loop through the rest of standard input

    size_t i(0);
    while (std::getline(std::cin, line) || buffer.size() > 0)
    {
        if (line[0] != '>' || buffer.size() > 0)
        {
            // if the buffer is not empty, then add the buffer to the front of the line
            if (buffer.size() > 0) { line = buffer + line; buffer.erase();}
            // iterate over the line
            for (char c : line)
            {
                // DEBUG START
                // if (c == '\n') { std::cout << "found line end.\n"; continue; }
                // DEBUG END
                i++;
                size_t h = hash_dna5(c);
                for (size_t j = 0; j < nk; j++)
                {
                    size_t k = kmers_array[j];
                    size_t k_base_length(wsize-k+1);
                    size_t position = calc_position(wsize,k,i,kmer_hashes_index[j]);
                    size_t last_position=calc_last_position(wsize,k,i,kmer_hashes_index[j]);
                    size_t old_hash = kmer_hashes[last_position];
                    bool last_hash_unique = count_equal_kmer_hashes(kmer_hashes,kmer_hashes_index[j],kmer_hashes_index[j]+wsize-k+1,kmer_hashes[position]) == 1;
                    current_unique_n_hashes[j] -= size_t(last_hash_unique);
                    size_t new_hash = extend_hash(old_hash,h,k);
                    kmer_hashes[position]=new_hash;
                    // update count of unique elements in k-mers
                    bool new_hash_unique = count_equal_kmer_hashes(kmer_hashes,kmer_hashes_index[j],kmer_hashes_index[j]+wsize-k+1,new_hash) == 1;
                    current_unique_n_hashes[j] += int(new_hash_unique);
                    // for (size_t hashval : kmer_hashes[j])
                    // {
                    //     std::cout << hashval << '\t';
                    // }
                    // std::cout << '\n';
                    // std::cout << "position: " << position << " last position: " << last_position << " kmer: " << k << " old hash: " << old_hash << " new hash: " << new_hash << " new letter: " << seqan3::to_char(sequence[i+wsize-1]) << '\n';
                }
            std::cout << product(current_unique_n_hashes,MAX_UNIQUE_HASHES,nk) << '\n';
            }
            buffer.clear();
        // std::cout << i+size_t(wsize/2) << '\n';
        // print_current_unique_hashes(current_unique_n_hashes);
        //results[i] = product(current_unique_n_hashes,MAX_UNIQUE_HASHES);  
        }
    }
    // print the last half of the window
    for (size_t i = 0; i < size_t(wsize/2); i++)
    {
        std::cout << product(current_unique_n_hashes,MAX_UNIQUE_HASHES,nk) << '\n';
    }
    return 0;
}


struct cmd_arguments
{
    size_t wsize{};
    std::vector<uint8_t> kmers{};
};
 
void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "VM";
    parser.info.short_description = "computes the sequence complexity for different k within a window w.";
    parser.info.version = "0.0.1";

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
    run_program(args.wsize, args.kmers);
 
    return 0;
}