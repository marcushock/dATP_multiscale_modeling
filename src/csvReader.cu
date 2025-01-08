#include "csvReader.h"  // Include the header file for the csvReader function

#include <iostream>     // Include the iostream library for input and output operations
#include <fstream>      // Include the fstream library for file stream operations
#include <string>       // Include the string library for string operations
#include <boost/tokenizer.hpp>  // Include the Boost tokenizer library for parsing CSV lines

// Function to read a CSV file and return a vector of pairs of floats
std::vector< std::pair<float, float> > csvReader(const char * name)
{
    using namespace std;  // Use the standard namespace
    using namespace boost;  // Use the Boost namespace
    vector< pair<float, float> > out;  // Vector to store the output pairs of floats

    ifstream in(name);  // Open the input file stream with the given file name
    if (!in.is_open()) return out;  // If the file cannot be opened, return an empty vector

    vector<string> lineVector;  // Vector to store the tokens of a line
    string line;  // String to store each line of the file

    // Read the file line by line
    while (getline(in, line))
    {
        // Tokenize the line using Boost tokenizer with escaped list separator
        tokenizer< escaped_list_separator<char> > tok(line);
        lineVector.assign(tok.begin(), tok.end());  // Assign the tokens to the lineVector

        if (lineVector.size() != 2) {
            cout<<"The actual number of elements in the line vector is "<<lineVector.size()<<endl;
            break; // If the line does not have exactly 2 tokens, break the loop

        }  

        // Convert the tokens to floats and add them as a pair to the output vector
        out.push_back(make_pair(atof(lineVector[0].c_str()), atof(lineVector[1].c_str())));
    }
    
    // Print the read data
    cout << "Read data: " << endl;
    for(int i = 0; i < out.size(); ++i)
        cout << out[i].first << ", " << out[i].second << endl;
    
    return out;  // Return the output vector
}
