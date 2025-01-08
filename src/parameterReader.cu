#include "parameterReader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cuda_runtime.h>

#define NUM_COLUMNS 35

// Define the kernel function to process the parameters with the data type given first (unlike void)
std::vector< std::vector<float> > parameterReader(const char * name)
{
    using namespace std;  // Use the standard namespace
    using namespace boost;  // Use the Boost namespace
    vector< vector<float> > params_out;  // Vector to store the output pairs of floats

    ifstream in(name);  // Open the input file stream with the given file name
    if (!in.is_open()) return params_out;  // If the file cannot be opened, return an empty vector

    vector<string> lineVector;  // Vector to store the tokens of a line
    string line;  // String to store each line of the file

    // Read the file line by line
    while (getline(in, line))
    {
        // Tokenize the line using Boost tokenizer with escaped list separator
        tokenizer< escaped_list_separator<char> > tok(line);
        lineVector.assign(tok.begin(), tok.end());  // Assign the tokens to the lineVector

        if (lineVector.size() != NUM_COLUMNS) {
            cout<< "Error: The number of columns in the file is not equal to 35" << endl;
            cout<< "Number of columns in the file: " << lineVector.size() << endl;
            break;  // If the line does not have exactly 35 tokens, break the loop
        }

        // Check if the line is a list of header variable names, or numeric 
        if (lineVector[0] != "protocol"){
            // Convert the tokens to floats and add them as a pair to the output vector
            std::vector<float> row(NUM_COLUMNS);
            for (int i = 0; i < NUM_COLUMNS; ++i) {
                row[i] = atof(lineVector[i].c_str());
            }
            params_out.push_back(row);
        }
    }
    
    // Print the read data
    cout << "Read parameters: " << endl;
    for(int i = 0; i < params_out.size(); ++i){
        for (int j = 0; j < NUM_COLUMNS; ++j){
            cout << params_out[i][j] << ", ";
        }
        cout << endl;
    }
    cout << endl;
    return params_out;  // Return the output vector
}

