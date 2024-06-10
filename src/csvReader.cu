#include "csvReader.h"

#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <string>
#include <boost/tokenizer.hpp>
std::vector< std::pair<float, float> > csvReader(const char * name)
{
    using namespace std;
    using namespace boost;
    vector< pair<float, float> > out;

    ifstream in(name);
    if (!in.is_open()) return out;

    vector<string> lineVector;
    string line;

    while (getline(in,line))
    {
        tokenizer< escaped_list_separator<char> > tok(line);
        lineVector.assign(tok.begin(),tok.end());

        if (lineVector.size() != 2) break;

        out.push_back(make_pair(atof(lineVector[0].c_str()), atof(lineVector[1].c_str())));
    }
    cout << "Read data: " << endl;
    for(int i = 0; i < out.size(); ++i)
        cout << out[i].first << ", " << out[i].second << endl;
    return out;
}
