#ifndef PARAMETERREADER_H
#define PARAMETERREADER_H

#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <string>
#include <boost/tokenizer.hpp>

std::vector< std::vector<float> > parameterReader(const char * name);


#endif // PARAMETERREADER_H
