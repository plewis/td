//
//  main.cpp
//  TreeDist project
//
//  Created by Paul O. Lewis on 2024-01-31.
//

#include <iostream>
#include <memory>
#include <regex>
#include "output_manager.hpp"

using namespace std;
using namespace treedist;

OutputManager om;

#include "node.hpp"
#include "lot.hpp"
#include "treedist.hpp"
#include "xtreedist.hpp"

Lot rng;

// static data member initializations
string     TreeDist::_program_name         = "treedist";
unsigned        TreeDist::_major_version        = 1;
unsigned        TreeDist::_minor_version        = 0;

const double    Node::_smallest_edge_length     = 1.0e-12;

int main(int argc, const char * argv[]) {

    TreeDist td;
    try {
        td.processCommandLineOptions(argc, argv);
        td.run();
    }
    catch(std::exception & x) {
        cerr << "Exception: " << x.what() << endl;
        cerr << "Aborted." << endl;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

    return 0;
}
