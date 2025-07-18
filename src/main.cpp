#include <iostream>
#include "op.hpp"

using namespace strom;

// static data member initializations
string  Strom::_program_name        = "kfdist";
unsigned     Strom::_major_version       = 1;
unsigned     Strom::_minor_version       = 0;
const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[]) {
    Strom strom;
    try {
        strom.processCommandLineOptions(argc, argv);
        strom.run();
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
