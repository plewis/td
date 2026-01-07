// Uncomment the line below to save incompatibility graphs as dot files for debugging purposes
//#define OP_SAVE_DOT_FILE

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include <numeric>
#include <thread>

using namespace std;

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/exception/all.hpp> // Include for boost::exception
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/range/adaptor/reversed.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>

using namespace boost;

#include "ncl/nxsmultiformat.h"

#include "conditionals.hpp"
#include "xop.hpp"
#include "lot.hpp"
#include "split.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "tree_manip.hpp"
#include "tree_summary.hpp"
#include "opvertex.hpp"
#include "lapjv.hpp"
#include "op.hpp"

using namespace op;

// static data member initializations
mutex OP::_mutex;
string  OP::_program_name        = "op";

// If these are changed, also change in OPDistTest.cpp
unsigned     OP::_major_version       = 1;
//string     OP::_minor_version       = "4.0"; // fixes major bug related to nested common edges
//string     OP::_minor_version       = "5.0"; // fixes bug resulting from polytomies in one of the two trees being compared
//string     OP::_minor_version       = "6.0"; // now saves conf file as a record of settings used and fixes bug in TreeManip::nexusTranslateCommand
//string     OP::_minor_version       = "7.0"; // keep setting added; refdist now saves sorted distances as comments in a nexus-formatted tree file
//string     OP::_minor_version       = "8.0"; // added randwalk setting
string     OP::_minor_version      = "8.1";  // fixed stack overflow in pairwise

bool         OP::_silent              = false; // set to true only for unit tests
const double Node::_smallest_edge_length = 1.0e-12;
vector<string> TreeManip::_taxon_names;
map<string, unsigned> TreeManip::_taxon_map;

#if defined(OP_SAVE_DOT_FILE)
unsigned OP::_graph_number          = 1;
#endif

int main(int argc, const char * argv[]) {
    int return_code = 0;
    try {
        OP strom;
        strom.processCommandLineOptions(argc, argv);
        strom.run();
    }
    catch(std::exception & x) {
        cerr << "Exception: " << x.what() << endl;
        cerr << "Aborted." << endl;
        return_code = 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        return_code = 2;
    }

    return return_code;
}
