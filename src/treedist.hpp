#pragma once

#include <chrono>
#include <stdexcept>
#include <vector>
#include <map>
#include <fstream>
#include <numeric>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include "conditionals.hpp"
#include "output_manager.hpp"
//#include "lot.hpp"
#include "tree_summary.hpp"
#include "xtreedist.hpp"

extern Lot rng;

using namespace std;
using namespace boost;
//using math::quadrature::trapezoidal;

namespace treedist {

    class TreeDist {
        public:
                                        TreeDist();
                                        ~TreeDist();

            void                        clear();
            void                        processCommandLineOptions(int argc, const char * argv[]);
            void                        run();
            
        private:
        
            bool                        _quiet;
            bool                        _debug;
            bool                        _deroot;
            string                      _tree_file;
            string                      _out_file;
            string                      _ref_file;
            unsigned                    _ref_tree;
            unsigned                    _skip;
            
            // Program name and version
            static string               _program_name;
            static unsigned             _major_version;
            static unsigned             _minor_version;
    };
    
    inline TreeDist::TreeDist() {
        clear();
    }

    inline TreeDist::~TreeDist() {
    }

    inline void TreeDist::clear() {
        _deroot   = false;
        _quiet    = false;
        _debug    = false;
        _ref_tree = 1;
        _skip = 0;
    }
    
    inline void TreeDist::processCommandLineOptions(int argc, const char * argv[]) {
        program_options::variables_map vm;
        program_options::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version,v", "show program version")
            ("quiet,q", program_options::bool_switch(&_quiet), "minimize output")
            ("debug,d", program_options::bool_switch(&_debug), "copious output for debugging purposes")
            ("treefile,t", program_options::value(&_tree_file), "name of file containing trees to compare")
            ("outfile,o", program_options::value(&_out_file), "name of file in which to store distances")
            ("reffile,r", program_options::value(&_ref_file), "name of file containing reference tree")
            ("reftree", program_options::value(&_ref_tree), "index of tree to be used as the reference tree (default 1, the first tree in the file after any skipped trees)")
            ("skip,s", program_options::value(&_skip), "number of trees at the beginning of treefile to skip (default is 0)")
            ("deroot", program_options::value(&_deroot), "calculate distances on derooted trees")
        ;
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        try {
            const program_options::parsed_options & parsed = program_options::parse_config_file< char >("treedist.conf", desc, false);
            program_options::store(parsed, vm);
        }
        catch(program_options::reading_file & x) {
            //::om.outputConsole("Note: configuration file (treedist.conf) not found\n");
        }
        program_options::notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            ::om.outputConsole(desc);
            ::om.outputNewline();
            exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            ::om.outputConsole(format("This is %s version %d.%d\n") % _program_name % _major_version % _minor_version);
            exit(1);
        }
        
        if (_debug) {
            _quiet = true;
        }
    }
        
    inline void TreeDist::run() {
        TreeSummary ts;
        if (_ref_file != "") {
            // Read reference tree from _ref_file first
            ts.readTreefile(_ref_file,  /*skip*/0, /*ref_tree*/_ref_tree, /*store_all*/false, _deroot, _debug);
            
            // Now read in trees to compare to reference tree
            ts.readTreefile(_tree_file, /*skip*/_skip, /*ref_tree*/0, /*store_all*/true, _deroot, _debug);
        }
        else {
            // All trees (including reference tree) are in _tree_file
            ts.readTreefile(_tree_file, /*skip*/_skip, /*ref_tree*/_ref_tree,        /*store_all*/true, _deroot, _debug);
        }
        
        vector<double> kfvect;
        vector<unsigned> rfvect;
        ts.calcDistances(kfvect, rfvect, _quiet, _debug);
        
        string outfname = (_out_file == "" ? "dists.txt" : _out_file);
        ofstream outf(outfname);
        outf << "tree\tname\tKF\tRF\n";
        unsigned n = (unsigned)kfvect.size();
        assert(n == rfvect.size());
        for (unsigned i = 0; i < n; i++) {
            outf << str(format("%d\t%s\t%.9f\t%d\n") % (i+1) % ts.getTreeName(i) % kfvect[i] % rfvect[i]);
        }
        outf.close();
    }

}   // namespace TreeDist

