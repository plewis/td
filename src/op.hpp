#pragma once

// Not sure why the following include is necessary (it is also in main.cpp above the line including op.hpp)
// but if it is omitted, I get the following fatal compiler error:
// "error: use of undeclared identifier 'filesystem'; did you mean 'std::__fs::filesystem'?"
// Could use std::filesystem if compiled under c++17 instead of c++11
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost;

namespace op {
    class OP {
    public:
        OP();
        ~OP();

        void                clear();
        void                processCommandLineOptions(int argc, const char * argv[]);
        //void                processFakeCommandLineOptions(const vector<string> & args);
        void                run();

    private:

        void readTrees();
        //void assertFileExists(string fn) const;
        void shrinkToHPD();
        void outputForGTP() const;
        void moveAlongPath();
        void calcDistanceToReference() const;
        void calcFrechetMean();
        void calcPairwise() const;
        inline void randomWalk() const;
        void outputVersionAndSettings() const;
        unsigned buildThreadSchedule();

#if defined(TESTKDE)
        void testKDE();
#endif
        void calcBandWidth(const vector<double> & sample);
        double kernelDensity(double x, const vector<double> & sample) const;
        void rPlotDists(vector<double> & dists, string & rscript, double & radius, double & hpd_lower, double & hpd_upper, unsigned hpd_level);
        bool findEmpiricalHPDwaterline(double & hpd_cutoff, double & min_log_posterior, double & max_log_posterior) const;
        void buildRefTree(TreeManip & tm) const;
        void buildTree(unsigned tree_index, TreeManip & tm) const;
        void calcBHVDistanceForPairs(
            unsigned start_pair,
            unsigned end_pair,
            vector<TreeManip> & mu,
            vector<double> & bhv_distances,
            vector<pair<Split::treeid_t, Split::treeid_t> > & in_pairs,
            vector<Split::treeid_pair_t> & ABpairs,
            vector<Split::split_pair_t> & commonPairs) const;
        double calcBHVDistance(
            TreeManip & starttm,
            TreeManip & endtm,
            vector<pair<Split::treeid_t, Split::treeid_t> > & in_pairs,
            vector<Split::treeid_pair_t> & ABpairs,
            vector<Split::split_pair_t> & commonPairs) const;
        double calcClusterDistance(
                TreeManip & starttm,
                TreeManip & endtm,
                const vector<pair<Split::treeid_t, Split::treeid_t> > & in_pairs,
                const vector<Split::split_pair_t> & commonPairs) const;
        double calcKFDistance(unsigned ref_index, unsigned test_index) const;
        unsigned chooseRandomTree(TreeManip & tm, Lot & lot) const;
        void displaceTreeAlongGeodesic(TreeManip & start_tree, TreeManip & end_tree, double displacement) const;
        bool frechetCloseEnough(vector<TreeManip> & mu, unsigned lower, unsigned upper, double epsilon) const;
        unsigned computeFrechetMean(TreeManip & mean_tree) const ;
        static double opCalcTreeIDLength(
            const Split::treeid_t & splits);
        double opCalcLeafContribution(
            const Split::treeid_t & Alvs,
            const Split::treeid_t & Blvs,
            vector<Split::split_pair_t> & commonPairs) const;
        double opFindCommonEdges(
            const Split::treeid_t & A,
            const Split::treeid_t & B,
            vector<Split::split_pair_t> & commonPairs) const;
//        void opSplitAtCommonEdges(
//            const vector<Split> & common_edges,
//            vector<pair<Split::treeid_t,Split::treeid_t> > & in_pairs) const;
        void opSplitAtCommonEdges(
            const vector<Split::split_pair_t> & commonPairs,
            vector<pair<Split::treeid_t,Split::treeid_t> > & in_pairs) const;
#if defined(OP_SAVE_DOT_FILE)
        static string opCreateVertexLabel(string name, string capacity, string edgelen, string bipartition);
        static string opCreateEdgeLabel(double capacity, double reverse_flow);
        static void opSaveIncompatibilityGraph(
            OPVertex & source,
            OPVertex & sink,
            vector<OPVertex> & avect,
            vector<OPVertex> & bvect);
#endif
        static void opEdmondsKarp(
            OPVertex & source,
            OPVertex & sink,
            vector<OPVertex> & avect,
            vector<OPVertex> & bvect,
            Split::treeid_t & C1,
            Split::treeid_t & C2,
            Split::treeid_t & D1,
            Split::treeid_t & D2,
            bool quiet);
        bool opRefineSupport(
            const Split::treeid_pair_t & AB,
            Split::treeid_pair_t & AB1,
            Split::treeid_pair_t & AB2) const;
        double opCalcGeodesicDist(
            const Split::treeid_pair_t & inpair) const;

        bool                    _noisy;
        bool                    _output_for_gtp;
        string                  _prefix;
        bool                    _pairwise;
        bool                    _snapshot;
        bool                    _random_walk;
        bool                    _gtp;
        bool                    _refdist;
        bool                    _frechet_mean;
        //bool                    _save_credible_set;
        bool                    _save_trees;
        unsigned                _nthreads;
        unsigned                _precision;
        unsigned                _random_number_seed;
        vector<unsigned>        _skip;
        vector<unsigned>        _stride;
        vector<int>             _keep;
        vector<int>             _subsample;
        vector<unsigned>        _subseed;
        vector<bool>            _rooted;
        vector<string>          _tree_file_names;
        unsigned                _nsteps;
        double                  _step_mu;
        double                  _step_sigma;
        string                  _ref_tree_filename;
        unsigned                _refskip;
        bool                    _refrooted;
        double                  _refscale;
        vector<double>          _lambda;
        vector<double>          _scale_by;
        TreeSummary::SharedPtr  _tree_summary;

        double                  _frechet_epsilon;
        unsigned                _frechet_n;
        unsigned                _frechet_k;

        bool                    _hpd_radius;
        unsigned                _radius_percent;

        double                  _kde_bandwidth;
        double                  _kde_sigma;
        double                  _kde_q25;
        double                  _kde_q75;

        vector<pair<unsigned,unsigned> > _thread_sched;
        static mutex            _mutex;

        static string           _program_name;
        static unsigned         _major_version;
        static unsigned         _minor_version;

        static bool             _silent;    // set to true for unit tests only

#if defined(OP_SAVE_DOT_FILE)
        static unsigned            _graph_number;
#endif

    };

    inline OP::OP() :
        _noisy(false),
        _prefix("outfile"),
        _pairwise(false),
        _snapshot(false),
        _random_walk(false),
        _gtp(false),
        _refdist(false),
        _frechet_mean(false),
        //_save_credible_set(false),
        _save_trees(false),
        _nthreads(1),
        _precision(9),
        _nsteps(0),
        _step_mu(-0.05),
        _step_sigma(0.05),
        _random_number_seed(1),
        _ref_tree_filename(""),
        _refskip(0),
        _refrooted(false),
        _refscale(1.0),
        _tree_summary(nullptr),
        _frechet_epsilon(0.00001),
        _frechet_n(10),
        _frechet_k(1000000),
        _hpd_radius(false),
        _radius_percent(95),
        _kde_bandwidth(0.0),
        _kde_sigma(0.0),
        _kde_q25(0.0),
        _kde_q75(0.0) {
        TreeManip::_taxon_names.clear();
        TreeManip::_taxon_map.clear();
    }

    inline OP::~OP() = default;

    // inline void OP::processFakeCommandLineOptions(const vector<string> & args) {
    //     // Create C-style argc and argv from a vector of strings
    //     char **custom_argv = new const char *[args.size() + 1];
    //     for (size_t i = 0; i < args.size(); ++i) {
    //         custom_argv[i] = new char[args[i].length() + 1];
    //         strcpy(custom_argv[i], args[i].c_str());
    //     }
    //     custom_argv[args.size()] = nullptr;
    //     int custom_argc = args.size();
//
    //     processCommandLineOptions(custom_argc, custom_argv);
//
    //     // Delete memory
    //     for (size_t i = 0; i < args.size(); ++i) {
    //         delete[] custom_argv[i];
    //     }
    //     delete[] custom_argv;
    // }

    inline void OP::outputVersionAndSettings() const {
        if (_silent)
            return;

        cout << str(format("%s version %d.%d") % OP::_program_name % OP::_major_version % OP::_minor_version) << endl;

        boost::filesystem::path currentPath = boost::filesystem::current_path();
        cerr << "Current working directory: " << currentPath << endl;

        string fn = _prefix + "-op.conf";
        cout << "Settings saved in file \"" << fn << "\"" << endl;
        cout << "Rename \"" << fn << " to \"op.conf\" to recreate this analysis." << endl;
        ofstream outf(fn);
        outf << str(format("# Configuration file for %s version %d.%d\n") % _program_name % _major_version % _minor_version);
        outf << "# https://github.com/plewis/op/\n";
        outf << str(format("prefix = %s\n") % _prefix);
        if (_noisy) {
            outf << "noisy = yes\n";
        }
        if (_hpd_radius) {
            outf << "hpdrad = yes\n";
        }
        outf << "radpct = " << _radius_percent << "\n";
        if (_save_trees) {
            outf << "savetrees = yes\n";
        }
        outf << str(format("seed = %d\n") % _random_number_seed);
        if (_random_walk) {
            outf << str(format("nsteps = %d\n") % _nsteps);
            outf << str(format("stepmu = %.9f\n") % _step_mu);
            outf << str(format("stepsigma = %.9f\n") % _step_sigma);
        }
        else if (_pairwise) {
            outf << "pairwise = yes\n";
        }
        else if (_refdist) {
            //outf << "refdist = yes\n";
            outf << "reftree = " << _ref_tree_filename << endl;
            outf << "refskip = " << _refskip << endl;
            outf << "refrooted = " << (_refrooted ? "yes" : "no") << endl;
            outf << str(format("refscale = %.9f") % _refscale) << endl;
        }
        else if (_gtp) {
            outf << "saveforgtp = yes\n";
        }
        else if (_frechet_mean) {
            outf << "frechetmean = yes\n";
            outf << str(format("frechet-e = %.9f\n") % _frechet_epsilon);
            outf << str(format("frechet-k = %d\n") % _frechet_k);
            outf << str(format("frechet-n = %d\n") % _frechet_n);
            outf << str(format("hpdrad = %s\n") % (_hpd_radius ? "yes" : "no"));
            outf << str(format("radpct = %d\n") % _radius_percent);
        }
        else {
            assert(_snapshot);
            for (auto & l : _lambda) {
                outf << str(format("lambda = %.9f\n") % l);
            }
        }
        for (unsigned i = 0; i < _tree_file_names.size(); ++i) {
            outf << "treefile = " << _tree_file_names[i] << endl;
            outf << "skip = " << _skip[i] << endl;
            outf << "stride = " << _stride[i] << endl;
            outf << "keep = " << _keep[i] << endl;
            outf << "subsample = " << _subsample[i] << endl;
            outf << "subseed = " << _subseed[i] << endl;
            outf << "rooted = " << (_rooted[i] ? "yes" : "no") << endl;
            outf << str(format("scale = %.9f\n") % _scale_by[i]);
        }
        outf.close();
    }

    inline void OP::processCommandLineOptions(int argc, const char * argv[]) {
        program_options::variables_map       vm;
        program_options::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version,v", "show program version")
            ("noisy",program_options::bool_switch(&_noisy), "show a lot of output (default: no)")
            ("reftree",  program_options::value(&_ref_tree_filename), "name of reference tree file in NEXUS format; output will be distances from this reference tree to all trees specified with --treefile (no default)")
            ("refskip", program_options::value(&_refskip), "number of trees to skip in specified reference tree file (default: 0)")
            ("refrooted", program_options::value(&_refrooted), "assume reference tree is rooted (default: no)")
            ("refscale", program_options::value(&_refscale), "rescale reference tree by this multiplicative factor (default: 1.0)")
            ("treefile,t",  program_options::value(&_tree_file_names), "name of tree file in NEXUS format (no default)")
            ("skip", program_options::value(&_skip), "number of trees to skip in specified treefile (default: 0)")
            ("stride", program_options::value(&_stride), "number of trees to skip before saving (default: 1)")
            ("keep", program_options::value(&_keep), "number of trees to retain from treefile (default: all)")
            ("subsample", program_options::value(&_subsample), "number of trees to randomly retain (default: no random subsampling)")
            ("subseed", program_options::value(&_subseed), "random number seed to use for subsampling (default: 1)")
            ("rooted", program_options::value(&_rooted), "assume trees are rooted (default: no)")
            ("precision", program_options::value(&_precision)->default_value(9), "number of digits precision to use in outputting distances (default: 9)")
            ("prefix", program_options::value(&_prefix), "filename prefix for output file name (default: 'outfile')")
            ("lambda", program_options::value(&_lambda), "specify a value in [0,1] to calculate tree at that point (assumes starting tree is first tree and ending tree is the second tree in the treefile)")
            ("nsteps", program_options::value(&_nsteps), "specifies number of steps to take in random walk")
            ("stepmu", program_options::value(&_step_mu), "specifies mean of normal variates used in random walk")
            ("stepsigma", program_options::value(&_step_sigma), "specifies standard deviation of normal variates used in random walk")
            ("pairwise", program_options::bool_switch(&_pairwise),"calculates pairwise distances (default: pairwise distances not calculated)")
            //("refdist", program_options::bool_switch(&_refdist),"calculates distance of the first tree to all other trees (default: distances not calculated)")
            ("frechetmean", program_options::bool_switch(&_frechet_mean),"calculate Frechet mean tree and variance (default: mean and variance not calculated)")
            ("frechet-e,e", program_options::value(&_frechet_epsilon), "successive Frechet mean approximations must all be this close to stop iterating (default: 0.00001)")
            ("frechet-n,n", program_options::value(&_frechet_n), "number of successive Frechet mean approximations to use for determining whether to stop iterating (default: 10)")
            ("frechet-k,k", program_options::value(&_frechet_k), "maximum number of Frechet mean iterations (default:1000000)")
            ("hpdrad",program_options::bool_switch(&_hpd_radius), "radius includes all trees within the radpct HPD credible set (default is that radius includes a percentage radpct of trees without regard to posterior)")
            ("radpct", program_options::value(&_radius_percent), "radius includes this fraction of sampled trees closest to mean (default: 0.95)")
            ("seed", program_options::value(&_random_number_seed), "pseudorandom number generator seed (used only when estimating mean tree)")
            ("scale", program_options::value(&_scale_by), "rescale all input trees by this multiplicative factor (default: 1.0)")
            ("nthreads", program_options::value(&_nthreads), "number of threads to use in determining Frechet mean (default: 1)")
            ("saveforgtp", program_options::bool_switch(&_gtp),"saves trees in format interpretable by Owen & Provan GTP java program)")
            ("savetrees","save trees after any HPD filtering and subsampling has taken place in nexus format with \"-op.tre\" extension (default: no)")
            //("credset","if posterior available (i.e. treefile from RevBayes), output credible set of trees and extreme trees (default: no)")
#if defined(TESTKDE)
            ("testkde", "test kernel density estimation")
#endif
;
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        try {
            const program_options::parsed_options & parsed = program_options::parse_config_file< char >("op.conf", desc, false);
            program_options::store(parsed, vm);
        }
        catch(program_options::reading_file &) {
            if (_noisy) {
                cout << "Note: configuration file (op.conf) not found" << endl;
            }
        }
        program_options::notify(vm);

        // If the user specified --help on the command line, output usage summary and quit
        if (vm.count("help") > 0) {
            cout << desc << "\n";
            exit(1);
        }

        // If the user specified --version on the command line, output the version and quit
        if (vm.count("version") > 0) {
            cout << str(format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << endl;
            exit(0);
        }

        // // If the user specified --noisy on the command line, set _noisy to true
        // if (vm.count("noisy") > 0) {
        //     _noisy = true;
        // }

#if defined(TESTKDE)
        // If the user specified --testkde on the command line, run testKDE() and quit
        if (vm.count("testkde") > 0) {
            testKDE();
            exit(0);
        }
#endif

        // If the user specified --prefix on the command line, set _prefix
        if (vm.count("prefix") > 0) {
            _prefix = vm["prefix"].as<string>();
        }

        // If the user specified --reftree on the command line, set _refdist to true
        if (vm.count("reftree") > 0) {
            // These will be reversed if nsteps > 0
            _refdist = true;
            _random_walk = false;
        }

        // If the user specified --nsteps on the command line, set _random_walk to true
        // and sanity-check the value of nsteps supplied
        if (vm.count("nsteps") > 0) {
            if (_nsteps > 0 && _refdist) {
                _random_walk = true;
                _refdist = false;
                if (_nsteps > 1000000) {
                    throw Xop(format("nsteps should be less than 1 million, but you specified %d") % _nsteps);
                }
                if (_step_sigma <= 0.0) {
                    throw Xop(format("stepsigma should be greater than zero, but you specified %g") % _step_sigma);
                }
            }
        }

        // If the user failed to specify --treefile on the command line, bail out because a treefile is needed for
        // anything except help and version
        if (vm.count("treefile") == 0 && !_random_walk) {
            cout << "You must specify a treefile if doing anything except --help, --version, or a random walk" << endl;
            exit(1);
        }

        // If the user specified --lambda on the command line, set _snapshot to true
        // and sanity-check the value of lambda supplied
        if (vm.count("lambda") > 0) {
            _snapshot = true;
            sort(_lambda.begin(), _lambda.end());
            for (auto & l : _lambda) {
                if (l < 0.0) {
                    throw Xop(format("Lambda must be greater than or equal to 0.0, but you specified %g") % l);
                } else if (l > 1.0) {
                    throw Xop(format("Lambda must be less than or equal to 1.0, but you specified %g") % l);
                }
            }
        }

        // // If the user specified --frechet on the command line, set _frechet_mean to true
        // if (vm.count("credset") > 0) {
        //     _save_credible_set = true;
        // }

        // If the user specified --frechet on the command line, set _frechet_mean to true
        if (vm.count("savetrees") > 0) {
            _save_trees = true;
        }

        if (_nthreads > 1)
            buildThreadSchedule();

        // Sanity check
        bool ok = ( _pairwise && !_frechet_mean && !_refdist && !_snapshot && !_gtp && !_random_walk);    // only pairwise chosen
        ok     |= (!_pairwise &&  _frechet_mean && !_refdist && !_snapshot && !_gtp && !_random_walk);    // only frechet mean chosen
        ok     |= (!_pairwise && !_frechet_mean &&  _refdist && !_snapshot && !_gtp && !_random_walk);    // only refdist chosen
        ok     |= (!_pairwise && !_frechet_mean && !_refdist &&  _snapshot && !_gtp && !_random_walk);    // only snapshot chosen
        ok     |= (!_pairwise && !_frechet_mean && !_refdist && !_snapshot &&  _gtp && !_random_walk);    // only gtp chosen
        ok     |= (!_pairwise && !_frechet_mean && !_refdist && !_snapshot && !_gtp &&  _random_walk);    // only random walk chosen

        if (!ok) {
            cerr << "Sorry, you can only do one of these things during a single run: " << endl;
            cerr << "(1) calculate pairwise distances (pairwise); " << endl;
            cerr << "(2) calculate distance from a reference tree to all other trees (refdist is yes if reftree is specified);" << endl;
            cerr << "(3) calculate frechet mean and variance (frechet); or" << endl;
            cerr << "(4) move along path a specified distance (lambda)" << endl;
            cerr << "(5) carry out random walk of length nsteps" << endl;
            cerr << "(6) save trees for use with Owen & Provan GTP software" << endl;
            cerr << "You specified: " << endl;
            cerr << "  pairwise    = " << (_pairwise ? "yes" : "no") << endl;
            cerr << "  frechetmean = " << (_frechet_mean ? "yes" : "no") << endl;
            cerr << "  refdist     = " << (_refdist ? "yes" : "no") << endl;
            cerr << "  lambda      = " << (_snapshot ? "yes" : "no") << endl;
            cerr << "  nsteps > 0  = " << (_random_walk ? "yes" : "no")  << endl;
            cerr << "  saveforgtp  = " << (_gtp ? "yes" : "no") << endl;
            throw Xop("Program aborted.");
        }

        if (!_random_walk) {
            if (_skip.empty()) {
                // no trees are skipped by default
                _skip.resize(_tree_file_names.size(), 0);
            }
            else if (_skip.size() != _tree_file_names.size()) {
                throw Xop("If specified for any treefile, a skip setting must be provided for every treefile");
            }

            if (_stride.empty()) {
                // stride is 1 by default
                _stride.resize(_tree_file_names.size(), 1);
            }
            else if (_stride.size() != _tree_file_names.size()) {
                throw Xop("If specified for any treefile, a stride setting must be provided for every treefile");
            }

            if (_keep.empty()) {
                // keep is -1 by default, which results in saving every tree
                _keep.resize(_tree_file_names.size(), -1);
            }
            else if (_keep.size() != _tree_file_names.size()) {
                throw Xop("If specified for any treefile, a keep setting must be provided for every treefile");
            }

            if (_subsample.empty()) {
                // subsample is -1 by default, which results in saving every tree (i.e. no subsampling)
                _subsample.resize(_tree_file_names.size(), -1);
            }
            else if (_subsample.size() != _tree_file_names.size()) {
                throw Xop("If specified for any treefile, a subsample setting must be provided for every treefile");
            }

            if (_subseed.empty()) {
                // subseed is 1 by default
                _subseed.resize(_tree_file_names.size(), 1);
            }
            else if (_subseed.size() != _tree_file_names.size()) {
                throw Xop("If specified for any treefile, a subseed setting must be provided for every treefile");
            }

            if (_rooted.empty()) {
                // trees are assumed to be unrooted by default
                _rooted.resize(_tree_file_names.size(), false);
            }
            else if (_rooted.size() != _tree_file_names.size()) {
                throw Xop("If specified for any treefile, rooted setting must be provided for every treefile");
            }

            if (_scale_by.empty()) {
                _scale_by.resize(_tree_file_names.size(), 1.0);
            }
            else if (_scale_by.size() != _tree_file_names.size()) {
                throw Xop("If specified, a scale setting must be provided for every treefile");
            }
        }
    }

    inline double OP::opCalcTreeIDLength(const Split::treeid_t & splits) {
        double length = 0.0;
        for (auto & split : splits) {
            length += pow(split.getEdgeLen(),2);
        }
        return sqrt(length);
    }

    inline double OP::opCalcLeafContribution(
                const Split::treeid_t & Alvs,
                const Split::treeid_t & Blvs,
                vector<Split::split_pair_t> & commonPairs) const {
        if (_noisy) {
            cout << "\nLeaves from starting tree:" << endl;
            for (auto & a : Alvs) {
                cout << "  " << a.createPatternRepresentation(true) << endl;
            }

            cout << "\nLeaves from ending tree:" << endl;
            for (auto & b : Blvs) {
                cout << "  " << b.createPatternRepresentation(true) << endl;
            }
        }

        // Compute leaf contribution
        double leaf_contribution_squared = 0.0;
        for (auto & b : Blvs) {
            auto it = find(Alvs.begin(), Alvs.end(), b);
            assert(it != Alvs.end());
            double leafa = it->getEdgeLen();
            double leafb = b.getEdgeLen();
            leaf_contribution_squared += pow(leafa-leafb, 2);
            commonPairs.emplace_back(*it, b);
        }

        if (_noisy)
            cout << str(format("\nLeaf contribution (squared) = %.9f") % leaf_contribution_squared) << endl;
        return leaf_contribution_squared;
    }

    inline double OP::opFindCommonEdges(
                const Split::treeid_t & A,
                const Split::treeid_t & B,
                vector<Split::split_pair_t> & commonPairs) const {
        // Find splits in the intersection of A and B
        vector<Split> common_edges;
        set_intersection(
            A.begin(), A.end(),
            B.begin(), B.end(),
            back_inserter(common_edges)
        );

        if (_noisy) {
            cout << "\nCommon edges:" << endl;
        }

        double common_edge_contribution_squared = 0.0;
        for (auto & s : common_edges) {
            auto itA = find(A.begin(), A.end(), s);
            double edgeA = itA->getEdgeLen();
            auto itB = find(B.begin(), B.end(), s);
            double edgeB = itB->getEdgeLen();
            common_edge_contribution_squared += pow(edgeA-edgeB, 2);
            commonPairs.emplace_back(*itA, *itB);
            if (_noisy) {
                double diff = edgeB-edgeA;
                cout << boost::str(format("%12.9f: %s\n") % diff % itA->createPatternRepresentation(false));
            }
        }

        if (_noisy) {
            cout << "\nCommon edge contribution (squared): " << setprecision(9) << common_edge_contribution_squared << endl;
        }

        // Count the number of splits in B compatible with each split in A (and vice versa)
        map<const Split *, unsigned> acompatibilities;
        map<const Split *, unsigned> bcompatibilities;
        for (auto & a : A) {
            for (auto & b : B) {
                if (a.compatibleWith(b)) {
                    acompatibilities[&a]++;
                    bcompatibilities[&b]++;
                }
            }
        }

        if (_noisy) {
            cout << "\nEdges in starting tree compatible with all splits in ending tree:" << endl;
        }

        // Add splits in A that are compatible with all splits in B to common_edges
        unsigned nAcompatB = 0;
        for (auto & apair : acompatibilities) {
            if (apair.second == B.size()) {
                // This split is compatible with every split in the other tree, so add it to the vector
                // of common edges if it is not already in the vector of common edges
                if (find(common_edges.begin(), common_edges.end(), *apair.first) == common_edges.end()) {
                    common_edges.push_back(*apair.first);
                    common_edge_contribution_squared += pow(apair.first->getEdgeLen(), 2);
                    commonPairs.emplace_back(*(apair.first), Split());
                    nAcompatB++;
                    if (_noisy) {
                        cout << apair.first->createPatternRepresentation(true) << endl;
                    }
                }
            }
        }

        if (_noisy) {
            if (nAcompatB == 0) {
                cout << "  None found." << endl;
            }
            cout << "\nEdges in ending tree compatible with all splits in starting tree:" << endl;
        }

        // Add splits in B that are compatible with all splits in A to common_edges
        unsigned nBcompatA = 0;
        for (auto & bpair : bcompatibilities) {
            if (bpair.second == A.size()) {
                // This split is compatible with every split in the other tree, so add it to the vector
                // of common edges if it is not already in the vector of common edges
                if (find(common_edges.begin(), common_edges.end(), *bpair.first) == common_edges.end()) {
                    common_edges.push_back(*bpair.first);
                    common_edge_contribution_squared += pow(bpair.first->getEdgeLen(), 2);
                    commonPairs.emplace_back(Split(), *(bpair.first));
                    nBcompatA++;
                    if (_noisy) {
                        cout << bpair.first->createPatternRepresentation(true) << endl;
                    }
                }
            }
        }

        if (_noisy) {
            if (nBcompatA == 0) {
                cout << "  None found." << endl;
            }
            if (nAcompatB + nBcompatA > 0) {
                cout << "\nCommon edge contribution (squared): " << setprecision(9) << common_edge_contribution_squared << endl;
            }
        }

        return common_edge_contribution_squared;
    }

    inline void OP::opSplitAtCommonEdges(const vector<Split::split_pair_t> & commonPairs, vector<pair<Split::treeid_t,Split::treeid_t> > & in_pairs) const {
        // Bail out now if there are no common edges
        if (commonPairs.size() == 0) {
            return;
        }

        // Expecting just one pair of split sets (treeIDs) to start
        assert(in_pairs.size() == 1);

        // Remove A and B from in_pairs
        vector<Split> A(in_pairs[0].first.begin(), in_pairs[0].first.end());
        vector<Split> B(in_pairs[0].second.begin(), in_pairs[0].second.end());
        in_pairs.clear();

        // Store common edges in a vector
        //TODO: should this be done from the start? That is, why have a vector of pairs when first == second?
        vector<Split> common_edges;
        for (auto & splitpair : commonPairs) {
            if (splitpair.first.getSize() > 0)
                common_edges.push_back(splitpair.first);
            else if (splitpair.second.getSize() > 0)
                common_edges.push_back(splitpair.second);
            else {
                throw Xop("Internal error: neither member of splitpair has size > 0");
            }
        }

        // For each common split, segregate all subsumed splits in A and B into separate in_pairs elements
        for (auto & common : common_edges) {
            unsigned nA = static_cast<unsigned>(A.size());
            unsigned nB = static_cast<unsigned>(B.size());
            pair<Split::treeid_t,Split::treeid_t> in_pair;
            stack<unsigned> erased;

            // Move splits in A to in_pair.first if they are subsumed in common
            for (unsigned i = 0; i < nA; i++) {
                if (A[i].subsumedIn(common)) {
                    // A[i] is subsumed in common
                    if (A[i] != common) {
                        in_pair.first.insert(A[i]);
                    }
                    erased.push(i);
                }
            }

            // Erase splits in common from A (note: starting at the end
            // so the indexing will not be changed by deletions)
            while (!erased.empty()) {
                unsigned i = erased.top();
                erased.pop();
                A.erase(A.begin() + i);
            }

            // Move splits in B to in_pair.second if they are subsumed in common
            for (unsigned i = 0; i < nB; i++) {
                if (B[i].subsumedIn(common)) {
                    // B[i] is subsumed in common
                    if (B[i] != common) {
                        in_pair.second.insert(B[i]);
                    }
                    erased.push(i);
                }
            }

            // Erase splits in common from B
            while (!erased.empty()) {
                unsigned i = erased.top();
                erased.pop();
                B.erase(B.begin() + i);
            }

            // Only save pair if at least one split is in each set
            if (in_pair.first.size() > 0 && in_pair.second.size() > 0)
                in_pairs.emplace_back(in_pair);
        }

        // Insert what remains of A and B into in_pairs if more than one split is in each set
        pair<Split::treeid_t,Split::treeid_t> last_pair;
        for (auto & split : A) {
            last_pair.first.insert(split);
        }
        for (auto & split : B) {
            last_pair.second.insert(split);
        }
        if (last_pair.first.size() > 0 && last_pair.second.size() > 0)
            in_pairs.emplace_back(last_pair);
    }

#if defined(OP_SAVE_DOT_FILE)
    inline string OP::opCreateVertexLabel(string name, string capacity, string edgelen, string bipartition) {
        string s = "";
        s += "<<table bgcolor=\"white\" border=\"0\">";
        s += "<tr><td><b><font color=\"blue\" face=\"Courier\" point-size=\"16\">%s</font></b></td></tr>";
        s += "<tr><td><b><font color=\"blue\" face=\"Courier\" point-size=\"16\">%s</font></b></td></tr>";
        s += "<tr><td><font color=\"black\" face=\"Courier\" point-size=\"10\">%s</font></td></tr>";
        s += "<tr><td><font color=\"black\" face=\"Courier\" point-size=\"10\">%s</font></td></tr>";
        s += "</table>>";
        return str(format(s) % name % capacity % bipartition % edgelen);
    }

    inline string OP::opCreateEdgeLabel(double capacity, double reverse_flow) {
        string s = "";
        s += "\tlabeldistance=8\n";
        s += "\tlabelangle=0\n";
        s += "\theadlabel=<<font color=\"red\" face=\"Verdana\" point-size=\"12\">%.3f</font>>\n";
        return str(format(s) % reverse_flow);
    }

    inline void OP::opSaveIncompatibilityGraph(OPVertex & source, OPVertex & sink, vector<OPVertex> & avect, vector<OPVertex> & bvect) {
        // Example of the kind of dot file generated by this function:
        // digraph G {
        //     rankdir=LR;
        //     graph [ranksep=2];
        //
        //     subgraph A_vertices {
        //         label="A";
        //             {
        //             rank=same;
        //             a1 [label = "0.349\na1\n----**-", shape = box, color = black];
        //             a2 [label = "0.100\na2\n--*-**-", shape = box, color = black];
        //             a3 [label = "0.240\na3\n-**-**-", shape = box, color = black];
        //             }
        //         a1 -> a2 [dir=none, style=invisible];
        //         a2 -> a3 [dir=none, style=invisible];
        //     }
        //
        //     subgraph B_vertices {
        //         label="B";
        //             {
        //             rank=same;
        //             b1 [label = "0.016\nb1\n--*-*--", shape = box, color = black];
        //             b2 [label = "0.524\nb2\n-**-*--", shape = box, color = black];
        //             b3 [label = "0.122\nb3\n-----**", shape = box, color = black];
        //             b4 [label = "0.339\nb4\n-**-***", shape = box, color = black];
        //             }
        //         b1 -> b2 [dir=none, style=invisible];
        //         b2 -> b3 [dir=none, style=invisible];
        //         b3 -> b4 [dir=none, style=invisible];
        //     }
        //
        //     a1 -> b1;
        //     a1 -> b2;
        //     a1 -> b3;
        //     a2 -> b2;
        //     a2 -> b3;
        //     a3 -> b3;
        //     a3 -> b4;
        // }
        //
        //

        auto asize = static_cast<unsigned>(avect.size());
        auto bsize = static_cast<unsigned>(bvect.size());
        unsigned minsize = min(asize, bsize);
        unsigned maxsize = max(asize, bsize);

        // Save all the entities needed
        // tuple key: <0> name, <1>capacity, <2>edgelen, <3>split, <4>shape, <5>color
        typedef std::tuple<string, string, string, string, string, string> gnode_element_t;
        typedef vector<gnode_element_t> gnode_t;

        // Save everything needed for the A nodes
        gnode_t anodes;
        for (unsigned i = 0; i < source._edges.size(); i++) {
            OPEdge * edge   = source._edges[i];
            string name     = edge->_to->_name;
            string capacity = (edge->_capacity == 0.0 ? "0" : str(format("%.3f") % edge->_capacity));
            string edgelen  = str(format("%.3f") % edge->_to->_split->getEdgeLen());
            string split    = str(format("%.3f") % edge->_to->_split->createPatternRepresentation(false));
            string shape    = (edge->_capacity == 0.0 ? "circle" : "box");
            string color    = (edge->_capacity == 0.0 ? "red" : "black");
            anodes.emplace_back(name, capacity, edgelen, split, shape, color);
        }
        if (asize < maxsize) {
            for (unsigned i = asize; i < maxsize; i++) {
                anodes.emplace_back(str(format("adummy%d") % (i+1)), "0", "0", "0", "box", "black");
            }
        }

        // Save everything needed for the B nodes
        gnode_t bnodes;
        for (unsigned i = 0; i < bsize; i++) {
            assert(bvect[i]._edges.size() == 1);
            OPEdge * edge   = bvect[i]._edges[0];
            string name     = edge->_from->_name;
            string capacity = (edge->_capacity == 0.0 ? "0" : str(format("%.3f") % edge->_capacity));
            string edgelen  = str(format("%.3f") % bvect[i]._split->getEdgeLen());
            string split    = str(format("%.3f") % bvect[i]._split->createPatternRepresentation(false));
            string shape    = (edge->_capacity == 0.0 ? "circle" : "box");
            string color    = (edge->_capacity == 0.0 ? "red" : "black");
            bnodes.emplace_back(name, capacity, edgelen, split, shape, color);
        }
        if (bsize < maxsize) {
            for (unsigned i = bsize; i < maxsize; i++) {
                bnodes.emplace_back(str(format("bdummy%d") % (i+1)), "0", "0", "0", "box", "black");
            }
        }

        // Create (or append to) a rundot.sh file containing commands to compile all the incompatibility graphs
        if (OP::_graph_number == 1) {
            ofstream shf("rundot.sh");
            shf << "#!/bin/bash\n\n";
            shf << "# This file requires previous installation of dot and is intended\n";
            shf << "# to be run on a Mac (because of the use of the open command)\n\n";
            shf << "dot -Tpng graph-1.dot > graph-1.png; open graph-1.png\n";
            shf.close();
        }
        else {
            ofstream shf("rundot.sh", ios::out | ios::app);
            shf << str(format("dot -Tpng graph-%d.dot > graph-%d.png; open graph-%d.png\n")
                % OP::_graph_number
                % OP::_graph_number
                % OP::_graph_number);
            shf.close();
        }

        // Open the dot file
        cout << "Saving incompatibility graph to file: " << str(format("graph-%d.dot") % OP::_graph_number) << endl;
        ofstream dotf(str(format("graph-%d.dot") % OP::_graph_number++));

        // Create the opening preamble
        dotf << "digraph G {\n";
        dotf << "\trankdir=LR;\n";
        dotf << "\tgraph [ranksep=2];\n\n";

        // Create a subgraph containing only the A vertices in a vertical rank
        dotf << "\tsubgraph Avertices {\n";
        dotf << "\t    label=\"A\";\n";
        dotf << "\t    {\n";
        dotf << "\t        rank=same;\n";
        for (unsigned i = 0; i < maxsize; i++) {
            if (i < asize) {
                // tuple key: <0> name, <1>capacity, <2>edgelen, <3>split, <4>shape, <5>color

                dotf << str(format("\t        %s [label = %s, shape = %s, color = %s];\n")
                    % get<0>(anodes[i])
                    % opCreateVertexLabel(get<0>(anodes[i]), get<1>(anodes[i]), get<2>(anodes[i]), get<3>(anodes[i]))
                    % get<4>(anodes[i])
                    % get<5>(anodes[i])
                );
            }
            else {
                dotf << str(format("\t        %s [style = invisible];\n")
                    % get<0>(anodes[i])
                );
            }
        }
        dotf << "\t    }\n";
        for (unsigned i = 0; i < maxsize - 1; i++) {
            dotf << str(format("\t    %s -> %s [dir=none, style=invisible];\n")
                    % get<0>(anodes[i])
                    % get<0>(anodes[i+1])
            );
        }
        dotf << "\t}\n";

        // Create a subgraph containing only the B vertices in a vertical rank
        dotf << "\tsubgraph Bvertices {\n";
        dotf << "\t    label=\"B\";\n";
        dotf << "\t    {\n";
        dotf << "\t        rank=same;\n";
        for (unsigned i = 0; i < maxsize; i++) {
            if (i < bsize) {
                dotf << str(format("\t        %s [label = %s, shape = %s, color = %s];\n")
                    % get<0>(bnodes[i])
                    % opCreateVertexLabel(get<0>(bnodes[i]), get<1>(bnodes[i]), get<2>(bnodes[i]), get<3>(bnodes[i]))
                    % get<4>(bnodes[i])
                    % get<5>(bnodes[i])
                );
            }
            else {
                dotf << str(format("\t        %s [style = invisible];\n")
                    % get<0>(bnodes[i])
                );
            }
        }
        dotf << "\t    }\n";
        for (unsigned i = 0; i < maxsize - 1; i++) {
            dotf << str(format("\t    %s -> %s [dir=none, style=invisible];\n")
                    % get<0>(bnodes[i])
                    % get<0>(bnodes[i+1])
            );
        }
        if (bsize < maxsize) {
            dotf << str(format("\t    %s -> bdummy%d [dir=none, style=invisible];\n")
                % get<0>(anodes[bsize-1])
                % (bsize+1));
            for (unsigned i = bsize; i < maxsize - 1; i++) {
                dotf << str(format("\t    bdummy%d -> bdummy%d [dir=none, style=invisible];\n") % (i+1) % (i+2));
            }
        }
        dotf << "\t}\n";

        // Create the edges connecting A with B vertices
        for (unsigned i = 0; i < asize; ++i) {
            for (unsigned j = 0; j < avect[i]._edges.size(); ++j) {
                OPEdge * edge = avect[i]._edges[j];
                dotf << str(format("\t\t%s -> %s [%s];\n") % edge->_from->_name % edge->_to->_name % opCreateEdgeLabel(edge->_capacity, edge->_reverse_flow));
            }
        }

        // dotf << "\t}\n";
        dotf << "}\n";
        dotf.close();
    }
#endif

    inline void OP::opEdmondsKarp(
            OPVertex & source,
            OPVertex & sink,
            vector<OPVertex> & avect,
            vector<OPVertex> & bvect,
            Split::treeid_t & C1,
            Split::treeid_t & C2,
            Split::treeid_t & D1,
            Split::treeid_t & D2,
            bool quiet) {
#if defined(OP_SAVE_DOT_FILE)
        opSaveIncompatibilityGraph(source, sink, avect, bvect);
#endif

        double cumulative_flow = 0.0;
        bool done_augmenting_path = false;
        while (!done_augmenting_path) {
            vector<OPVertex *> route;

            // All vertices begin as unvisited, have residual capacity 0, and all "A" to "B" edges begin as non-reversed
            source._parent_edge = nullptr;
            for (auto & a : avect) {
                a._parent_edge = nullptr;
                a._residual_capacity = 0.0;
                for (auto & a_edge : a._edges) {
                    a_edge->_edge_is_reversed = false;
                }
            }
            for (auto & b : bvect) {
                b._parent_edge = nullptr;
            }
            sink._parent_edge = nullptr;

            // Add the source to the route
            route.push_back(&source);
            unsigned route_cursor = 0;

            // Find the next augmenting route
            bool sink_found = false;
            while (!sink_found) {
                OPVertex * current = route[route_cursor];
                for (auto & edge : current->_edges) {
                    if (edge->_capacity > 0.0 && edge->_to->_parent_edge == nullptr) {
                        edge->_to->_parent_edge = edge;

                        route.push_back(edge->_to);
                        if (edge->_to == &sink) {
                            sink_found = true;
                            break;
                        }
                    }
                    else {
                        bool already_visited = static_cast<bool>(edge->_to->_parent_edge != nullptr);
                        if (!already_visited) {
                            // If it is a "B" vertex that has zero capacity, then the route cannot
                            // go from this "B" vertex to the sink. It may, however, be able to go back to
                            // an "A" vertex if that "A" vertex was not accessible from the source
                            // and there is residual flow on the edge.
                            if (edge->_to == &sink) {
                                OPVertex * B_vertex = edge->_from;
                                for (auto s_to_A_edge : source._edges) {
                                    // edge from source to an "A" vertex
                                    if (s_to_A_edge->_capacity == 0) {
                                        // See if any of the "A" vertex edges lead to the "B" vertex in question
                                        OPVertex * A_vertex = s_to_A_edge->_to;
                                        for (unsigned k = 0; k < A_vertex->_edges.size(); k++) {
                                            OPEdge * A_to_B_edge = A_vertex->_edges[k];  // edge from "A" to a "B" vertex
                                            if (A_vertex->_parent_edge == nullptr && A_to_B_edge->_to == B_vertex && A_to_B_edge->_reverse_flow > 0.0) {
                                                // If we are here, we know that:
                                                // 1. the "A" vertex has not been visited
                                                // 2. the "A" vertex is joined to the "B" vertex in question
                                                // 3. the "A" vertex is not accessible from the source, and
                                                // 4. there is residual flow on the A_to_B_edge

                                                A_vertex->_parent_edge = A_to_B_edge;
                                                A_to_B_edge->_edge_is_reversed = true;
                                                A_vertex->_residual_capacity = A_to_B_edge->_reverse_flow;
                                                route.push_back(A_vertex);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                route_cursor++;
                if (route_cursor >= route.size()) {
                    // Not able to reach sink, so all augmenting routes have been found
                    done_augmenting_path = true;
                    break;
                }
            }

            if (!done_augmenting_path) {
                // Follow _from pointers from sink back to source to identify
                // the bottleneck and determine flow through this route
                double min_capacity = 1.0;
                OPVertex * current = &sink;
                while (current != &source) {
                    OPEdge * edge = current->_parent_edge;
                    if (edge->_edge_is_reversed) {
                        if (edge->_reverse_flow < min_capacity) {
                            min_capacity = edge->_reverse_flow;
                        }
                        current = edge->_to;
                    }
                    else {
                        if (edge->_capacity < min_capacity) {
                            min_capacity = edge->_capacity;
                        }
                        current = edge->_from;
                    }
                }

                // Adjust capacities and flows along the route
                cumulative_flow += min_capacity;
                current = &sink;
                while (current != &source) {
                    OPEdge * edge = current->_parent_edge;
                    if (edge->_edge_is_reversed) {
                        edge->_reverse_flow -= min_capacity;
                        if (fabs(edge->_reverse_flow) < 1e-10)
                            edge->_reverse_flow = 0.0;
                        edge->_capacity += min_capacity;
                        current = edge->_to;
                    }
                    else {
                        edge->_capacity -= min_capacity;
                        if (fabs(edge->_capacity) < 1e-10)
                            edge->_capacity = 0.0;
                        edge->_reverse_flow += min_capacity;
                        current = edge->_from;
                    }
                }
            }
#if defined(OP_SAVE_DOT_FILE)
            opSaveIncompatibilityGraph(source, sink, avect, bvect);
#endif
        }   // while (!done_augmenting_path)

        // Identify C1, C2, D1, and D2
        // C1 and D2 compose the min weight vertex cover
        // C2 and D` compose the independent set
        for (auto & source_edge : source._edges) {
            OPVertex * avertex = source_edge->_to;
            if (source_edge->_capacity > 0.0 || avertex->_residual_capacity > 0.0) {
                // Because this source edge has remaining capacity, its distal vertex is part of the independent set
                C2.insert(*(avertex->_split));

                // This A vertex allows access to the B side, so any connected B vertices with
                // zero capacity are part of the vertex cover
                for (const auto central_edge : avertex->_edges) {
                    OPVertex * bvertex = central_edge->_to;
                    assert(bvertex->_edges.size() == 1);
                    OPEdge * sink_edge = bvertex->_edges[0];
                    if (sink_edge->_capacity == 0.0) {
                        D2.insert(*(bvertex->_split));
                    }
                }
            }
            else {
                // Because this A-vertex has zero capacity, it is part of the vertex cover
                C1.insert(*(avertex->_split));

                // No need to consider connected B vertices because this A-vertex already
                // covers all connected edges
            }
        }

        // D1 includes every split not in D2
        D1.clear();
        for (auto & b : bvect) {
            if (D2.count(*(b._split)) == 0) {
                D1.insert(*(b._split));
            }
        }
    }

    inline bool OP::opRefineSupport(const Split::treeid_pair_t & AB, Split::treeid_pair_t & AB1, Split::treeid_pair_t & AB2) const {
        // Create a vector of incompatibility graph vertices
        vector<OPVertex> avect(AB.first.size());
        vector<OPVertex> bvect(AB.second.size());

        // Calculate weights for the "A" vertices.
        // For example, if A = {a1,a2,a3}, then the weight of a1 is
        // weight a1 = a1^2 / (a1^2 + a2^2 + a3^2)
        unsigned aindex = 0;
        double asum = 0.0;
        for (auto & a : AB.first) {
            asum += pow(a.getEdgeLen(),2);
        }
        for (auto & a : AB.first) {
            avect[aindex]._split = &a;
            avect[aindex]._weight = pow(a.getEdgeLen(),2)/asum;
            aindex++;
        }

        // Calculate weights for the "B" vertices.
        // For example, if B = {b1,b2,b3}, then the weight of b1 is
        // weight b1 = b1^2 / (b1^2 + b2^2 + b3^2)
        unsigned bindex = 0;
        double bsum = 0.0;
        for (auto & b : AB.second) {
            bsum += pow(b.getEdgeLen(),2);
        }
        for (auto & b : AB.second) {
            bvect[bindex]._split = &b;
            bvect[bindex]._weight = pow(b.getEdgeLen(),2)/bsum;
            bindex++;
        }

        // Create the incompatibility graph
        // A vertices go on the left, B vertices go on the right, and edges connect an A vertex
        // to a B vertex only if the two vertices are incompatible.
        vector<OPEdge> all_edges;
        all_edges.reserve(avect.size() * bvect.size() + avect.size() + bvect.size());
        unsigned nincompatibilities = 0;
        auto asize = static_cast<unsigned>(avect.size());
        auto bsize = static_cast<unsigned>(bvect.size());

        // Create edges from source to the A-vertices
        OPVertex source;
        source._name = "source";
        for (unsigned i = 0; i < asize; i++) {
            // Assign a name to avect[i]
            avect[i]._name = str(format("a%d") % i);

            // Create the forward edge
            all_edges.emplace_back();
            OPEdge & source_forward_edge = all_edges.back();
            source_forward_edge._from = &source;
            source_forward_edge._to = &avect[i];
            source_forward_edge._capacity = avect[i]._weight;
            source_forward_edge._flow = 0.0;
            source_forward_edge._reverse_flow = 0.0;
            source_forward_edge._open = true;
            source._edges.push_back(&source_forward_edge);
        }

        // Create edges from A-vertices to B-vertices
        for (unsigned i = 0; i < asize; i++) {
            // Get split associated with avect[i]
            const Split * a = avect[i]._split;
            assert(a);
            for (unsigned j = 0; j < bsize; j++) {
                // Get split associated with bvect[j]
                const Split * b = bvect[j]._split;
                assert(b);

                if (!a->compatibleWith(*b)) {
                    nincompatibilities++;

                    // Create a forward edge in the incompatibility graph
                    all_edges.emplace_back();
                    OPEdge & forward_edge = all_edges.back();
                    forward_edge._from = &avect[i];
                    forward_edge._to = &bvect[j];
                    forward_edge._capacity = 1.0;
                    forward_edge._flow = 0.0;
                    forward_edge._reverse_flow = 0.0;
                    forward_edge._open = true;
                    avect[i]._edges.emplace_back(&forward_edge);
                }
            }
        }

        // Create edges from the B-vertices to the sink
        OPVertex sink;
        sink._name = "sink";
        for (unsigned j = 0; j < bsize; j++) {
            // Assign a name to bvect[j]
            bvect[j]._name = str(format("b%d") % j);

            // Create the edge
            all_edges.emplace_back();
            OPEdge & sink_forward_edge = all_edges.back();
            sink_forward_edge._from = &bvect[j];
            sink_forward_edge._to = &sink;
            sink_forward_edge._capacity = bvect[j]._weight;
            sink_forward_edge._flow = 0.0;
            sink_forward_edge._reverse_flow = 0.0;
            sink_forward_edge._open = true;
            bvect[j]._edges.emplace_back(&sink_forward_edge);
        }

        bool success = false;
        if (nincompatibilities < asize*bsize) {
            // At least one independent pair of vertices exists
            // Carry out Edmonds-Karp algorithm to find min-weight vertex cover (identifies max weight independent set)
            // In Owens-Provan terminology,
            //   C1 = A's contribution to vertex cover     (equals AB1.first)
            //   C2 = A's contribution to independent set  (equals AB2.first)
            //   D1 = B's contribution to independent set  (equals AB1.second)
            //   D2 = B's contribution to vertex cover     (equals AB2.second)
            // A1 and A2 must represent a non-trivial partition of A
            // B1 and B2 must represent a non-trivial partition of B
            // C2 and D1 (AB2.first and AB1.second) are compatible sets of splits
            // C1 and D2 (AB1.first and AB2.second) compose the minimum weight vertex cover
            // ||C1||/||D1|| < ||C2||/||D2|| must be true
            // If all of the above conditions hold, then the support can be refined:
            // A = (A1,A2) and B = (B1,B2)
            // where A1 = C1, A2 = C2, B1 = D1, and B2 = D2
            // Length of this segment of the geodesic is
            //   L = sqrt{ (||A1|| + ||B1||)^2 +  (||A2|| + ||B2||)^2 }
            // Orthants crossed:
            //   start:        A1, A2
            //      edges in A1 decreasing, edges in B1 increasing
            //   intermediate: B1, A2
            //      edges in A2 decreasing, edges in B2 increasing
            //   finish:       B1, B2
            // The point at which A1 edges become 0 is
            //   ||A1||/(||A1|| + ||B1||)
            // The point at which A2 edges become 0 is
            //   ||A2||/(||A2|| + ||B2||)
            opEdmondsKarp(source, sink, avect, bvect, AB1.first, AB2.first, AB1.second, AB2.second, _noisy);

            // if (!quiet) {
            //     cout << "\nResults:" << endl;
            //     cout << str(format("  ||C1|| = %.9f") % C1len) << endl;
            //     cout << str(format("  ||C2|| = %.9f") % C2len) << endl;
            //     cout << str(format("  ||D1|| = %.9f") % D1len) << endl;
            //     cout << str(format("  ||D2|| = %.9f") % D2len) << endl;
            //     cout << "  Check whether ||C1||/||D1|| < ||C2||/||D2||" << endl;
            //     cout << str(format("    %.9f < %.9f") % (C1len/D1len) % (C2len/D2len));
            //     cout << endl;
            // }

            // Check conditions for successful refinement
            bool A_is_trivial = (AB1.first.empty()) || (AB2.first.empty());
            bool B_is_trivial = (AB1.second.empty())|| (AB2.second.empty());
            double C1len = opCalcTreeIDLength(AB1.first);
            double C2len = opCalcTreeIDLength(AB2.first);
            double D1len = opCalcTreeIDLength(AB1.second);
            double D2len = opCalcTreeIDLength(AB2.second);
            //double fraction1 = C1len/(C1len + D1len);
            //double fraction2 = C2len/(C2len + D2len);
            bool P2 = C1len/D1len < C2len/D2len;
            success = !A_is_trivial && !B_is_trivial && P2;

            if (_noisy) {
                if (success) {
                    cout << "\nSuccessfully refined support:" << endl;
                }
                else {
                    cout << "\nSupport not refined because:" << endl;
                    if (A_is_trivial) {
                        cout << "  A is a trivial partition:" << endl;
                    }
                    if (B_is_trivial) {
                        cout << "  B is a trivial partition:" << endl;
                    }
                    if (!P2) {
                        cout << "  ||C1||/||D1|| is not strictly less than ||C2||/||D2||" << endl;
                    }
                }
                // cout << "  Input A vertices:" << endl;
                // for (auto & a : AB.first) {
                //     cout << "    " << a.createPatternRepresentation() << endl;
                // }
                // cout << "  Input B vertices:" << endl;
                // for (auto & b : AB.second) {
                //     cout << "    " << b.createPatternRepresentation() << endl;
                // }
                cout << "  A1 vertices:" << endl;
                if (AB1.first.empty()) {
                    cout << "    empty set" << endl;
                }
                else {
                    for (auto & a : AB1.first) {
                        cout << "    " << a.createPatternRepresentation(true) << endl;
                    }
                }
                cout << "  A2 vertices:" << endl;
                if (AB2.first.empty()) {
                    cout << "    empty set" << endl;
                }
                else {
                    for (auto & a : AB2.first) {
                        cout << "    " << a.createPatternRepresentation(true) << endl;
                    }
                }
                cout << "  B1 vertices:" << endl;
                if (AB1.second.empty()) {
                    cout << "    empty set" << endl;
                }
                else {
                    for (auto & b : AB1.second) {
                        cout << "    " << b.createPatternRepresentation(true) << endl;
                    }
                }
                cout << "  B2 vertices:" << endl;
                if (AB2.second.empty()) {
                    cout << "    empty set" << endl;
                }
                else {
                    for (auto & b : AB2.second) {
                        cout << "    " << b.createPatternRepresentation(true) << endl;
                    }
                }
                cout << endl;
            }
        }
        return success;
    }

    inline double OP::opCalcGeodesicDist(const Split::treeid_pair_t & inpair) const {
        // Assumes a_splits and b_splits have no common edges
        vector<Split::treeid_pair_t> ABpairs;
        ABpairs.push_back(inpair);
        vector<Split::treeid_pair_t> support;
        bool done = false;
        while (!done) {
            unsigned nrefinements = 0;
            for (auto & ABpair : ABpairs) {
                Split::treeid_pair_t AB1;
                Split::treeid_pair_t AB2;
                bool success = opRefineSupport(ABpair, AB1, AB2);
                if (success) {
                    // ABpair was successfully refined, so add AB1 and AB2 to support
                    support.push_back(AB1);
                    support.push_back(AB2);
                    nrefinements++;
                }
                else {
                    // ABpair was not successfully refined, so add ABpair to support
                    support.push_back(ABpair);
                }
            }
            done = (nrefinements == 0);
            if (!done) {
                ABpairs = support;
                support.clear();
            }
        }

        // Calculate geodesic distance
        unsigned ratio_index = 1;
        double geodesic_distance = 0.0;
        for (auto & AB : support) {
            double dropped_length = opCalcTreeIDLength(AB.first);
            double added_length   = opCalcTreeIDLength(AB.second);
            double fraction = dropped_length/(dropped_length + added_length);
            double ratio = dropped_length/added_length;
            geodesic_distance += pow(dropped_length + added_length, 2);

            if (_noisy) {
                cout << str(format("\nRatio %d: %.9f (lambda = %.9f)") % ratio_index % ratio % fraction) << endl;
                cout << "  Edges dropped:" << endl;
                for (auto & a : AB.first) {
                    cout << "    " << a.createPatternRepresentation() << endl;
                }
                cout << "  Edges added:" << endl;
                for (auto & b : AB.second) {
                    cout << "    " << b.createPatternRepresentation() << endl;
                }
            }

            ++ratio_index;
        }

        geodesic_distance = sqrt(geodesic_distance);
        return geodesic_distance;
    }

    inline void OP::buildRefTree(TreeManip & tm) const {
        string newick = _tree_summary->getRefNewick();
        assert(!newick.empty());

        bool isrooted = _tree_summary->isRefRooted();
        tm.setIsRooted(isrooted);
        try {
            tm.buildFromNewick(newick, /*rooted*/isrooted, /*allow_polytomies*/true);
        } catch (Xop & e) {
            cerr << "Could not build reference tree" << endl;
            throw e;
        }
    }

    inline void OP::buildTree(unsigned tree_index, TreeManip & tm) const {
        string newick = _tree_summary->getNewick(tree_index);

        bool isrooted = _tree_summary->isRooted(tree_index);
        tm.setIsRooted(isrooted);
        //if (!isrooted) {
        //    throw Xop("Trees must be rooted in this version");
        //}
        try {
            tm.buildFromNewick(newick, /*rooted*/isrooted, /*allow_polytomies*/true);
        } catch (Xop & e) {
            cerr << "Could not build tree with index " << tree_index << endl;
            throw e;
        }
        //tm.setLeafNames(_taxon_labels);
    }

    inline double OP::calcBHVDistance(
            TreeManip & starttm,
            TreeManip & endtm,
            vector<pair<Split::treeid_t, Split::treeid_t> > & in_pairs,
            vector<Split::treeid_pair_t> & ABpairs,
            vector<Split::split_pair_t> & commonPairs) const {
        if (_noisy) {
            cout << "\nKey to taxa:" << endl;
            for (unsigned i = 0; i < TreeManip::_taxon_names.size(); i++) {
                cout << str(format("%12d %s") % (i+1) % TreeManip::_taxon_names[i]) << endl;
            }
            cout << "\nStarting tree description:" << endl;
            cout << "  " << starttm.makeNewick(9,false) << endl;
            cout << "\nEnding tree description:" << endl;
            cout << "  " << endtm.makeNewick(9,false) << endl;
        }

        ABpairs.clear();

        // Store splits from the starting tree
        Split::treeid_t A0;
        Split::treeid_t Alvs;
        starttm.storeSplits(A0, Alvs);

        // Only save splits with non-zero edge length
        Split::treeid_t A;
        for (auto & a : A0) {
            if (a.getEdgeLen() > 0.0) {
                A.insert(a);
            }
            else {
                // Drop this split because it has an edge length of zero
                starttm.dropSplit(a); //TODO: is this OK?
            }
        }

        if (_noisy) {
            cout << "\nInternal splits from starting tree:" << endl;
            for (const auto& a : A) {
                cout << "  " << a.createPatternRepresentation(true) << endl;
            }
        }

        // Store splits from the ending tree
        Split::treeid_t B0;
        Split::treeid_t Blvs;
        endtm.storeSplits(B0, Blvs);

        // Only save splits with non-zero edge length
        Split::treeid_t B;
        for (auto & b : B0) {
            if (b.getEdgeLen() > 0.0) {
                B.insert(b);
            }
            else {
                // Drop this split because it has an edge length of zero
                endtm.dropSplit(b); //TODO: is this OK?
            }
        }

        if (_noisy) {
            cout << "\nInternal splits from ending tree:" << endl;
            for (const auto& b : B) {
                cout << "  " << b.createPatternRepresentation(true) << endl;
            }
        }

        // Find common edges and calculate the contribution of common edge lengths to the geodesic
        double common_edge_contribution_squared = opFindCommonEdges(A, B, commonPairs);

        in_pairs.emplace_back(A,B);
        opSplitAtCommonEdges(commonPairs, in_pairs);

        // Calculate the contribution of leaf edges to the geodesic
        double leaf_contribution_squared = opCalcLeafContribution(Alvs, Blvs, commonPairs);

        unsigned pair_index = 1;
        vector<double> geodesic_distances;
        for (const auto & inpair : in_pairs) {
            if (_noisy)
                cout << str(format("\nTree pair %d (of %d)") % pair_index % in_pairs.size()) << endl;

            ABpairs.push_back(inpair);

            if (_noisy) {
                cout << "  A splits:" << endl;
                for (const auto& a : inpair.first) {
                    cout << "    " << a.createPatternRepresentation(true) << endl;
                }

                cout << "  B splits:" << endl;
                for (const auto& b : inpair.second) {
                    cout << "    " << b.createPatternRepresentation(true) << endl;
                }
            }

            double L = opCalcGeodesicDist(inpair);

            if (_noisy)
                cout << str(format("  L for tree pair %d = %.9f") % pair_index % L) << endl;

            geodesic_distances.push_back(L);
            ++pair_index;
        }

        if (_noisy)
            cout << endl;

        // Calculate total geodesic distance
        double total_geodesic_distance = 0.0;
        for (double x : geodesic_distances) {
            total_geodesic_distance += pow(x, 2);
        }
        total_geodesic_distance += leaf_contribution_squared;
        total_geodesic_distance += common_edge_contribution_squared;
        total_geodesic_distance = sqrt(total_geodesic_distance);

        if (_noisy)
            cout << str(format("Total geodesic distance = %.9f") % total_geodesic_distance) << endl;

        return total_geodesic_distance;
    }

#if defined(CLUSTER_DISTANCE)
    inline double OP::calcClusterDistance(
        TreeManip & starttm,
        TreeManip & endtm,
        const vector<pair<Split::treeid_t, Split::treeid_t> > & in_pairs,
        const vector<Split::split_pair_t> & commonPairs) const {
        // commonPairs are the splits common to both input trees (including trivial splits)
        // in_pairs.first holds the unique splits in the first input tree
        // in_pairs.second holds the unique splits in the second input tree
        assert(in_pairs.size() == 1);
        if (_noisy) {
            cerr << "\n********** calcClusterDistance **********" << endl;
        }

        // Calculate mci_score of common splits
        double mci_score = 0.0;
        double total_entropy = 0.0;
        for (const auto & s : commonPairs) {
            unsigned nbits = s.first.getNumBitsSet();
            unsigned nlvs = s.first.getSize();
            if (nbits > 1 && nbits < nlvs - 1) {
                double h = 2.0*s.first.entropy();
                total_entropy += h;
                double I = s.first.mutualClusteringInfo(s.first);
                mci_score += I;
                if (_noisy) {
                    cerr << "  common pair:  " << s.first.createPatternRepresentation() << " (h = " << h << ", I = " << I << ", mci_score = " << mci_score << ")" << endl;
                }
            }
        }

        // Create a cost matrix with A splits as rows and B splits as columns
        auto Asplits = in_pairs[0].first;
        auto Bsplits = in_pairs[0].second;

        auto nAsplits = static_cast<unsigned>(Asplits.size());
        auto nBsplits = static_cast<unsigned>(Bsplits.size());
        assert(nAsplits == nBsplits);
        LAPJV lapjv(nAsplits);
        int i = 0;
        int j = 0;
        double maxuint = pow(2,22);
        vector<double> entropyAsplits(nAsplits);
        vector<double> entropyBsplits(nBsplits);
        for (const auto & b : Bsplits) {
            entropyBsplits[j++] = b.entropy();
        }
        for (const auto & a : Asplits) {
            entropyAsplits[i] = a.entropy();
            j = 0;
            for (const auto & b : Bsplits) {
                double mci = a.mutualClusteringInfo(b);
                double mci_truncated = floor(mci*maxuint);
                lapjv.assignCost(i, j, -mci_truncated);
                ++j;
            }
            ++i;
        }
        double total_cost = lapjv.lap()/maxuint;
        mci_score -= total_cost;

        vector<unsigned> optimal_pairings;
        lapjv.getOptimalPairings(optimal_pairings);
        for (unsigned ii = 0; ii < optimal_pairings.size(); ++ii) {
            unsigned jj = optimal_pairings[ii];
            total_entropy += entropyAsplits[ii];
            total_entropy += entropyBsplits[jj];
        }

        double max_value = total_entropy/2.0;
        double d = (max_value - mci_score)/max_value;
        if (_noisy) {
            cerr << "max_value = " << max_value << endl;
            cerr << "mci_score = " << mci_score << endl;
            cerr << "distance = " << d << endl;
        }
        return d;
    }
#endif

#if defined(KUHNER_FELSENSTEIN_DISTANCE)
    inline double OP::calcKFDistance(unsigned ref_index, unsigned test_index) const {
        // Get the reference tree
        string ref_newick = _tree_summary->getNewick(ref_index);
        bool ref_isrooted = _tree_summary->isRooted(ref_index);

        // Get the test tree
        string test_newick = _tree_summary->getNewick(test_index);
        bool test_isrooted = _tree_summary->isRooted(ref_index);

        // Ensure both trees are rooted
        if (!ref_isrooted || !test_isrooted) {
            throw Xop(format("Trees must be rooted in this version of %s") % OP::_program_name);
        }

        // Build the reference tree
        TreeManip reftm;
        reftm.buildFromNewick(ref_newick, /*rooted*/ref_isrooted, /*allow_polytomies*/true);
        //TODO: get rooted status from treeManip object

        // Store splits from the reference tree
        Split::treeid_t refsplits;
        Split::treeid_t reflvs;
        reftm.storeSplits(refsplits, reflvs);

        // Build the test tree
        TreeManip testtm;
        testtm.buildFromNewick(test_newick, /*rooted*/test_isrooted, /*allow_polytomies*/false);

        // Store splits from the reference tree
        Split::treeid_t testsplits;
        Split::treeid_t testlvs;
        testtm.storeSplits(testsplits, testlvs);

        // Store union of refsplits and testsplits in allsplits
        Split::treeid_t allsplits;
        set_union(
            refsplits.begin(), refsplits.end(),
            testsplits.begin(), testsplits.end(),
            inserter(allsplits, allsplits.begin()));

        // Traverse allsplits, storing squared branch length differences in KLinternals
        vector<double> KLinternals(allsplits.size());
        unsigned i = 0;
        for (auto s : allsplits) {
            Node * nd0 = reftm.getNodeWithSplit(s);
            Node * nd  = testtm.getNodeWithSplit(s);
            assert(!(nd0 == nullptr && nd == nullptr));
            if (nd0 == nullptr) {
                double edge_length = nd->getEdgeLength();
                double square = pow(edge_length, 2.0);
                KLinternals[i++] = square;
            }
            else if (nd == nullptr) {
                double edge_length = nd0->getEdgeLength();
                double square = pow(edge_length, 2.0);
                KLinternals[i++] = square;
            }
            else {
                double edge_length0 = nd0->getEdgeLength();
                double edge_length  = nd->getEdgeLength();
                double square = pow(edge_length0 - edge_length, 2.0);
                KLinternals[i++] = square;
            }
        }

        // Create the map in which keys are taxon names and values are Node pointers
        // for the reference tree
        map<string, Node *> leafmap0;
        reftm.createLeafNodeMap(leafmap0);

        // Create a map in which keys are taxon names and values are Node pointers
        // for the test tree
        map<string, Node *> leafmap;
        testtm.createLeafNodeMap(leafmap);

        // The two trees should have the same number of leaves
        assert(leafmap0.size() == leafmap.size());

        // Get taxon names from the reference tree (assuming the taxon names
        // in the test tree are the same)
        vector<string> names;
        names.reserve(leafmap0.size());
    for (const auto& p : leafmap0) {
            names.push_back(p.first);
        }
        sort(names.begin(), names.end());

        // Now calculate squares for leaf nodes, storing in KLleaves
        vector<double> KLleaves(names.size());
        i = 0;
        for (const auto& nm : names) {
            Node * nd0 = leafmap0[nm];
            Node * nd  = leafmap[nm];
            double edge_length0 = nd0->getEdgeLength();
            double edge_length  = nd->getEdgeLength();
            double square = pow(edge_length0 - edge_length, 2.0);
            KLleaves[i++] = square;
        }

        // Calculate KL distance
        double KLdist = 0.0;
        for (auto square : KLinternals) {
            KLdist += square;
        }
        for (auto square : KLleaves) {
            KLdist += square;
        }

        return KLdist;
    }
#endif

    inline unsigned OP::chooseRandomTree(TreeManip & tm, Lot & lot) const {
        int n = static_cast<int>(_tree_summary->getNumTrees());
        auto index = static_cast<unsigned>(lot.randint(0, n-1));
        string newick = _tree_summary->getNewick(index);
        bool rooted = _tree_summary->isRooted(index);
        tm.buildFromNewick(newick, rooted, /*allow_polytomies*/true);
        return index;
    }

    inline void OP::displaceTreeAlongGeodesic(TreeManip & start_tree, TreeManip & end_tree, double displacement) const {
        // Move start_tree a distance displacement along the geodesic from start_tree to end_tree
        // First, get geodesic
        vector<Split::treeid_pair_t> ABpairs;
        vector<Split::split_pair_t> commonPairs;
        vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
        calcBHVDistance(start_tree, end_tree, in_pairs, ABpairs, commonPairs);

        // Support A = (A1, A2, ..., Ak) and B = (B1, B2, ..., Bk)
        // Path has i-1,2,...,k "legs"
        //    G0 if lambda/(1-lambda) <= length(A1)/length(B1)
        //    Gi if length(Ai)/length(Bi) = lambda/(1-lambda) <= length(Ai+1)/length(Bi+1)
        //    Gk if length(Ak)/length(Bk) <= lambda/(1-lambda)
        // and tree Ti along path has edge sets
        //    B1 U ... U Bi U Ai+1 U ... U Ak
        // and edge sets
        //    length edge e in Ti = [(1-lambda) length(Aj) - lambda length(Bj)]/length(Aj) if e in Aj
        //    length edge e in Ti = [lambda length(Bj) - (1-lambda) length(Aj)]/length(Bj) if e in Bj

        // Precalculate lengths of all segments
        auto support_size = static_cast<unsigned>(ABpairs.size());
        vector<double> lenA(support_size);
        vector<double> lenB(support_size);
        for (unsigned i = 0; i < support_size; ++i) {
            lenA[i] = opCalcTreeIDLength(ABpairs[i].first);
            lenB[i] = opCalcTreeIDLength(ABpairs[i].second);
        }

        // Walk down ABpairs, dropping and adding splits as needed from the start_tree until we arrive at the destination leg
        double lambda = displacement;
        unsigned leg = 0;
        while (true) {
            double lambda_leg = 1.0;
            if (leg < support_size) {
                double ratio_leg = lenA[leg]/lenB[leg];
                lambda_leg = ratio_leg/(1.0 + ratio_leg);
            }
            if (lambda <= lambda_leg) {
                break;
            }

            for (auto & asplit : ABpairs[leg].first) {
                // Drop asplit from the start_tree
                start_tree.dropSplit(asplit);
            }

            for (auto & bsplit : ABpairs[leg].second) {
                // Add bsplit to the start_tree
                start_tree.addSplit(bsplit);
            }

            ++leg;
        }

        // Modify edge lengths
        for (unsigned i = 0; i < leg; ++i) {
            for (auto & bsplit : ABpairs[i].second) {
                double edgelen_multiplicative_factor = (lambda*lenB[i] - (1.0 - lambda)*lenA[i])/lenB[i];
                double edgelen = bsplit.getEdgeLen();
                start_tree.setEdgeLength(bsplit, edgelen_multiplicative_factor*edgelen);
            }
        }
        for (unsigned i = leg; i < ABpairs.size(); ++i) {
            for (auto & asplit : ABpairs[i].first) {
                double edgelen_multiplicative_factor = ((1 - lambda)*lenA[i] - lambda*lenB[i])/lenA[i];
                double edgelen = asplit.getEdgeLen();
                start_tree.setEdgeLength(asplit, edgelen_multiplicative_factor*edgelen);
            }
        }

        for (auto & commonPair : commonPairs) {
            double from_edgelen = commonPair.first.getEdgeLen();
            double to_edgelen = commonPair.second.getEdgeLen();
            double new_edgelen = from_edgelen + lambda*(to_edgelen - from_edgelen);
            start_tree.setEdgeLength(commonPair.first, new_edgelen);
        }
    }

    inline unsigned OP::buildThreadSchedule() {
        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads > 0) {
            std::cout << "Number of hardware threads available: " << num_threads << std::endl;
        }

        // Let n be the number of pairwise comparisons
        unsigned n = _frechet_n*(_frechet_n - 1)/2;
        assert(n > 0);
        _thread_sched.clear();
        // Example: frechet-n = 10 -> n = 10*9/2=45 pairwise comparisons, _nthreads = 4
        // entities_per_thread = 45/4=11, remainder = 45-4*11=1
        // thread 0: begin = 0, end = 0+11+1=12, remainder = 1-1=0
        // thread 1: begin = 12, end = 12+11+0=23, remainder = 0
        // thread 2: begin = 23, end = 23+11+0=34, remainder = 0
        // thread 3: begin = 34, end = 34+11+0=45, remainder = 0
        unsigned entities_per_thread = (unsigned)floor(1.0*n/_nthreads);
        unsigned remainder = n - _nthreads*entities_per_thread;
        unsigned begin = 0;
        for (unsigned i = 0; i < _nthreads; i++) {
            unsigned end = begin + entities_per_thread;
            if (remainder > 0) {
                // Add another job to early threads to use up remainder
                end++;
                remainder--;
            }
            _thread_sched.push_back(make_pair(begin,end));
            begin = end;
        }
        return n;
    }

    inline void OP::calcBHVDistanceForPairs(
            unsigned start_pair,
            unsigned end_pair,
            vector<TreeManip> & mu,
            vector<double> & bhv_distances,
            vector<pair<Split::treeid_t, Split::treeid_t> > & in_pairs,
            vector<Split::treeid_pair_t> & ABpairs,
            vector<Split::split_pair_t> & commonPairs) const {
        // Find the mean trees in the range of pairs starting at start_pair and ending at end_pair
        for (unsigned k = start_pair; k < end_pair; ++k) {
            // Find indices that go with pair k
            unsigned i = static_cast<unsigned>(floor(0.5*(1.0 + sqrt(1.0 + 8.0*k))));
            unsigned j = k - i*(i - 1.0)/2.0;
            TreeManip tm0(mu[i].getTree());
            TreeManip tm1(mu[j].getTree());
            double bhvdist = calcBHVDistance(tm0, tm1, in_pairs, ABpairs, commonPairs);
            {
                lock_guard<mutex> guard(_mutex);
                assert(k < bhv_distances.size());
                bhv_distances[k] = bhvdist;
            }
        }
    }

    inline bool OP::frechetCloseEnough(vector<TreeManip> & mu, unsigned lower, unsigned upper, double epsilon) const {
        // Compute pairwise distances between trees in mu with index >= lower and index < upper and return
        // true iff all pairwise distances are less than epsilon
        bool is_close_enough = true;
        vector<Split::treeid_pair_t> ABpairs;
        vector<Split::split_pair_t> commonPairs;
        vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
        assert(lower < upper);
        assert (upper <= mu.size());
        if (_nthreads > 1) {
            // Assign each pair to a different thread until all are done
            unsigned npairs = _frechet_n*(_frechet_n - 1)/2;
            vector<thread> threads;
            vector<double> bhv_distances(npairs, 0.0);
            for (unsigned i = 0; i < _nthreads; i++) {
                unsigned start_pair = _thread_sched[i].first;
                unsigned end_pair = _thread_sched[i].second;
                threads.emplace_back(thread(&OP::calcBHVDistanceForPairs,
                    this,
                    start_pair,
                    end_pair,
                    std::ref(mu),
                    std::ref(bhv_distances),
                    std::ref(in_pairs),
                    std::ref(ABpairs),
                    std::ref(commonPairs))
                );
            }

            // The join function causes this loop to pause until the ith thread finishes
            for (unsigned i = 0; i < threads.size(); i++) {
                threads[i].join();
            }

            // Find largest bhv distance
            double largest = *max_element(bhv_distances.begin(), bhv_distances.end());
            if (largest > epsilon) {
                is_close_enough = false;
            }
        }
        else {
            // first pair tried is mu[lower] vs. mu[upper-1], which should be the most distant pair
            // increasing the chances of breaking the loop early if possible
            for (unsigned i = lower; i < upper - 1; ++i) {
                for (unsigned j = upper - 1; j > i; --j) {
                    in_pairs.clear();
                    double bhvdist = calcBHVDistance(mu[i-1], mu[j-1], in_pairs, ABpairs, commonPairs);
                    if (bhvdist > epsilon) {
                        is_close_enough = false;
                        break;
                    }
                }
                if (!is_close_enough)
                    break;
            }
        }
        return is_close_enough;
    }

    inline unsigned OP::computeFrechetMean(TreeManip & mean_tree) const {
        // Returns the number of iterations required to compute the Frechet mean tree
        double   epsilon = _frechet_epsilon; // successive mean estimates must be at least this close to stop iterating
        unsigned N = _frechet_n;   // number of previous mean estimates that must be as close as epsilon
        unsigned K = _frechet_k; // maximum number of iterations
        unsigned k = 1;   // keeps track of iterations
        vector<TreeManip> mu;
        //mu.reserve(K);// the trail of estimated mean trees (always has length k)
        Lot lot;
        lot.setSeed(_random_number_seed);
        mu.emplace_back();
        //@ unsigned prev_index = chooseRandomTree(mu[k-1], lot);
        chooseRandomTree(mu[k-1], lot);
        bool done = false;
        //@ unsigned curr_index = 0;
        while (!done) {
            ++k;
            mu.emplace_back();
            assert(mu.size() == k);

            //@ curr_index = chooseRandomTree(mu[k-1], lot);
            chooseRandomTree(mu[k-1], lot);

            displaceTreeAlongGeodesic(mu[k-1], mu[k-2], 1.0*k/(k+1));
            //@ prev_index = curr_index;
            if (k >= K) {
                done = true;
            }
            if (k > N) {
                done = frechetCloseEnough(mu, k-N, k, epsilon);
            }
        }
        mean_tree.setTree(mu[k-1].getTree());
        return k;
    }

#if defined(TESTKDE)
    inline void OP::testKDE() {
        // Initialize pseudorandom number generator
        Lot lot;
        lot.setSeed(12345);

        // Simulate Gamma(shape, scale) data
        double shape = 2.0;
        double scale = 3.0;
        unsigned sample_size = 1000;
        vector<double> sample(sample_size);
        for (unsigned i = 0; i < sample_size; ++i) {
            sample[i] = lot.gamma(shape,scale);
        }

        // Calculate _kde_sigma and _kde_band_width
        sort(sample.begin(), sample.end());
        calcBandWidth(sample);

        // Save x, true density of x, KDE of x to a file suitable for use with both R and Tracer
        double min_xvalue = 0.0;
        double max_xvalue = 10.0;
        unsigned n = 10000;
        double incr = (max_xvalue - min_xvalue)/n;
        vector<string> x(n);
        vector<string> ptrue(n);
        vector<string> pkde(n);
        double xvalue = 0.0;
        for (unsigned i = 0; i < n; ++i) {
            double kde = kernelDensity(xvalue, sample);
            double true_density = 0.0;
            if (xvalue > 0.0) {
                double logf = (shape - 1.0)*log(xvalue) - xvalue/scale - shape*log(scale) - boost::math::lgamma(shape);
                true_density = exp(logf);
            }
            x[i] = str(format("%0.9f") % xvalue);
            ptrue[i] = str(format("%0.9f") % true_density);
            pkde[i] = str(format("%0.9f") % kde);
            xvalue += incr;
        }

        vector<string> xsampled(sample_size);
        for (unsigned i = 0; i < sample_size; ++i) {
            xsampled[i] = str(format("%0.9f") % sample[i]);
        }

        ofstream outf("kde_test.R");
        outf << "x <- c(" << boost::join(x, ", ") << ")\n";
        outf << "ptrue <- c(" << boost::join(ptrue, ", ") << ")\n";
        outf << "pkde <- c(" << boost::join(pkde, ", ") << ")\n";
        outf << "xsampled <- c(" << boost::join(xsampled, ", ") << ")\n";
        outf << "plot(x, ptrue, type=\"l\", lty=\"solid\", lwd=1, col=\"navy\", xlab=\"x\", ylab=\"density\", main=\"KDE test\")\n";
        outf << "lines(x, pkde, lty=\"solid\", lwd=5, col=\"red\")\n";
        outf << "lines(density(xsampled), lty=\"solid\", lwd=3, col=\"green\")\n";
        outf << "quantile(xsampled, probs = c(0.25, 0.50, 0.75), type = 7)\n";
        outf.close();

        cerr << "File \"kde_test.R\" has been saved." << endl;
    }
#endif

    inline void OP::calcBandWidth(const vector<double> & sample) {
        // Assumes sample is already sorted from lowest to highest value
        // Calculate IQR (InterQuartile Range)
        auto n = static_cast<unsigned>(sample.size());

        // Calculate the 0.25 quantile
        double p = 0.25;
        double k = p*(n - 1.0) + 1.0;
        auto ki = static_cast<unsigned>(floor(k));
        _kde_q25 = sample[ki - 1];
        if (k != ki) {
            // k is not an integer, so interpolate
            unsigned kiplus1 = ki + 1;
            _kde_q25 = (sample[ki-1] + (k - ki)*(sample[kiplus1-1] - sample[ki-1]));
        }

        // Calculate the 0.75 quantile
        p = 0.75;
        k = p*(n - 1.0) + 1.0;
        ki = static_cast<unsigned>(floor(k));
        _kde_q75 = sample[ki - 1];
        if (k != ki) {
            // k is not an integer, so interpolate
            unsigned kiplus1 = ki + 1;
            _kde_q75 = (sample[ki-1] + (k - ki)*(sample[kiplus1-1] - sample[ki-1]));
        }
        double IQR = _kde_q75 - _kde_q25;

        // Calculate the sample standard deviation
        double sumx = 0.0;
        double sumxsq = 0.0;
        for (double xi : sample) {
            sumx += xi;
            sumxsq += pow(xi, 2.0);
        }
        double meanx = sumx/n;
        double varx = (sumxsq - n*pow(meanx, 2.0))/(n-1);
        _kde_sigma = sqrt(varx);

        // Calculate the rule-of-thumb bandwidth
        double A = min(_kde_sigma, IQR/1.34);
        _kde_bandwidth = 0.9*A*pow(n, -0.2);
    }

    inline double OP::kernelDensity(double x, const vector<double> & sample) const {
        if (x <= 0.0)
            return 0.0;
        auto n = static_cast<double>(sample.size());
        double sigma = 1.0; //_kde_sigma;
        double numer = 0.0;
        for (double xi : sample) {
            double term = (x - xi)/(_kde_bandwidth*sigma);
            numer += exp(-0.5*pow( term, 2 ));
        }
        double denom = n*_kde_bandwidth*sigma*sqrt(2.0*M_PI);
        double f = numer/denom;
        return f;
    }

    inline void OP::rPlotDists(vector<double> & dists, string & rscript, double & radius, double & hpd_lower, double & hpd_upper, unsigned hpd_level) {
        assert(!dists.empty());

        // Create a data set mirrored across the lower boundary 0.0
        auto n = static_cast<unsigned>(dists.size());
        vector<double> mirrored(2*n);
        for (unsigned i = 0; i < n; ++i) {
            mirrored[i] = dists[i];
            mirrored[i+n] = -1.0*dists[i];
        }

        // Calculate the bandwidth for kernel density estimation
        sort(mirrored.begin(), mirrored.end());
        calcBandWidth(mirrored);

        // Create a vector containing (density, dist) tuples for every element of dists
        vector<pair<double, double> > density_dist;
        for (double x : dists) {
            double d1 = kernelDensity(x, mirrored);
            double d2 = kernelDensity(-x, mirrored);
            density_dist.emplace_back(d1+d2, x);
        }

        // Sort unmirrored vector from the highest density to the lowest
        sort(density_dist.begin(), density_dist.end(), greater<pair<double, double> >());

        // Find cutoff such that indices 0 to cutoff constitute the HPD
        auto cutoff = static_cast<unsigned>(ceil(static_cast<double>(0.01*hpd_level*density_dist.size())));
        assert(cutoff <= density_dist.size());

        // Initialize hpd_lower to largest distance
        hpd_lower = *max_element(dists.begin(), dists.end());

        // Initialize hpd_upper to smallest distance
        hpd_upper = *min_element(dists.begin(), dists.end());

        // Find hpd_lower and hpd_upper
        for (unsigned i = 0; i < cutoff; ++i) {
            if (density_dist[i].second < hpd_lower) {
                hpd_lower = density_dist[i].second;
            }
            if (density_dist[i].second > hpd_upper) {
                hpd_upper = density_dist[i].second;
            }
        }

        // Sort the unmirrored vector now by the lowest to the highest dist rather than the highest to the lowest density
        sort(density_dist.begin(), density_dist.end(), [](const pair<double, double> & left, const pair<double, double> & right) {
            return left.second < right.second;
        });
        radius = density_dist[cutoff - 1].second;

        rscript = "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
        rscript += "setwd(cwd)\n";
        rscript += "pdf(\"density.pdf\")\n";

        // Save distances to R variable d
        rscript += str(format("d <- c(%.9f") % density_dist[0].second);
        for (unsigned i = 1; i < density_dist.size(); ++i) {
            rscript += str(format(", %.9f") % density_dist[i].second);
        }
        rscript += ");\n";

        // Save kernel densities to R variable fd
        rscript += str(format("fd <- c(%.9f") % density_dist[0].first);
        for (unsigned i = 1; i < density_dist.size(); ++i) {
            rscript += str(format(", %.9f") % density_dist[i].first);
        }
        rscript += ");\n";

        rscript += "dd <- density(d)\n";
        rscript += "plot(dd, type=\"l\", bty=\"n\", xlab=\"Distance\", ylab=\"Density\", main=\"Distribution of distances\");\n";
        rscript += "#points(d, fd, col=\"red\");\n";
        rscript += str(format("included <- dd$x > %g & dd$x < %g\n") % hpd_lower % hpd_upper);
        rscript += str(format("polygon(c(%g, dd$x[included], %g), c(0, dd$y[included], 0), density=NULL, col=\"lavender\");\n") % hpd_lower % hpd_upper);
        rscript += "cat(sprintf(\"bandwidth = %g\\n\", dd$bw))\n";
        rscript += "dev.off()\n";
    }

    inline void OP::calcPairwise() const {
        unsigned ntrees = _tree_summary->getNumTrees();

        if (_noisy)
            cout << "Computing pairwise distance matrix..." << endl;

        double dist_matrix[ntrees][ntrees];
        TreeManip starttm;
        TreeManip endtm;
        for (unsigned i = 0; i < ntrees - 1; i++) {
            buildTree(i, starttm);
            for (unsigned j = i+1; j < ntrees; j++) {
                buildTree(j, endtm);
                vector<Split::treeid_pair_t> ABpairs;
                vector<Split::split_pair_t> commonPairs;
                vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
                double bhvdist = calcBHVDistance(starttm, endtm, in_pairs, ABpairs, commonPairs);
                dist_matrix[i][j] = bhvdist;
                dist_matrix[j][i] = bhvdist;
            }
        }

        // Save the distance matrix to an R file
        string fn = _prefix + ".R";
        ofstream distf(fn);

        // Save the lower triangle of the pairwise distance matrix by column.
        // For example, this matrix
        // 0 0 0 0
        // 1 0 0 0
        // 2 4 0 0
        // 3 5 6 0
        // would be saved as
        // c(1,2,3,4,5,6)
        // so the indexes saved would be
        // (1,0), (2,0), (3,0), (2,1), (3,1), (3,2)
        distf << "lower_triangle_by_column <- c(";
        bool first = true;
        for (unsigned j = 0; j < ntrees-1; j++) {
            for (unsigned i = j+1; i < ntrees; i++) {
                if (first) {
                    distf << dist_matrix[i][j];
                    first = false;
                }
                else
                    distf << ", " << dist_matrix[i][j];
            }
        }
        distf << ")" << endl;
        distf << boost::str(format("m <- matrix(rep(0, %d), nrow = %d, ncol = %d, dimnames = list(paste0('t', 1:%d), paste0('t', 1:%d)))")
            % (ntrees*ntrees) % ntrees % ntrees % ntrees % ntrees) << endl;
        distf << "m[lower.tri(m, diag = FALSE)] <- lower_triangle_by_column" << endl;
        distf << "d <- as.dist(m)" << endl;
        distf.close();

        // string gofn = _prefix + "-go.R";
        // ofstream gof(gofn);
        // gof << "setwd(getSrcDirectory(function(){})[1])" << endl;
        // gof << "source(\"" << fn << "\")" << endl;
        // gof << "mapping <- cmdscale(d, 2)" << endl;
        // gof << "n <- nrow(mapping)" << endl;
        // gof << "cols <- c(rep(\"navy\", n-1),\"red\")" << endl;
        // gof << "plot(mapping, asp=1, ann=F, axes=F, col=cols, pch=16)" << endl;
        // gof.close();

        if (_noisy)
            cout << "Done." << endl;
    }

    inline void OP::randomWalk() const {
        // Carry out random walk starting from reference tree and walking for _nsteps steps
        // Each step involves adding a normal variate (mean = _step_mu, sd = _step_sigma)
        // to each internal edge length. If any edge length becomes negative, choose randomly
        // among alternative orthants.
        if (_noisy)
            cout << "Conducting random walk " << _nsteps << " long starting from reference tree" << endl;

        // Create a pseudorandom number generator
        Lot lot;
        lot.setSeed(_random_number_seed);

        // Build starting tree
        TreeManip starttm;
        buildRefTree(starttm);

        // Build end tree
        TreeManip endtm;
        buildRefTree(endtm);

        // Create vectors that must be provided to calcBHVDistance function
        vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
        vector<Split::treeid_pair_t> ABpairs;
        vector<Split::split_pair_t> commonPairs;

        // Create vector to store trees along path
        // tuple members: bhvdist, nNNIs, tree index, is rooted, newick
        vector< std::tuple<double, unsigned, unsigned, bool, string> > tvect;

        // Add starting tree to tvect
        string newick = endtm.makeNewick(9, true);
        tvect.emplace_back(0.0, 0, 0, starttm.getIsRooted(), newick);

        // Conduct random walk
        unsigned tree_index = 1;
        for (unsigned step = 0; step < _nsteps; ++step) {
            // //temporary!
            // cerr << endtm.makeNewick(5, false) << endl;

            // Take a random step
            unsigned nNNIs = endtm.randomStep(lot, _step_mu, _step_sigma);

            // Calculate BHV distance from start tree to end tree
            in_pairs.clear();
            ABpairs.clear();
            commonPairs.clear();
            double bhvdist = calcBHVDistance(starttm, endtm, in_pairs, ABpairs, commonPairs);

            // Get newick tree description of end tree
            newick = endtm.makeNewick(9, true);

            // Add end tree to tvect
            tvect.emplace_back(bhvdist, nNNIs, tree_index++, endtm.getIsRooted(), newick);
        }

        // Write trees to a file with BHV distance inserted as a nexus comment for each tree
        string fn = _prefix + ".tre";
        ofstream outf(fn);
        outf << "#nexus\n\n";
        outf << "begin trees;\n";
        for (const auto & t : tvect) {
            double dist         = get<0>(t);
            unsigned nNNIs      = get<1>(t);
            unsigned tree_index = get<2>(t);
            bool is_rooted      = get<3>(t);
            string newick       = get<4>(t);
            outf << "  tree t" << tree_index << " [bhvdist = " << setprecision(9) << dist << ", nNNIs = " << nNNIs << "] = " << (is_rooted ? "[&R] " : "[&U] ") << newick << ";\n";
        }
        outf << "end;" << endl;
        outf.close();

        // Create R file that plots BHV distance through time
        string Rfn = _prefix + ".R";
        ofstream Rf(Rfn);
        Rf << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
        Rf << "setwd(cwd)\n";
        Rf << "pdf(\"bhv_along_path.pdf\")\n";
        vector<string> bhv_str_vect;
        for (const auto & t : tvect) {
            double dist = get<0>(t);
            bhv_str_vect.push_back(str(format("%g") % dist));
        }
        Rf << str(format("bhvdist <- c(%s)\n") % boost::join(bhv_str_vect, ", "));
        Rf << "plot(bhvdist, type=\"l\", lwd=2, col=\"navy\", xlab=\"Step\", ylab=\"BHV distance\", main=\"BHV distance through time\")\n";
        Rf << "dev.off()\n";
        Rf.close();
    }

    inline void OP::calcDistanceToReference() const {
        string fn = _prefix + ".tre";
        if (_noisy)
            cout << "Writing geodesic distances (sorted) to file \"" << fn << "\"" << endl;

        int ndecimals = static_cast<int>(_precision);
        unsigned ntrees = _tree_summary->getNumTrees();
        TreeManip starttm;
        buildRefTree(starttm);
        vector< std::tuple<double, unsigned, bool, string> > tvect;
        for (unsigned i = 0; i < ntrees; i++) {
            TreeManip endtm;
            buildTree(i, endtm);
            vector<Split::treeid_pair_t> ABpairs;
            vector<Split::split_pair_t> commonPairs;
            vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
            double bhvdist = calcBHVDistance(starttm, endtm, in_pairs, ABpairs, commonPairs);
            string newick = endtm.makeNewick(9, true);
            tvect.emplace_back(bhvdist, i, endtm.getIsRooted(), newick);
        }

        // Sort trees from lowest to highest BHV distance
        sort(tvect.begin(), tvect.end());

        // Write trees to a file with BHV distance inserted as a nexus comment for each tree
        ofstream outf(fn);
        outf << "#nexus\n\n";
        outf << "begin trees;\n";
        for (const auto & t : tvect) {
            double dist         = get<0>(t);
            unsigned tree_index = get<1>(t);
            bool is_rooted      = get<2>(t);
            string newick    = get<3>(t);
            outf << "  tree t" << tree_index << " [bhvdist = " << setprecision(9) << dist << "] = " << (is_rooted ? "[&R] " : "[&U] ") << newick << ";\n";
        }
        outf << "end;" << endl;
        outf.close();
    }

#if 0
    inline void OP::calcDistanceToReference() const {
        string fn = _prefix + ".txt";
        if (_noisy)
            cout << "Writing geodesic distances (sorted) to file \"" << fn << "\"" << endl;

        int ndecimals = static_cast<int>(_precision);
        unsigned ntrees = _tree_summary->getNumTrees();
        ofstream outf(fn);
#   if defined(CLUSTER_DISTANCE)
        outf << "tree\tgeodesic\tcluster\tnewick" << endl;
#   else
        outf << "tree\tgeodesic\tnewick" << endl;
#   endif
        TreeManip starttm;
        buildTree(0, starttm);
        for (unsigned i = 1; i < ntrees; i++) {
            TreeManip endtm;
            buildTree(i, endtm);
            vector<Split::treeid_pair_t> ABpairs;
            vector<Split::split_pair_t> commonPairs;
            vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
            double bhvdist = calcBHVDistance(starttm, endtm, in_pairs, ABpairs, commonPairs);
#   if defined(CLUSTER_DISTANCE)
            double clusterdist = calcClusterDistance(starttm, endtm, in_pairs, commonPairs);
            outf << (i+1) << '\t' << setprecision(ndecimals) << bhvdist << '\t' << setprecision(ndecimals) << clusterdist << '\n';
#   else
            outf << (i+1) << '\t' << setprecision(ndecimals) << bhvdist << '\n';
#   endif
        }
        outf.close();
    }
#endif

    inline bool OP::findEmpiricalHPDwaterline(double & hpd_cutoff, double & min_log_posterior, double & max_log_posterior) const {
        hpd_cutoff = numeric_limits<double>::lowest();
        auto log_posteriors = _tree_summary->getLogPosteriors();
        bool log_posteriors_available = !log_posteriors.empty();
        if (log_posteriors_available) {
            unsigned nposteriors = static_cast<unsigned>(log_posteriors.size());

            // Sort log posteriors from highest to lowest
            vector<double> log_posteriors_sorted(log_posteriors.begin(), log_posteriors.end());
            sort(log_posteriors_sorted.begin(), log_posteriors_sorted.end(), greater<double>());

            // Empirical HPD interval is obtained from trees having log posterior > hpd_cutoff
            unsigned last = static_cast<unsigned>(floor(0.01*_radius_percent*nposteriors));
            hpd_cutoff = log_posteriors_sorted[last];

            min_log_posterior = *log_posteriors_sorted.rbegin();
            max_log_posterior = *log_posteriors_sorted.begin();
        }
        return log_posteriors_available;
    }

    inline void OP::calcFrechetMean() {
        int ndecimals = static_cast<int>(_precision);
        unsigned ntrees = _tree_summary->getNumTrees();
        TreeManip mean_tree;
        if (_noisy)
            cout << "Computing Frechet mean tree..." << endl;

        // Compute the mean tree
        unsigned number_of_iterations = computeFrechetMean(mean_tree);

        // log_posteriors will only be available if treefile was in RevBayes format
        double empirical_hpd_cutoff = numeric_limits<double>::lowest();
        double min_log_posterior = numeric_limits<double>::max();;
        double max_log_posterior = numeric_limits<double>::lowest();;
        bool log_posteriors_available = findEmpiricalHPDwaterline(empirical_hpd_cutoff, min_log_posterior, max_log_posterior);

        // Compute the variance
        double variance = 0.0;
        vector<double> bhvdists(ntrees, 0.0);
        double empirical_hpd_lower = numeric_limits<double>::max(); // smallest distance inside HPD credible set of trees
        double empirical_hpd_upper = 0.0;                           // largest distance inside HPD credible set of trees
        double smallest_distance = numeric_limits<double>::max();
        double largest_distance = 0.0;
        auto log_posteriors = _tree_summary->getLogPosteriors();

        // // These only used if log posteriors are available and _save_credible_set is true
        // TreeManip furthest_tree_inside_HPD_credible_set;
        // double furthest_distance_inside_HPD_credible_set = 0.0;
        // double furthest_inside_log_posterior = 0.0;
        // TreeManip closest_tree_not_inside_HPD_credible_set;
        // double closest_distance_not_inside_HPD_credible_set = numeric_limits<double>::max();
        // double closest_outside_log_posterior = 0.0;

        for (unsigned i = 0; i < ntrees; i++) {
            TreeManip tm;
            buildTree(i, tm);
            vector<Split::treeid_pair_t> ABpairs;
            vector<Split::split_pair_t> commonPairs;
            vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
            double bhvdist = calcBHVDistance(mean_tree, tm, in_pairs, ABpairs, commonPairs);
            variance += bhvdist*bhvdist;
            bhvdists[i] = bhvdist;
            if (bhvdist < smallest_distance) {
                smallest_distance = bhvdist;
            }
            if (bhvdist > largest_distance) {
                largest_distance = bhvdist;
            }
            // if (log_posteriors_available && _save_credible_set) {
            //     double logp = log_posteriors[i];
            //     if (logp > empirical_hpd_cutoff) {
            //         if (bhvdist > furthest_distance_inside_HPD_credible_set) {
            //             furthest_tree_inside_HPD_credible_set.setTree(tm.getTree());
            //             furthest_inside_log_posterior = logp;
            //         }
            //         if (bhvdist < empirical_hpd_lower) {
            //             empirical_hpd_lower = bhvdist;
            //         }
            //         if (bhvdist > empirical_hpd_upper) {
            //             empirical_hpd_upper = bhvdist;
            //         }
            //     }
            //     else {
            //         if (bhvdist < closest_distance_not_inside_HPD_credible_set) {
            //             closest_tree_not_inside_HPD_credible_set.setTree(tm.getTree());
            //             closest_outside_log_posterior = logp;
            //         }
            //     }
            // }
            //outf << (i+1) << '\t' << setprecision(static_cast<int>(ndecimals)) << bhvdist << '\n';
        }
        variance /= (ntrees - 1);

        // Save the mean tree, variance, HPD interval, and radius
        string nxsfn = _prefix + ".tre";
        ofstream mean_file_nexus(nxsfn);
        mean_file_nexus << "#nexus\n\n";
        mean_file_nexus << "begin trees;\n";
        mean_file_nexus << TreeManip::nexusTranslateCommand();
        mean_file_nexus << "  tree meantree = [&U] " << mean_tree.makeNewick(9, false) << ";\n";
        mean_file_nexus << "end;\n";
        mean_file_nexus.close();

        string fn = _prefix + ".R";
        ofstream mean_file(fn);
        mean_file << "# " << mean_tree.makeNewick(9, true) << endl;
        mean_file << boost::str(boost::format("# sample size = %d\n") % ntrees);
        mean_file << boost::str(boost::format("# variance = %.9f\n") % variance);
        mean_file << boost::str(boost::format("# tree length = %.9f\n") % mean_tree.calcTreeLength());
        mean_file << boost::str(boost::format("# iterations = %d\n") % number_of_iterations);
        mean_file << boost::str(boost::format("# frechet-e = %g\n") % _frechet_epsilon);
        mean_file << boost::str(boost::format("# frechet-k = %d\n") % _frechet_k);
        mean_file << boost::str(boost::format("# frechet-n = %d\n") % _frechet_n);
        mean_file << boost::str(boost::format("# pseudorandom number seed = %d\n") % _random_number_seed);
        string rscript;
        double hpd_lower;
        double hpd_upper;
        double radius;

        // Calculate kernel density and use it to fill in hpd_lower, hpd_higher, and radius
        rPlotDists(bhvdists, rscript, radius, hpd_lower, hpd_upper, _hpd_radius ? 100 : _radius_percent);

        mean_file << boost::str(boost::format("# q25 = %.9f\n") % _kde_q25);
        mean_file << boost::str(boost::format("# q75 = %.9f\n") % _kde_q75);
        mean_file << boost::str(boost::format("# IQR = %.9f\n") % (_kde_q75 - _kde_q25));
        mean_file << boost::str(boost::format("# KDE sigma = %.9f\n") % _kde_sigma);
        mean_file << boost::str(boost::format("# KDE bandwidth = %.9f\n") % _kde_bandwidth);
        mean_file << boost::str(boost::format("# %d%% HPD lower = %.9f\n") % _radius_percent % hpd_lower);
        mean_file << boost::str(boost::format("# %d%% HPD upper = %.9f\n") % _radius_percent % hpd_upper);
        mean_file << boost::str(boost::format("# %d%% radius = %.9f\n") % _radius_percent % radius);\
        // if (log_posteriors_available && _save_credible_set) {
        //     mean_file << boost::str(boost::format("# smallest log posterior = %.9f\n") % min_log_posterior);
        //     mean_file << boost::str(boost::format("# largest log posterior = %.9f\n") % max_log_posterior);
        //     mean_file << boost::str(boost::format("# %d%% empirical HPD cutoff = %.9f\n") % _radius_percent % empirical_hpd_cutoff);
        //     mean_file << boost::str(boost::format("# %d%% empirical HPD lower = %.9f\n")  % _radius_percent % empirical_hpd_lower);
        //     mean_file << boost::str(boost::format("# %d%% empirical HPD upper = %.9f\n")  % _radius_percent % empirical_hpd_upper);
        // }
        mean_file << boost::str(boost::format("# smallest distance = %.9f\n") % smallest_distance);
        mean_file << boost::str(boost::format("# largest distance = %.9f\n") % largest_distance);
        mean_file << "# The R code below creates a plot showing the kernel density for distance from mean\n";
        mean_file << boost::str(boost::format("# tree and the %d%% highest probability density credible interval\n") % _radius_percent) << rscript;
        mean_file << "\n" << rscript << endl;
        mean_file.close();

        // if (log_posteriors_available && _save_credible_set) {
        //     // Save furthest tree from mean that is nevertheless inside the HPD credible set
        //     string xfn = _prefix + "-extreme.tre";
        //     ofstream xf(xfn);
        //     xf << "#nexus\n\nbegin trees;\n";
        //     xf << "  tree mean_tree = [&R] " << mean_tree.makeNewick(9, true) << ";\n";
        //     xf << "  tree furthest_inside = [&R] [log-posterior = " << setprecision(9) << furthest_inside_log_posterior << "] " << furthest_tree_inside_HPD_credible_set.makeNewick(9, true) << ";\n";
        //     xf << "  tree closest_outside = [&R] [log-posterior = " << setprecision(9) << closest_outside_log_posterior << "] " << closest_tree_not_inside_HPD_credible_set.makeNewick(9, true) << ";\n";
        //     xf << "end;\n";
        //     xf.close();
        //
        //     string hpdfn = _prefix + "-hpd-credible-set.tre";
        //     ofstream hpdf(hpdfn);
        //     hpdf << "#nexus\n\nbegin trees;\n";
        //     for (unsigned i = 0; i < ntrees; i++) {
        //         string newick = _tree_summary->getNewick(i);
        //         double logp = log_posteriors[i];
        //         if (logp > empirical_hpd_cutoff) {
        //             hpdf << "  tree t" << (i+1) << " = [&R] [log-posterior = " << setprecision(9) << logp << "] " << newick << ";\n";
        //         }
        //     }
        //     hpdf << "end;\n";
        //     hpdf.close();
        // }

        if (_noisy) {
            cout << boost::str(format("%d iterations required (of %d max. iterations):") % number_of_iterations % _frechet_k) << endl;
            cout << "Mean tree:" << endl;
            cout << mean_tree.makeNewick(9, true) << endl;
            cout << "Variance = " << variance << endl;
            cout << "Tree length = " << setprecision(ndecimals) << mean_tree.calcTreeLength() << endl;
        }

        if (_noisy)
            cout << "Done." << endl;
    }

   // inline void OP::assertFileExists(string fn) const {
   //     // Example relative path
   //     filesystem::path relativePath = "my_file.txt";

   //     // Get the current working directory
   //     filesystem::path currentPath = filesystem::current_path();

   //     // Combine the current path and the relative path to get the absolute path
   //     filesystem::path absolutePath = currentPath / relativePath;
   // }

    inline void OP::readTrees() {
        // Create a master TreeSummary object
        _tree_summary = std::make_shared<TreeSummary>();

        // Read in reference tree if specified
        if (_refdist || _random_walk) {
            // Assume that the tree file is in Nexus format
            _tree_summary->readTreefile(_ref_tree_filename,
                true,
                _refskip,
                1,
                _refrooted,
                _refscale,
                1,
                -1,
                1);
            if (!_tree_summary->isRefTree()) {
                // Failed to read any trees, so the file may not be in Nexus format.
                // Assume next that the file is in RevBayes format.
                _tree_summary->readRevBayesTreefile(_ref_tree_filename,
                    true,
                    false,
                    _radius_percent,
                    _refskip,
                    1,
                    _refrooted,
                    _refscale,
                    -1,
                    -1,
                    1);
            }

            if (!_tree_summary->isRefTree()) {
                // Could not obtain tree from the specified _ref_tree_filename
                throw Xop("Could not obtain a tree from the specified reference tree file");
            }
        }

//        unsigned radpct = _radius_percent;
//        if (_hpd_radius) {
//            radpct = 100;
//            shrinkToHPD();
//        }

        // Read in trees
        unsigned which_treefile = 0;
        for (const auto & fn : _tree_file_names) {
            if (_hpd_radius) {
                // Assume RevBayes tree file format because _hpd_radius requires log-posterior for each tree
                TreeSummary ts;
                ts.readRevBayesTreefile(fn,
                    false,
                    _hpd_radius,
                    _radius_percent,
                    _skip[which_treefile],
                    _stride[which_treefile],
                    _rooted[which_treefile],
                    _scale_by[which_treefile],
                    _keep[which_treefile],
                    _subsample[which_treefile],
                    _subseed[which_treefile]);
                unsigned ts_ntrees = ts.getNumTrees();
                if (ts_ntrees > 0)
                    _tree_summary->absorbTreeSummary(ts);
                else
                    throw Xop(format("Number of trees saved from file \"%s\" was zero") % fn);
                if (_save_trees) {
                    string basefn = boost::filesystem::path(fn).filename().string();
                    ts.saveTrees(str(format("%s-saved.tre") % _prefix));
                }
            }
            else {
                // Assume first that the tree file is in Nexus format
                TreeSummary ts;
                ts.readTreefile(fn,
                    false,
                    _skip[which_treefile],
                    _stride[which_treefile],
                    _rooted[which_treefile],
                    _scale_by[which_treefile],
                    _keep[which_treefile],
                    _subsample[which_treefile],
                    _subseed[which_treefile]);
                unsigned ts_ntrees = ts.getNumTrees();
                if (ts_ntrees == 0) {
                    // Failed to read any trees, so the file may not be in Nexus format.
                    // Assume next that the file is in RevBayes format.
                    ts.readRevBayesTreefile(fn,
                        false,
                        false,
                        _radius_percent,
                        _skip[which_treefile],
                        _stride[which_treefile],
                        _rooted[which_treefile],
                        _scale_by[which_treefile],
                        _keep[which_treefile],
                        _subsample[which_treefile],
                        _subseed[which_treefile]);
                    ts_ntrees = ts.getNumTrees();
                }
                if (ts_ntrees > 0)
                    _tree_summary->absorbTreeSummary(ts);
                else
                    throw Xop(format("Number of trees saved from file \"%s\" was zero") % fn);
                if (_save_trees) {
                    string basefn = boost::filesystem::path(fn).filename().string();
                    ts.saveTrees(str(format("%s-saved.tre") % _prefix));
                }
            }
            ++which_treefile;
        }

        // If the tree file just read was in RevBayes format, append log_posteriors to a file
        auto log_posteriors = _tree_summary->getLogPosteriors();
        if (!log_posteriors.empty()) {
            string fn = _prefix + "-logposteriors.txt";
            ofstream outf(fn, ios::out | ios::app);
            for (const auto & log_posterior : log_posteriors) {
                outf << log_posterior << endl;
            }
            outf.close();
        }

        if (_noisy) {
            cout << "\nRead " << _tree_summary->getNumTrees() << " trees from " << _tree_file_names.size() << " tree files." << endl;
        }
    }

    inline void OP::moveAlongPath() {
        int ndecimals = static_cast<int>(_precision);
        if (_tree_summary->getNumTrees() < 2) {
            throw Xop("Must input at least 2 trees to compute tree at a particular lambda value");
        }

        // Save the tree at _lambda from the starting tree toward the ending tree
        string fn = boost::str(format("%s-lambda.tre") % _prefix);
        ofstream middle_file(fn);
        middle_file << "#nexus\n\n";
        middle_file << "begin trees;\n";
        unsigned t = 0;
        for (auto & l : _lambda) {
            if (_noisy)
                cout << "Computing tree at lambda = " << setprecision(ndecimals) << l << "..." << endl;

            TreeManip starttm;
            buildTree(0, starttm);

            TreeManip endtm;
            buildTree(1, endtm);

            displaceTreeAlongGeodesic(starttm, endtm, l);

            middle_file << str(format("  tree t%d [lambda = %.9f] = %s;\n") % (t++) % l % starttm.makeNewick(9, true));
            if (_noisy) {
                cout << "Tree at lambda = " << setprecision(ndecimals) << l << ":" << endl;
                cout << starttm.makeNewick(9, true) << endl;
            }
        }
        middle_file << "end;\n";
        middle_file.close();
    }

    inline void OP::outputForGTP() const {
        string fn1 = _prefix + "-gtp.txt";
        string fn2 = _prefix + "-gtp-taxon-mapping.txt";
        if (_noisy) {
            cout << "Writing trees in newick format to file \"" << fn1 << "\"\n";
            cout << "Writing mapping of taxon numbers to names to file \"" << fn2 << "\"" << endl;
        }

        ofstream gtpf(fn1);
        for (unsigned i = 0; i < _tree_summary->getNumTrees(); ++i) {
            TreeManip tm;
            bool isrooted = _tree_summary->isRooted(i);
            tm.buildFromNewick(_tree_summary->getNewick(i), isrooted, true);
            gtpf << tm.makeNewick(9, false) << ";\n";
        }
        gtpf.close();

        // Save the mapping of taxon numbers used in trees to taxon names
        ofstream taxmapf(fn2);
        for (unsigned i = 0; i < TreeManip::_taxon_names.size(); ++i) {
            taxmapf << str(format("%12d %s\n") % (i+1) % TreeManip::_taxon_names[i]);
        }
        taxmapf.close();
    }

    inline void OP::shrinkToHPD() {
        auto log_posteriors = _tree_summary->getLogPosteriors();
        if (log_posteriors.empty()) {
            throw Xop("Cannot compute HPD radius without log posteriors (the --hpdrad option should only be used with RevBayes tree files that record the log posterior for each tree)");
        }

        // Retain only trees inside the _radius_percent HPD credible set
        _tree_summary->retainHPDCredibleSet(_radius_percent, log_posteriors);
    }

#if defined(OLD_KF_CODE)
    if (!_quiet)
        cout << "Writing KF distances to file \"kfdists.txt\"" << endl;
    ofstream outf("kfdists.txt");
    outf << "tree distance to tree 1" << endl;
    for (unsigned i = 1; i < ntrees; i++) {
        double kfss = calcKFDistance(0, i);
        double kfdist = sqrt(kfss);
        outf << str(format("%d\t%.5f") % (i+1) % kfdist) << endl;
    }
    outf.close();
#endif

    inline void OP::run() {
        outputVersionAndSettings();
        readTrees();

        if (_gtp) {
            outputForGTP();
            return;
        }

        if (_random_walk) {
            randomWalk();
            return;
        }

        if (_snapshot) {
            moveAlongPath();
            return;
        }

        if (_pairwise) {
            calcPairwise();
            return;
        }

        if (_frechet_mean) {
            calcFrechetMean();
            return;
        }

        // Compute the geodesic distance between the first tree and all others
        if (_refdist) {
            calcDistanceToReference();
            return;
        }

        cerr << "OP was not asked to do anything with the trees in the specified treefile!" << endl;
    }
}

