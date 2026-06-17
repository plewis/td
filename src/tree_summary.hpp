#pragma once
#include "lot.hpp"

namespace op {

    class TreeSummary {
        public:
            enum class TreeFileType {
                NEXUS,
                REVBAYES,
                BPP,
                OTHER
            };
                                        TreeSummary();
                                        ~TreeSummary();

#if defined(ALWAYS_ROOTED)
            static string               scaleEdgeLengths(const string & newick, double scaler);
            void                        readBPPTreefile(bool reftree, unsigned skip, unsigned stride, double scaler, int keep, int subsample, unsigned subseed);
            void                        readRevBayesTreefile(const string & filename, bool reftree, bool hpdrad, double radpct, unsigned skip, unsigned stride, double scaler, int keep, int subsample, unsigned subseed);
            void                        readTreefile(const string & filename, bool reftree, unsigned skip, unsigned stride, double scaler, int keep, int subsample, unsigned subseed);
#else
            static string               scaleEdgeLengths(const string & newick, bool rooted, double scaler);
            void                        readBPPTreefile(bool reftree, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed);
            void                        readRevBayesTreefile(const string & filename, bool reftree, bool hpdrad, double radpct, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed);
            void                        readTreefile(const string & filename, bool reftree, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed);
#endif
            void                        saveTrees(const string & filename) const;
            bool                        readFileIntoBuffer(const string & filename);

            // Utility functions
            TreeFileType                treeFileTypeFromBuffer() const;
            void                        stripThetasFromBPPTree(string & newick) const;
            void                        retainHPDCredibleSet(unsigned hpdpct, const vector<double> & log_posteriors);
            void                        randomlySample(unsigned sample_size, unsigned rnseed);
            void                        absorbTreeSummary(const TreeSummary & other);
            void                        clear();

            // Accessors
            const vector<double> &      getLogPosteriors() const;
            unsigned                    getNumTrees() const;
            string                      getRefNewick() const;
            Tree::SharedPtr             getRefTree() const;
            bool                        isRefTree() const;
#if defined(ALWAYS_ROOTED)
#else
            bool                        isRefRooted() const;
#endif
            Tree::SharedPtr             getTree(unsigned index) const;
            string                      getNewick(unsigned index) const;
            bool                        isIncludedInHPD(unsigned index) const;
#if defined(ALWAYS_ROOTED)
#else
            bool                        isRooted(unsigned index) const;
#endif

        private:

            //Split::treemap_t            _treeIDs;
            stringstream                _buffer;
            string                      _reftree;
            vector<string>              _newicks;
            vector<double>              _log_posteriors;
#if defined(ALWAYS_ROOTED)
#else
            bool                        _is_refrooted;
            vector<bool>                _is_rooted;
#endif
        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
        };

#if defined(ALWAYS_ROOTED)
    inline TreeSummary::TreeSummary() {
        _reftree = "";
    }
#else
    inline TreeSummary::TreeSummary() : _is_refrooted(false) {
        _reftree = "";
    }
#endif

    inline TreeSummary::~TreeSummary() {
    }

    inline unsigned TreeSummary::getNumTrees() const {
        return static_cast<unsigned>(_newicks.size());
    }

    inline const vector<double> & TreeSummary::getLogPosteriors() const {
        return _log_posteriors;
    }

    inline Tree::SharedPtr TreeSummary::getRefTree() const {
        TreeManip tm;

        // build the tree
#if defined(ALWAYS_ROOTED)
        tm.buildFromNewick(_reftree, true, true);
#else
        tm.buildFromNewick(_reftree, _is_refrooted, true);
#endif

        return tm.getTree();
    }

    inline Tree::SharedPtr TreeSummary::getTree(unsigned index) const {
#if defined(ALWAYS_ROOTED)
#else
        assert(_is_rooted.size() == _newicks.size());
#endif
        if (index >= _newicks.size())
            throw Xop("getTree called with index >= number of stored trees");

        TreeManip tm;

        // build the tree
#if defined(ALWAYS_ROOTED)
        tm.buildFromNewick(_newicks[index], true, true);
#else
        tm.buildFromNewick(_newicks[index], _is_rooted[index], true);
#endif
        return tm.getTree();
    }

    inline bool TreeSummary::isRefTree() const {
        return !_reftree.empty();
    }

#if defined(ALWAYS_ROOTED)
#else
    inline bool TreeSummary::isRefRooted() const {
        return _is_refrooted;
    }
#endif

#if defined(ALWAYS_ROOTED)
#else
    inline bool TreeSummary::isRooted(unsigned index) const {
        if (index >= _is_rooted.size())
            throw Xop("isRooted called with index >= number of stored trees");

        return _is_rooted[index];
    }
#endif

    inline string TreeSummary::getRefNewick() const {
        return _reftree;
    }

    inline string TreeSummary::getNewick(unsigned index) const {
        if (index >= _newicks.size())
            throw Xop("getNewick called with index >= number of stored trees");

        return _newicks[index];
    }

    inline void TreeSummary::clear() {
        //_treeIDs.clear();
        _reftree = "";
        _newicks.clear();
        _log_posteriors.clear();
#if defined(ALWAYS_ROOTED)
#else
        _is_refrooted = false;
        _is_rooted.clear();
#endif
    }

#if defined(ALWAYS_ROOTED)
    inline string TreeSummary::scaleEdgeLengths(const string & newick, double scaler) {
        TreeManip tm;
        tm.buildFromNewick(newick, true, true);
        if (scaler != 1.0)
            tm.scaleAllEdgeLengths(scaler);
        return tm.makeNewick(9, false);
    }
#else
    inline string TreeSummary::scaleEdgeLengths(const string & newick, bool rooted, double scaler) {
        TreeManip tm;
        tm.buildFromNewick(newick, rooted, true);
        if (scaler != 1.0)
            tm.scaleAllEdgeLengths(scaler);
        return tm.makeNewick(9, false);
    }
#endif

    inline TreeSummary::TreeFileType TreeSummary::treeFileTypeFromBuffer() const {
        // Get reference to data stored in buffer
        const string & s = _buffer.str();

        // The _buffer should not be empty when this function is called
        assert(!s.empty());

        // Copy the first 50 characters of s
        string first_fifty = s.substr(0, 50);

        //temporary!
        cerr << "\n********** length of s = " << s.size() << " **********" << endl;
        cerr << "\n********** first_fifty = " << first_fifty << " **********" << endl;

        // If the file is in NEXUS format, the first 50 characters should start with "#nexus" (case insensitive)
        regex nexus_re(R"(#[Nn][Ee][Xx][Uu][Ss][\s\S]+)");
        if (regex_match(first_fifty, nexus_re)) {
            return TreeFileType::NEXUS;
        }

        // If the file is in RevBayes format, the first 50 characters should start
        // "Iteration\tPosterior\tLikelihood\tPrior\tpsi"
        regex revbayes_re(R"(Iteration\s+Posterior\s+Likelihood\s+Prior\s+psi[\s\S]+)");
        if (regex_match(first_fifty,  revbayes_re)) {
            return TreeFileType::REVBAYES;
        }

        // If the file is in BPP format, there should be theta values following the # character
        regex bpp_re(R"([\s\S]+?\s+[#][.0-9]+[:]\s+[.0-9]+[,)][\s\S]+)");
        if (regex_match(first_fifty, bpp_re)) {
            return TreeFileType::BPP;
        }

        return TreeFileType::OTHER;
    }

    inline void TreeSummary::stripThetasFromBPPTree(string & newick) const {
        // Define a regular expression that matches BPP theta values
        regex bpp_theta_re(R"((\s+#[.0-9]+?)(?=[:;]))");
        string bpp_theta_replacement = "";

        // Use regex_replace to create a new string with all matches replaced
        newick = regex_replace(newick, bpp_theta_re, bpp_theta_replacement);

        regex bpp_space_re(R"([:]\s+([.0-9]+))");
        string bpp_space_replacement = ":$1";

        // Use regex_replace to create a new string with all matches replaced
        newick = regex_replace(newick, bpp_space_re, bpp_space_replacement);
    }

    inline bool TreeSummary::readFileIntoBuffer(const string & filename) {
        // Read entire tree file into a buffer
        ifstream inf(filename.c_str());
        if (!static_cast<bool>(inf)) {
            throw Xop(format("File \"%s\" does not exist.") % filename);
        }
        else {
            //temporary!
            cerr << "File \"" << filename.c_str() << "\" exists." << endl;
        }
        _buffer.clear();
        _buffer << inf.rdbuf();
        inf.close();

        // //temporary!
        // string result = _buffer.str();
        // std::cerr << "Resulting string: " << result << endl;

        return true;
    }

#if defined(ALWAYS_ROOTED)
    inline void TreeSummary::readBPPTreefile(bool reftree, unsigned skip, unsigned stride, double scaler, int keep, int subsample, unsigned subseed)
#else
    inline void TreeSummary::readBPPTreefile(bool reftree, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed)
#endif
        {
        // Specifying stride, keep, subsample, and subseed makes no sense if reading in reference tree
        bool default_stride    = (stride == 1);
        bool default_keep      = (keep == -1);
        bool default_subsample = (subsample == -1);
        bool default_subseed   = (subseed == 0);
        assert( !reftree || ( default_stride && default_keep && default_subsample && default_subseed ) );

        // Object should be empty when this function is called
        // To accumulate trees from different treefiles, use absorbTreeSummary function
        // to pull trees from a different TreeSummary object into this object
        assert(_newicks.empty());
        assert(_log_posteriors.empty());
#if defined(ALWAYS_ROOTED)
#else
        assert(_is_rooted.empty());
#endif

        // Negative keep says to keep all trees
        if (keep < 0) {
            keep = std::numeric_limits<int>::max();
        }

        // Go through buffer line by line
        unsigned t = 0; // indexes all trees in the file, even if they are not saved
        int nkept = 0; // keeps track of the number of trees from the file that are saved
        string newick;
        while (getline(_buffer, newick)) {
            // Determine whether this line holds a tree that should be saved
            bool do_sample = (t >= skip);
            do_sample = do_sample && (((t - skip) % stride) == 0);
            do_sample = do_sample && (nkept < keep);
            if (do_sample) {
                // // Check whether this line has the expected number of tab-separated fields
                if (reftree) {
#if defined(ALWAYS_ROOTED)
                    stripThetasFromBPPTree(newick);
                    _reftree = scaleEdgeLengths(newick, scaler);
#else
                    _is_refrooted = rooted;
                    stripThetasFromBPPTree(newick);
                    _reftree = scaleEdgeLengths(newick, rooted, scaler);
#endif
                    return;
                }
                else {
                    // Store rooting information for this tree
#if defined(ALWAYS_ROOTED)
                    // Store newick tree description for this tree
                    stripThetasFromBPPTree(newick);
                    _newicks.emplace_back(scaleEdgeLengths(newick, scaler));
#else
                    _is_rooted.push_back(rooted);
                    // Store newick tree description for this tree
                    stripThetasFromBPPTree(newick);
                    _newicks.emplace_back(scaleEdgeLengths(newick, rooted, scaler));
#endif
                    nkept++;
                }
            }
            t++;
        }

        // If subsampling, preserve only subsample randomly-chosen trees
        if (subsample != -1) {
            randomlySample(subsample, subseed);
        }
    }

#if defined(ALWAYS_ROOTED)
    inline void TreeSummary::readRevBayesTreefile(const string & filename, bool reftree, bool hpdrad, double radpct, unsigned skip, unsigned stride, double scaler, int keep, int subsample, unsigned subseed)
#else
    inline void TreeSummary::readRevBayesTreefile(const string & filename, bool reftree, bool hpdrad, double radpct, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed)
#endif
    {
        // Specifying stride, keep, subsample, and subseed makes no sense if reading in reference tree
        assert(!reftree || ((stride == 1) && (keep == -1 || keep == 1) && (subsample == -1) && (subseed == 0)));

        // Object should be empty when this function is called
        // To accumulate trees from different treefiles, use absorbTreeSummary function
        // to pull trees from a different TreeSummary object into this object
        assert(_newicks.empty());
        assert(_log_posteriors.empty());
#if defined(ALWAYS_ROOTED)
#else
        assert(_is_rooted.empty());
#endif

        // Negative keep says to keep all trees
        if (keep < 0) {
            keep = std::numeric_limits<int>::max();
        }

        // Read entire tree file into a buffer
        ifstream inf(filename.c_str());
        stringstream buffer;
        buffer << inf.rdbuf();

        // Get header
        string line;
        getline(buffer, line);

        // Split header at tabs
        vector<string> parts;
        boost::algorithm::split(parts, line, boost::is_any_of("\t"));
        auto nparts = static_cast<unsigned>(parts.size());

        // Assuming this is indeed a RevBayes tree file
        bool is_combined_treefile = false;
        if (nparts == 5) {
            if (parts[4] != "psi") {
                throw Xop("Expecting 5th column of RevBayes tree file to be labeled \"psi\"");
            }
        }
        else if (nparts == 6) {
            if (parts[5] != "psi") {
                throw Xop("Expecting 6th column of RevBayes combined tree file to be labeled \"psi\"");
            }
            is_combined_treefile = true;
        }

        // Go through buffer line by line
        unsigned t = 0; // indexes all trees in the file, even if they are not saved
        int nkept = 0; // keeps track of the number of trees from the file that are saved
        while (getline(buffer, line)) {
            // Determine whether this line holds a tree that should be saved
            bool do_sample = (t >= skip);
            do_sample = do_sample && (((t - skip) % stride) == 0);
            do_sample = do_sample && (nkept < keep);
            if (do_sample) {
                // Check whether this line has the expected number of tab-separated fields
                boost::algorithm::split(parts, line, boost::is_any_of("\t"));
                nparts = static_cast<unsigned>(parts.size());
                assert((is_combined_treefile && nparts == 6) || nparts == 5);
                if (reftree) {
#if defined(ALWAYS_ROOTED)
                    string & newick = parts[is_combined_treefile ? 5 : 4];
                    _reftree = scaleEdgeLengths(newick, scaler);
#else
                    _is_refrooted = rooted;
                    string & newick = parts[is_combined_treefile ? 5 : 4];
                    _reftree = scaleEdgeLengths(newick, rooted, scaler);
#endif
                    return;
                }
                else {
                    // Store rooting information for this tree
#if defined(ALWAYS_ROOTED)
#else
                    _is_rooted.push_back(rooted);
#endif

                    // Store log-posterior for this tree
                    double log_posterior = 0.0;
                    try {
                        log_posterior = lexical_cast<double>(parts[is_combined_treefile ? 2 : 1]);
                    } catch (boost::bad_lexical_cast &) {
                        throw Xop(format("Error parsing log-posterior from RevBayes tree file \"%s\" line %d: \"%s\"") % filename % (t+1) % parts[is_combined_treefile ? 2 : 1]);
                    }
                    _log_posteriors.emplace_back(log_posterior);

                    // Store newick tree description for this tree
                    string & newick = parts[is_combined_treefile ? 5 : 4];
#if defined(ALWAYS_ROOTED)
                    _newicks.emplace_back(scaleEdgeLengths(newick, scaler));
#else
                    _newicks.emplace_back(scaleEdgeLengths(newick, rooted, scaler));
#endif
                    nkept++;
                }
            }
            t++;
        }

        // If hpdrad, delete any trees that are not in the radpct HPD credible set
        if (hpdrad) {
            vector<double> logposteriors(_log_posteriors.begin(), _log_posteriors.end());
            retainHPDCredibleSet(radpct, logposteriors);
        }

        // If subsampling, preserve only subsample randomly-chosen trees
        if (subsample != -1) {
            randomlySample(subsample, subseed);
        }
    }

#if defined(ALWAYS_ROOTED)
    inline void TreeSummary::readTreefile(const string & filename, bool reftree, unsigned skip, unsigned stride, double scaler, int keep, int subsample, unsigned subseed)
#else
    inline void TreeSummary::readTreefile(const string & filename, bool reftree, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed)
#endif
    {
        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

        // Specifying stride, keep, subsample, and subseed makes no sense if reading in reference tree
        assert(!reftree || ((stride == 1) && (keep == -1 || keep == 1) && (subsample == -1) && (subseed == 1)));

        // Negative keep says to keep all trees
        if (keep < 0) {
            keep = std::numeric_limits<int>::max();
        }

        //MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
        MultiFormatReader nexusReader(-1, NxsReader::IGNORE_WARNINGS);

        // Both of these needed to suppress "storing read block" messages
        // see NxsReader::statusMessage in nxsreader.cpp
        nexusReader.SetAlwaysReportStatusMessages(false);
        nexusReader.SetWarningOutputLevel(NxsReader::SUPPRESS_WARNINGS_LEVEL);

        try {
            nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(...) {
            nexusReader.DeleteBlocksFromFactories();
            throw;
        }

        int nkept = 0;
        auto numTaxaBlocks = static_cast<int>(nexusReader.GetNumTaxaBlocks());
        for (int i = 0; i < numTaxaBlocks; ++i) {
            NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
            //string taxaBlockTitle = taxaBlock->GetTitle();

            // taxon_index_map will remain empty unless the ordering of taxa in taxaBlock
            // differs from that in TreeManip::_taxon_names, in which case taxon_index_map
            // translates the index in taxaBlock to the index implied
            // by TreeManip::_taxon_names.
            map<unsigned, unsigned> taxon_index_map;

            if (TreeManip::_taxon_names.empty()) {
                // Copy taxon labels into the TreeManip::_taxon_names vector
                TreeManip::_taxon_names.resize(taxaBlock->GetNumTaxonLabels());
                TreeManip::_taxon_map.clear();
                for (unsigned j = 0; j < TreeManip::_taxon_names.size(); ++j) {
                    string taxon_name = taxaBlock->GetTaxonLabel(j);
                    TreeManip::_taxon_names[j] = taxon_name;
                    TreeManip::_taxon_map[taxon_name] = j;
                }
                assert(taxaBlock->GetNumActiveTaxa() == TreeManip::_taxon_names.size());
            }
            else {
                // Check that taxon labels in taxa block match those in TreeManip::_taxon_names
                vector<string> old_names(TreeManip::_taxon_names.begin(), TreeManip::_taxon_names.end());
                vector<string> new_names;
                for (unsigned j = 0; j < taxaBlock->GetNumTaxonLabels(); ++j) {
                    new_names.push_back(taxaBlock->GetTaxonLabel(j));
                }
                bool ok = true;
                for (unsigned j = 0; j < taxaBlock->GetNumTaxonLabels(); ++j) {
                    if (old_names[j] != new_names[j]) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) {
                    cerr << boost::str(format("File name: \"%s\"") % filename) << endl;
                    cerr << boost::str(format("%12s %12s %12s") % "taxon" % "prev" % "new") << endl;
                    for (unsigned j = 0; j < taxaBlock->GetNumTaxonLabels(); ++j) {
                        cerr << boost::str(format("%12d %12s %12s") % (j+1) % old_names[j] % new_names[j]) << endl;
                    }

                    // Create a mapping of taxon index in new taxa block to taxon index in prev taxa block
                    // Note that taxon_index_map stores indices of taxa, not taxon numbers (the taxon
                    // number is its index plus 1).
                    taxon_index_map.clear();
                    for (unsigned j = 0; j < taxaBlock->GetNumTaxonLabels(); ++j) {
                        string no_underscores = std::regex_replace(new_names[j], std::regex("_"), " ");
                        taxon_index_map[j] = TreeManip::_taxon_map.at(no_underscores) ;
                    }

                    // Check whether the names are the same, just in a different order
                    sort(old_names.begin(), old_names.end());
                    sort(new_names.begin(), new_names.end());
                    bool ok2 = true;
                    for (unsigned j = 0; j < taxaBlock->GetNumTaxonLabels(); ++j) {
                        if (old_names[j] != new_names[j]) {
                            ok2 = false;
                            break;
                        }
                    }

                    if (ok2) {
                        cerr << "Taxon labels in this taxa block are identical to those from a previous taxa block" << endl;
                        cerr << "but are in a different order. Correcting newick descriptions to match previous taxa block." << endl;
                    }
                    else {
                        throw Xop("Taxon labels in taxa block do not match those from a previous taxa block");
                    }
                }
            }

            const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
            for (unsigned j = 0; j < nTreesBlocks; ++j) {
                const NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, j);
                unsigned ntrees = treesBlock->GetNumTrees();
                if (skip < ntrees) {
                    //cout << "Trees block contains " << ntrees << " tree descriptions.\n";
                    for (unsigned t = skip; t < ntrees; ++t) {
                        // NxsFullTreeDescription::TreeDescFlags possibilities
                        // (stored in flags data member):
                        // NXS_IS_ROOTED_BIT                    = 0x0001,
                        // NXS_HAS_SOME_EDGE_LENGTHS_BIT        = 0x0002,
                        // NXS_MISSING_SOME_EDGE_LENGTHS_BIT    = 0x0004,
                        // NXS_EDGE_LENGTH_UNION                = 0x0006,
                        // NXS_INT_EDGE_LENGTHS_BIT             = 0x0008,
                        // NXS_HAS_ALL_TAXA_BIT                 = 0x0010,
                        // NXS_HAS_NHX_BIT                      = 0x0020,
                        // NXS_HAS_DEG_TWO_NODES_BIT            = 0x0040,
                        // NXS_HAS_POLYTOMY_BIT                 = 0x0080,
                        // NXS_HAS_INTERNAL_NAMES_BIT           = 0x0100,
                        // NXS_HAS_NEW_INTERNAL_NAMES_BIT       = 0x0200,
                        // NXS_KNOWN_INTERNAL_NAMES_BIT         = 0x0400,
                        // NXS_SOME_ZERO_EDGE_LEN_BIT           = 0x0800,
                        // NXS_SOME_NEGATIVE_EDGE_LEN_BIT       = 0x1000,
                        // NXS_TREE_PROCESSED                   = 0x2000

                        bool do_sample =  (((t - skip) % stride) == 0);
                        do_sample = do_sample && (nkept < keep);
                        if (do_sample) {
                            const NxsFullTreeDescription & d = treesBlock->GetFullTreeDescription(t);

                            // If the full tree description is "processed," then node indices will be 1 + index of taxon in
                            // taxa block
                            assert(d.IsProcessed());

#if defined(ALWAYS_ROOTED)
#else
                            bool is_rooted = d.IsRooted();
                            if (is_rooted != rooted) {
                                throw Xop(format("Trees in \"%s\" were %s but you specified \"rooted = %s\" (possibly by default)") % filename % (is_rooted ? "rooted" : "unrooted") % (rooted ? "yes" : "no"));
                            }
#endif

                            // store the newick tree description
                            string newick = d.GetNewick();;
                            if (!taxon_index_map.empty()) {
                                TreeManip::stripComments(newick);
                                newick = TreeManip::renumberNewick(newick, taxon_index_map);
                            }

                            if (reftree) {
#if defined(ALWAYS_ROOTED)
                                _reftree = scaleEdgeLengths(newick, scaler);
#else
                                _is_refrooted = is_rooted;
                                _reftree = scaleEdgeLengths(newick, is_rooted, scaler);
#endif
                                break;
                            }
                            else {
#if defined(ALWAYS_ROOTED)
                                // Sampling all trees in file, go ahead and save to _newicks
                                _newicks.push_back(scaleEdgeLengths(newick, scaler));
#else
                                _is_rooted.push_back(is_rooted);
                                // Sampling all trees in file, go ahead and save to _newicks
                                _newicks.push_back(scaleEdgeLengths(newick, is_rooted, scaler));
#endif
                                nkept++;
                            }
                        } // do_sample
                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop

            // No longer need taxon_index_map if it was created
            taxon_index_map.clear();
        } // TAXA block loop

        if (subsample != -1) {
            randomlySample(subsample, subseed);
        }

        // No longer any need to store raw data from the nexus file
        nexusReader.DeleteBlocksFromFactories();
    }

    inline void TreeSummary::randomlySample(unsigned sample_size, unsigned rnseed) {
        // Initialize the pseudorandom number generator
        Lot lot;
        lot.setSeed(rnseed);

        // Do not randomly sample if there are fewer trees than sample size
        unsigned ntrees = static_cast<unsigned>(_newicks.size());
        if (ntrees <= sample_size) {
            return;
        }

        // Move _newicks to a new home
        vector<string> newicks;
        newicks.insert(newicks.end(), make_move_iterator(_newicks.begin()), make_move_iterator(_newicks.end()));
        _newicks.erase(_newicks.begin(), _newicks.end());

        // Move _log_posteriors to a new home
        vector<double> log_posteriors;
        log_posteriors.insert(log_posteriors.end(), make_move_iterator(_log_posteriors.begin()), make_move_iterator(_log_posteriors.end()));
        _log_posteriors.erase(_log_posteriors.begin(), _log_posteriors.end());

#if defined(ALWAYS_ROOTED)
        for (unsigned i = 0; i < sample_size; ++i) {
            unsigned j = lot.randint(0, newicks.size()-1);

            _newicks.emplace_back(newicks[j]);
            _log_posteriors.push_back(log_posteriors[j]);
        }
#else
        // Move _is_rooted to a new home
        vector<bool> is_rooted;
        is_rooted.insert(is_rooted.end(), make_move_iterator(_is_rooted.begin()), make_move_iterator(_is_rooted.end()));
        _is_rooted.erase(_is_rooted.begin(), _is_rooted.end());

        for (unsigned i = 0; i < sample_size; ++i) {
            unsigned j = lot.randint(0, newicks.size()-1);

            _newicks.emplace_back(newicks[j]);
            _log_posteriors.push_back(log_posteriors[j]);
            _is_rooted.push_back(is_rooted[j]);
        }
#endif
    }

    inline void TreeSummary::retainHPDCredibleSet(unsigned hpdpct, const vector<double> & log_posteriors) {
        // Create a vector of indices for all included trees
        unsigned nlogposteriors = static_cast<unsigned>(log_posteriors.size());
        vector<unsigned> indices(nlogposteriors);
        for (unsigned i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }

        // Sort the vector of indices by descending posterior probability
        sort(indices.begin(), indices.end(), [&log_posteriors](unsigned i1, unsigned i2) {
            return log_posteriors[i1] > log_posteriors[i2];
        });

        // Determine which element of indices is the last one we will consider
        unsigned last = static_cast<unsigned>(ceil(0.01*hpdpct*nlogposteriors));
        assert (last <= nlogposteriors);

        // Make a temporary copy of _newicks and _log_posteriors
        vector<string> orig_newicks(_newicks.begin(), _newicks.end());
        vector<double> orig_log_posteriors(_log_posteriors.begin(), _log_posteriors.end());

        // Clear _newicks and _log_posteriors
        _newicks.clear();
        _newicks.resize(last);
        _log_posteriors.clear();
        _log_posteriors.resize(last);

        for (unsigned i = 0; i < last; ++i) {
            unsigned j = indices[i];
            _newicks[i] = orig_newicks[j];
            _log_posteriors[i] = orig_log_posteriors[j];
        }
    }

    inline void TreeSummary::absorbTreeSummary(const TreeSummary & other) {
        // Appends _newicks, _log_posteriors, and _is_rooted from other to this
        // Throw exception if both this and other have a reference tree defined because
        // there is no way to append a reference tree
        if (!_reftree.empty() && !other._reftree.empty()) {
            throw Xop("Error in absorbTreeSummary: reference trees exist in both source and destination TreeSummary objects");
        }
        _newicks.insert(_newicks.end(), other._newicks.begin(), other._newicks.end());
        _log_posteriors.insert(_log_posteriors.end(), other._log_posteriors.begin(), other._log_posteriors.end());
#if defined(ALWAYS_ROOTED)
#else
        _is_rooted.insert(_is_rooted.end(), other._is_rooted.begin(), other._is_rooted.end());
#endif
    }

    inline void TreeSummary::saveTrees(const string & fn) const {
        // Saves trees (and log posteriors, if available) to a file named fn
        size_t ntrees = _newicks.size();
#if defined(ALWAYS_ROOTED)
#else
        assert(_is_rooted.size() == ntrees);
#endif
        assert(_log_posteriors.empty() || _log_posteriors.size() == ntrees);
        ofstream hpdf(fn);
        hpdf << "#nexus\n\n";
        hpdf << "begin trees;\n";
        hpdf << TreeManip::nexusTranslateCommand();
        for (unsigned i = 0; i < ntrees; i++) {
            string newick = getNewick(i);
#if defined(ALWAYS_ROOTED)
            if (_log_posteriors.empty()) {
                hpdf << "  tree t" << (i+1) << " = [&R] " << newick << ";\n";
            }
            else {
                double logp = _log_posteriors[i];
                hpdf << "  tree t" << (i+1) << " = [&R] [log-posterior = " << setprecision(9) << logp << "] " << newick << ";\n";
            }
#else
            if (_log_posteriors.empty()) {
                hpdf << "  tree t" << (i+1) << " = [&" << (_is_rooted[i] ? "R" : "U") << "] " << newick << ";\n";
            }
            else {
                double logp = _log_posteriors[i];
                hpdf << "  tree t" << (i+1) << " = [&" << (_is_rooted[i] ? "R" : "U") << "] [log-posterior = " << setprecision(9) << logp << "] " << newick << ";\n";
            }
#endif
        }
        hpdf << "end;\n";
        hpdf.close();

    }
}
