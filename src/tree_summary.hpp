#pragma once
#include "lot.hpp"

namespace op {

    class TreeSummary {
        public:
                                        TreeSummary();
                                        ~TreeSummary();

            static string               scaleEdgeLengths(const string & newick, bool rooted, double scaler);
            void                        readRevBayesTreefile(const string & filename, bool reftree, bool hpdrad, double radpct, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed);
            void                        readTreefile(const string & filename, bool reftree, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed);
            void                        saveTrees(const string & filename) const;

            // Utility functions
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
            bool                        isRefRooted() const;
            Tree::SharedPtr             getTree(unsigned index) const;
            string                      getNewick(unsigned index) const;
            bool                        isRooted(unsigned index) const;
            bool                        isIncludedInHPD(unsigned index) const;

        private:

            //Split::treemap_t            _treeIDs;
            string                      _reftree;
            bool                        _is_refrooted;
            vector<string>              _newicks;
            vector<double>              _log_posteriors;
            vector<bool>                _is_rooted;

        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
        };

    inline TreeSummary::TreeSummary() : _is_refrooted(false) {
        _reftree = "";
    }

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
        tm.buildFromNewick(_reftree, _is_refrooted, true);

        return tm.getTree();
    }

    inline Tree::SharedPtr TreeSummary::getTree(unsigned index) const {
        assert(_is_rooted.size() == _newicks.size());
        if (index >= _newicks.size())
            throw Xop("getTree called with index >= number of stored trees");

        TreeManip tm;

        // build the tree
        tm.buildFromNewick(_newicks[index], _is_rooted[index], true);

        return tm.getTree();
    }

    inline bool TreeSummary::isRefTree() const {
        return !_reftree.empty();
    }

    inline bool TreeSummary::isRefRooted() const {
        return _is_refrooted;
    }

    inline bool TreeSummary::isRooted(unsigned index) const {
        if (index >= _is_rooted.size())
            throw Xop("isRooted called with index >= number of stored trees");

        return _is_rooted[index];
    }

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
        _is_refrooted = false;
        _newicks.clear();
        _log_posteriors.clear();
        _is_rooted.clear();
    }

    inline string TreeSummary::scaleEdgeLengths(const string & newick, bool rooted, double scaler) {
        TreeManip tm;
        tm.buildFromNewick(newick, rooted, true);
        if (scaler != 1.0)
            tm.scaleAllEdgeLengths(scaler);
        return tm.makeNewick(9, false);
    }

    inline void TreeSummary::readRevBayesTreefile(const string & filename, bool reftree, bool hpdrad, double radpct, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed) {
        // Specifying stride, keep, subsample, and subseed makes no sense if reading in reference tree
        assert(!reftree || ((stride == 1) && (keep == -1 || keep == 1) && (subsample == -1) && (subseed == 0)));

        // Object should be empty when this function is called
        // To accumulate trees from different treefiles, use absorbTreeSummary function
        // to pull trees from a different TreeSummary object into this object
        assert(_newicks.empty());
        assert(_log_posteriors.empty());
        assert(_is_rooted.empty());

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
        unsigned nkept = 0; // keeps track of the number of trees from the file that are saved
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
                    _is_refrooted = rooted;
                    string & newick = parts[is_combined_treefile ? 5 : 4];
                    _reftree = scaleEdgeLengths(newick, rooted, scaler);
                    return;
                }
                else {
                    // Store rooting information for this tree
                    _is_rooted.push_back(rooted);

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
                    _newicks.emplace_back(scaleEdgeLengths(newick, rooted, scaler));
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

    inline void TreeSummary::readTreefile(const string & filename, bool reftree, unsigned skip, unsigned stride, bool rooted, double scaler, int keep, int subsample, unsigned subseed) {
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

        unsigned nkept = 0;
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

                            bool is_rooted = d.IsRooted();
                            if (is_rooted != rooted) {
                                throw Xop(format("Trees in \"%s\" were %s but you specified \"rooted = %s\" (possibly by default)") % filename % (is_rooted ? "rooted" : "unrooted") % (rooted ? "yes" : "no"));
                            }

                            // store the newick tree description
                            string newick = d.GetNewick();;
                            if (!taxon_index_map.empty()) {
                                newick = TreeManip::renumberNewick(newick, taxon_index_map);
                            }

                            if (reftree) {
                                _is_refrooted = is_rooted;
                                _reftree = scaleEdgeLengths(newick, is_rooted, scaler);
                                break;
                            }
                            else {
                                _is_rooted.push_back(is_rooted);
                                // Sampling all trees in file, go ahead and save to _newicks
                                _newicks.push_back(scaleEdgeLengths(newick, is_rooted, scaler));
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

        vector<bool> is_rooted;
        is_rooted.insert(is_rooted.end(), make_move_iterator(_is_rooted.begin()), make_move_iterator(_is_rooted.end()));
        _is_rooted.erase(_is_rooted.begin(), _is_rooted.end());

        //temporary!
        if (true) {
            // Randomly sample subsample trees from tmp_newicks
            ofstream tmpf(str(format("selection-%d.txt") % rnseed));

            for (unsigned i = 0; i < sample_size; ++i) {
                unsigned j = lot.randint(0, newicks.size()-1);

                //temporary!
                tmpf << j << "\n";

                _newicks.push_back(newicks[j]);
                _log_posteriors.push_back(log_posteriors[j]);
                _is_rooted.push_back(is_rooted[j]);
            }

            //temporary!
            tmpf.close();
        }
        else {
            // Create a vector of (uniform-random-variate, index) pairs for all included trees
            vector<pair<double,unsigned> > indices(ntrees);
            for (unsigned i = 0; i < ntrees; ++i) {
                indices[i] = make_pair(lot.uniform(), i);
            }

            // Sort the vector to randomize the indices
            sort(indices.begin(), indices.end());

            for (unsigned i = 0; i < sample_size; ++i) {
                unsigned j = indices[i].second;
                _newicks.emplace_back(newicks[j]);
                _log_posteriors.emplace_back(log_posteriors[j]);
                _is_rooted.push_back(is_rooted[j]);
            }
        }
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
        _is_rooted.insert(_is_rooted.end(), other._is_rooted.begin(), other._is_rooted.end());
    }

    inline void TreeSummary::saveTrees(const string & fn) const {
        // Saves trees (and log posteriors, if available) to a file named fn
        size_t ntrees = _newicks.size();
        assert(_is_rooted.size() == ntrees);
        assert(_log_posteriors.empty() || _log_posteriors.size() == ntrees);
        ofstream hpdf(fn);
        hpdf << "#nexus\n\n";
        hpdf << "begin trees;\n";
        hpdf << TreeManip::nexusTranslateCommand();
        for (unsigned i = 0; i < ntrees; i++) {
            string newick = getNewick(i);
            if (_log_posteriors.empty()) {
                hpdf << "  tree t" << (i+1) << " = [&" << (_is_rooted[i] ? "R" : "U") << "] " << newick << ";\n";
            }
            else {
                double logp = _log_posteriors[i];
                hpdf << "  tree t" << (i+1) << " = [&" << (_is_rooted[i] ? "R" : "U") << "] [log-posterior = " << setprecision(9) << logp << "] " << newick << ";\n";
            }
        }
        hpdf << "end;\n";
        hpdf.close();

    }
}
