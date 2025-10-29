#pragma once

namespace op {

    class TreeSummary {
        public:
                                        TreeSummary() = default;
                                        ~TreeSummary() = default;

            static string               scaleEdgeLengths(const string & newick, bool rooted, double scaler);
            void                        readRevBayesTreefile(const string & filename, unsigned skip, unsigned stride, bool rooted, double scaler);
            void                        readTreefile(const string & filename, unsigned skip, unsigned stride, bool rooted, double scaler);
            //void                        showSummary() const;
            const vector<double> &      getLogPosteriors() const;
            unsigned                    getNumTrees() const;
            Tree::SharedPtr             getTree(unsigned index) const;
            string                      getNewick(unsigned index);
            bool                        isRooted(unsigned index);
            void                        clear();

        private:

            //Split::treemap_t            _treeIDs;
            vector<string>              _newicks;
            vector<double>              _log_posteriors;
            vector<bool>                _is_rooted;

        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
        };

    inline unsigned TreeSummary::getNumTrees() const {
        return static_cast<unsigned>(_newicks.size());
    }

    inline const vector<double> & TreeSummary::getLogPosteriors() const {
        return _log_posteriors;
    }

    inline Tree::SharedPtr TreeSummary::getTree(unsigned index) const {
        if (index >= _newicks.size())
            throw Xop("getTree called with index >= number of stored trees");

        TreeManip tm;

        // build the tree
        tm.buildFromNewick(_newicks[index], false, true);

        return tm.getTree();
    }

    inline bool TreeSummary::isRooted(unsigned index) {
        if (index >= _is_rooted.size())
            throw Xop("isRooted called with index >= number of stored trees");

        return _is_rooted[index];
    }

    inline string TreeSummary::getNewick(unsigned index) {
        if (index >= _newicks.size())
            throw Xop("getNewick called with index >= number of stored trees");

        return _newicks[index];
    }

    inline void TreeSummary::clear() {
        //_treeIDs.clear();
        _newicks.clear();
        }

    inline string TreeSummary::scaleEdgeLengths(const string & newick, bool rooted, double scaler) {
        TreeManip tm;
        tm.buildFromNewick(newick, rooted, true);
        tm.scaleAllEdgeLengths(scaler);
        return tm.makeNewick(9, false);
    }

    inline void TreeSummary::readRevBayesTreefile(const string & filename, unsigned skip, unsigned stride, bool rooted, double scaler) {
        ifstream inf(filename.c_str());
        stringstream buffer;
        buffer << inf.rdbuf();
        string line;

        // Get header
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

        unsigned t = 0;
        while (getline(buffer, line)) {
            bool do_sample = (t >= skip);
            do_sample = do_sample && (((t - skip) % stride) == 0);
            if (do_sample) {
                boost::algorithm::split(parts, line, boost::is_any_of("\t"));
                nparts = static_cast<unsigned>(parts.size());
                assert((is_combined_treefile && nparts == 6) || nparts == 5);
                _is_rooted.push_back(rooted);
                auto log_posterior = lexical_cast<double>(parts[is_combined_treefile ? 2 : 1]);
                _log_posteriors.emplace_back(log_posterior);
                string & newick = parts[is_combined_treefile ? 5 : 4];

                bool do_scale = static_cast<bool>(scaler != 1.0);
                if (do_scale)
                    _newicks.emplace_back(scaleEdgeLengths(newick, rooted, scaler));
                else
                    _newicks.emplace_back(newick);
            }
            t++;
        }
    }

    inline void TreeSummary::readTreefile(const string & filename, unsigned skip, unsigned stride, bool rooted, double scaler) {
        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

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

        auto numTaxaBlocks = static_cast<int>(nexusReader.GetNumTaxaBlocks());
        for (int i = 0; i < numTaxaBlocks; ++i) {
            NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
            //string taxaBlockTitle = taxaBlock->GetTitle();

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
                sort(old_names.begin(), old_names.end());
                sort(new_names.begin(), new_names.end());
                bool ok = true;
                for (unsigned j = 0; j < taxaBlock->GetNumTaxonLabels(); ++j) {
                    if (old_names[j] != new_names[j]) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) {
                    cerr << boost::str(format("%12s %12s %12s") % "taxon" % "prev" % "new") << endl;
                    for (unsigned j = 0; j < taxaBlock->GetNumTaxonLabels(); ++j) {
                        cerr << boost::str(format("%12d %12s %12s") % (j+1) % old_names[j] % new_names[j]) << endl;
                    }
                    throw Xop("Taxon labels in taxa block do not match those from a previous taxa block");
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
                        if (do_sample) {
                            const NxsFullTreeDescription & d = treesBlock->GetFullTreeDescription(t);

                            // If the full tree description is "processed," then node indices will be 1 + index of taxon in
                            // taxa block
                            assert(d.IsProcessed());

                            bool is_rooted = d.IsRooted();
                            if (is_rooted != rooted) {
                                throw Xop(format("Tree in \"%s\" was %s but you led me to expect they would be %s") % filename % (is_rooted ? "rooted" : "unrooted") % (rooted ? "rooted" : "unrooted"));
                            }
                            _is_rooted.push_back(is_rooted);

                            // store the newick tree description
                            string newick = d.GetNewick();;

                            bool do_scale = static_cast<bool>(scaler != 1.0);
                            if (do_scale)
                                _newicks.push_back(scaleEdgeLengths(newick, is_rooted, scaler));
                            else
                                _newicks.push_back(newick);
                        } // do_sample
                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop
        } // TAXA block loop

        // No longer any need to store raw data from the nexus file
        nexusReader.DeleteBlocksFromFactories();
    }

}
