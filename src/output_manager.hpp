#pragma once

#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace boost;

namespace treedist {

    class OutputManager {
        public:
            typedef std::shared_ptr< OutputManager > SharedPtr;

                  OutputManager();
                  ~OutputManager();
            
            void   outputNewline(bool flush = false) const;
            void   outputConsole(const string & s, bool flush = false) const;
            void   outputConsole(const format & fmt, bool flush = false) const;
            void   outputConsole(const program_options::options_description & description) const;
    };
    
    inline OutputManager::OutputManager() {
    }

    inline OutputManager::~OutputManager() {
    }

    inline void OutputManager::outputNewline(bool flush) const {
        cout << endl;
        if (flush)
            cout.flush();
    }
    
    inline void OutputManager::outputConsole(const string & s, bool flush) const {
        cout << s;
        if (flush)
            cout.flush();
    }
    
    inline void OutputManager::outputConsole(const format & fmt, bool flush) const {
        cout << str(fmt);
        if (flush)
            cout.flush();
    }
    
    inline void OutputManager::outputConsole(const program_options::options_description & description) const {
        cout << description << endl;
    }
    
}
