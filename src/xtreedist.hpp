#pragma once

using namespace std;
using namespace boost;

#include <boost/format.hpp>

namespace treedist {

    class XTreeDist : public std::exception {
        public:
                                XTreeDist() throw() {}
                                XTreeDist(const string s) throw() : _msg() {_msg = s;}
                                XTreeDist(const format & f) throw() : _msg() {_msg = str(f);}
            virtual             ~XTreeDist() throw() {}
            const char *        what() const throw() {return _msg.c_str();}

        private:

            string         _msg;
    };

}
