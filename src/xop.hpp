#pragma once

using namespace std;
using namespace boost;

class Xop : public std::runtime_error
    {
    public:
                    Xop() : std::runtime_error("error") {}
                    Xop(const format & f) : runtime_error(str(f)) {}
        explicit    Xop(const string & s) : runtime_error(s) {}
                    ~Xop() override = default;
    };
