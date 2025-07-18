#pragma once

class XStrom : public std::runtime_error
    {
    public:
                    XStrom() : std::runtime_error("error") {}
        explicit    XStrom(const string & s) : std::runtime_error(s) {}
                    ~XStrom() override = default;
    };
