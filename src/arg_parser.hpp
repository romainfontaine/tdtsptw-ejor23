#ifndef TDTSPTW_ARG_PARSER_HPP
#define TDTSPTW_ARG_PARSER_HPP

#include <iostream>
#include <list>
#include <unordered_set>
#include <sstream>
#include <unordered_map>
#include <string>
#include <optional>
#include <functional>
#include <filesystem>

using namespace std;

namespace low_level_parser {

    template<typename T>
    bool parse(const string &v, const string &key, T &v_, string &error_msg) = delete;

    template<>
    bool parse<string>(const string &v, const string &key, string &v_, string &error_msg) {
        v_ = v;
        return true;
    }

    inline char stochar(const string &s, size_t *idx) {
        if (s.size() == 1) {
            *idx = 1;
            return s[0];
        }
        *idx = s.size() ^ 1; // set error
        return 0; // does not matter
    }

    inline bool stobool(const string &s, size_t *idx) {
        if (s == "0"s || s == "1"s || s == "true"s || s == "false"s) {
            *idx = s.size();
            return s == "1"s || s == "true"s;
        }
        *idx = s.size() ^ 1; // set error
        return false; // does not matter
    }

#define PARSER(T, PARSE_FUNC)                       \
template<>                                          \
inline bool parse<T>(const string &v,               \
        const string &key, T &v_,                   \
        string &error_msg) {                        \
    size_t idx;                                     \
    v_ = PARSE_FUNC(v, &idx);                       \
    if (idx != v.size()){                           \
        ostringstream oss;                          \
        oss << "Expecting a '"#T"' value for '"     \
        << key << "'='" << v << "'";                \
        error_msg = oss.str();                      \
        return false;                               \
    }                                               \
    return true;                                    \
}

    PARSER(uint, stoul)

    PARSER(int, stoi)

    PARSER(double, stod)

    PARSER(float, stof)

    PARSER(char, stochar)

    PARSER(bool, stobool)

#undef PARSER
}

class ArgParser {
    string name_, version_, help_, bin_;
    list <string> positional;
    unordered_map<string, string> arguments;

    bool keyExists(const string &k) const {
        return arguments.count(k) != 0;
    }

    string getArg(const string &key) {
        // Precondition : key exists in arguments
        string v(move(arguments[key]));
        arguments.erase(key);
        return v;
    }

    template<typename T>
    void parsingErrorSet(const T &v, const string &name, const unordered_set<T> &s) {
        ostringstream oss;
        oss << "value '" << v << "' of " << name << " should belong to set {";
        bool first = true;
        for (const auto &s_: s) {
            if (!first)
                oss << ", ";
            else
                first = false;
            oss << "'" << s_ << "'";
        }
        oss << "}";
        parsingError(oss.str());
    }

    template<typename T>
    T parseString(const string &value, const string &key) {
        T v;
        string err;
        if (!low_level_parser::parse<T>(value, key, v, err))
            parsingError(err);
        return v;
    }

public:
    void parsingError(const string &msg) {
        cerr << "Parsing error: " << msg << "." << endl << endl;
        help();
        exit(1);
    }

    ArgParser(const string &name, const string &version, const string &help, const string &bin) : name_(name),
                                                                                                  version_(version),
                                                                                                  help_(help),
                                                                                                  bin_(bin) {}

    ArgParser(const ArgParser &) = delete;

    void help(const bool &full = false) const {
        cout << name_ << " version " << version_ << endl;
	cout << help_ << endl;
    }

    void parse(int argc, char *argv[]) {
        for (uint i = 0; i < static_cast<uint>(argc); i++) {
            string s(argv[i]);
            if (s[0] == '-') {
                const auto &eq_it = s.find('=');
                const string &key = s.substr(0, eq_it);
                if (key.size() <= 1)
                    parsingError("empty argument");
                string val;
                if (eq_it == string::npos) {
                    val = "true";
                } else {
                    val = s.substr(eq_it + 1);
                    if (!val.size())
                        parsingError("empty value for argument '" + key + "'");
                }
                if (arguments.count(key))
                    parsingError("argument '" + key + "' was provided more than once");
                arguments[key] = val;
            } else {
                positional.push_back(move(s));
            }
        }
    }

    template<typename T>
    T getPositional(const string &arg_name, function<bool(const T &)> fn, const string &error_msg) {
        const T &v = getPositional<T>(arg_name);
        if (!fn(v))
            parsingError(error_msg);
        return v;
    }

    template<typename T>
    T getPositional(const string &arg_name, const initializer_list<T> &s) {
        return getPositional<T>(arg_name, unordered_set<T>(s));
    }

    template<typename T>
    T getPositional(const string &arg_name, const unordered_set<T> &s) {
        const T v(getPositional<T>(arg_name));
        if (s.count(v) == 0)
            parsingErrorSet(v, arg_name, s);
        return v;
    }

    void ensureFileExists(const string &filename, const string &arg_name) {
        if (!filesystem::exists(filename))
            parsingError("file '" + filename + "' for '" + arg_name + "' does not exist");
    }

    string getPositionalFile(const string &arg_name) {
        const string &filename = getPositional<string>(arg_name);
        ensureFileExists(filename, arg_name);
        return filename;
    }

    template<typename T>
    T getPositional(const string &arg_name) {
        if (positional.empty())
            parsingError("missing positional argument '" + arg_name + "'");
        const string v_(move(positional.front()));
        positional.pop_front();
        return parseString<T>(v_, arg_name);
    }

    uint remainingPositionals() const {
        return positional.size();
    }

    template<typename T>
    T getArgument(const string &key, const T &default_v, const bool &required,
                  const initializer_list<T> &s) {
        return getArgument<T>(key, default_v, required, unordered_set<T>(s));
    }

    template<typename T>
    T getArgument(const string &key, const T &default_v, const bool &required,
                  const unordered_set<T> &s) {
        T v(getArgument<T>(key, default_v, required));
        if (s.count(v) == 0)
            parsingErrorSet(v, key, s);
        return v;
    }

    template<typename T>
    T getArgument(const string &key, const T &default_v, const bool &required,
                  function<bool(const T &)> fn, const string &error_msg) {
        T v(getArgument<T>(key, default_v, required));
        if (!fn(v))
            parsingError(error_msg);
        return v;
    }

    template<typename T>
    T getArgument(const string &key, const T &default_v, const bool &required = false) {
        if (!keyExists(key)) {
            if (!required)
                return default_v;
            else
                parsingError("missing required argument '" + key + "'");
        }
        return parseString<T>(getArg(key), key);
    }

    bool hasArgument(const string &key) const {
        return keyExists(key);
    }

    void done() {
        // If some args have not been used, crash.
        if (!positional.empty() || !arguments.empty()) {
            string msg("too many arguments: ");
            if (!positional.empty()) {
                for (const auto &a: positional)
                    msg += "'" + a + "', ";
            }
            if (!arguments.empty()) {
                for (const auto &[a, b]: arguments)
                    msg += "'" + a + "=" + b + "', ";
            }
            msg.erase(msg.length() - 2);
            parsingError(msg);
        }
    }

    void dump() {
        cout << "Positional : " << endl;
        for (const auto &a: positional)
            cout << "'" << a << "', ";
        cout << endl;
        cout << "Arguments : " << endl;
        for (const auto &[a, b]: arguments)
            cout << "'" << a << "'='" << b << "', ";
        cout << endl;
    }
};

#endif //TDTSPTW_ARG_PARSER_HPP

