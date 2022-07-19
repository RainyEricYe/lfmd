
#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <getopt.h>
#include <algorithm>

#include "option.h"
using namespace std;

typedef unsigned long ulong;
typedef map<char,double> mCharDouble;
typedef map<char, ulong> mCharUlong;
typedef map<string,ulong> mStrUlong;
typedef map<string, vector<string> > mStrStrV;

typedef pair<char,ulong> pCharUlong;
typedef pair<string,ulong> pStrUlong;
typedef pair<double, set<char> > pDoubleCharSet;
typedef pair<char,double> pCharDouble;

typedef vector< set<char> > vCharSet;
typedef vector<string> vString;

// descending sort pair
inline bool _cmpByFirst(const pDoubleCharSet &a, const pDoubleCharSet &b) { return a.first > b.first; }
inline bool _cmpBySecond(const pCharUlong &a, const pCharUlong &b) { return a.second > b.second; }
inline bool _cmpBySecond_StrUlong(const pStrUlong &a, const pStrUlong &b) { return a.second > b.second; }
inline bool _cmpBySecond_CharDouble(const pCharDouble &a, const pCharDouble &b) { return a.second > b.second; }

void replace (string &str, const string &from, const string &to, size_t more=0 );
mStrUlong fetchInDel(string &seq, char type);
vector< pair<string, ulong> > selectInDel( const mStrUlong &m );
string adjust_p(const string &qs, const Option &opt);

inline double errorRate(const char &q, const Option &opt)
{
    return pow( 10, (double)(opt.phredOffset-(int)q)/10 );
}

inline vector<double> quaToErrorRate(const string &qs, const Option &opt)
{
    vector<double> eV;
    for ( size_t i(0); i != qs.size(); i++ )  eV.push_back( errorRate(qs[i], opt) );
    return eV;
}

inline string _itoa( size_t &i )
{
    ostringstream osm;
    osm << i;
    return osm.str();
}

inline ulong countN(const string &s)
{
    ulong n(0);
    for ( size_t i(0); i != s.size(); i++ ) {
        if ( s[i] == 'N' ) n++;
    }
    return n;
}

inline bool lowQuality(const char &q, const Option &opt)
{
    return ( int(q) - opt.phredOffset < opt.baseQuaCutoff );
}

inline string reverseComplement(const string &seq)
{
    string str("");
    for ( string::const_reverse_iterator it = seq.rbegin(); it != seq.rend(); it++) {
        switch (*it) {
            case 'A': str += "T";    break;
            case 'C': str += "G";    break;
            case 'G': str += "C";    break;
            case 'T': str += "A";    break;
            case 'N':
            default:  str += "N";    break;
        }
    }

    return str;
}

inline void convertBase(string &seq, char &ref)
{
    for ( auto &i : seq ) {
        switch(i) {
            case 'a': i = 'A'; break;
            case 'c': i = 'C'; break;
            case 'g': i = 'G'; break;
            case 't': i = 'T'; break;
            case 'n': i = 'N'; break;
            case '.':
            case ',': i = ref; break;
            default:           break;
        }
    }
}

inline mCharUlong countBaseNum(const string &s)
{
    mCharUlong m;
    for ( auto &i : s ) m[i]++;
    return m;
}

inline char errorRateToChar(const double &q, const Option &opt)
{
    return char( opt.phredOffset - (int)(-10 * log(q) / log(10.0) ) );
}
#endif // MAIN_H_
