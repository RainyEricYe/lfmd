#include "main.h"

void replace (string &str, const string &from, const string &to, size_t more )
{
    size_t pos;
    size_t offset(0);
    const size_t increment = to.size();

    while ((pos = str.find(from, offset)) != string::npos) {
        str.replace(pos, from.size()+more, to);
        offset = pos + increment;
    }
}

mStrUlong fetchInDel(string &s, char type)
{
    mStrUlong m;
    size_t p(0);

    while ( (p=s.find(type,p)) != string::npos ) {
        string len("");
        size_t offset(0);
        size_t length;

        for ( size_t i(p+1); i != s.size(); i++ ) {
            if ( isdigit( s[i] ) ) {
                len += s[i];
            }
            else {
                offset = i;
                break;
            }
        }

        length = atoi(len.c_str());
        string indel( s.begin()+offset, s.begin()+offset+length );
        m[indel]++;

        s.replace(p, offset-p+length, "");
  //      cout << offset << ' '<< p << ' ' << length << ' ' << offset-p+length << ' ' << indel
    //       << "\nnew seq: " << s << endl;

    }

    return m;
}

vector<pStrUlong> selectInDel( const mStrUlong &m )
{
//    vector<pStrUlong> v( m.begin(), m.end() );
    vector<pStrUlong> v;

    for ( auto &p : m ) {
        if ( countN(p.first) / (double)p.first.size() > 0.5 ) continue;
        v.push_back(p);
    }

    if ( v.size() > 1 )
        sort( v.begin(), v.end(), _cmpBySecond_StrUlong );
/*
    vector< vector<pStrUlong>::iterator > fv;

    for ( vector<pStrUlong>::iterator it = v.begin(); it != v.end(); it++ ) {
        if ( countN(it->first) / (double)it->first.size() > 0.5 ) fv.push_back(it);
    }

    for ( auto & f : fv ) v.erase(f);
*/
    return v;
}

string adjust_p(const string &qs, const Option &opt)
{
    ostringstream o;

    map<double, vector<char> > m;
    vector<double> v = quaToErrorRate(qs, opt);

    if ( v.size() > 1 ) {
        sort( v.begin(), v.end() );

        ulong total( v.size() );
        for ( size_t i(0); i != total; i++ ) {
            m[ v[i] ].push_back( errorRateToChar(v[i], opt) );
        }

        for ( auto & q : qs) {
            double e = errorRate(q, opt);
            o << m[ e ].back();
            m[ e ].pop_back();
        }

        return o.str();
    }
    else {
        return qs;
    }
}
