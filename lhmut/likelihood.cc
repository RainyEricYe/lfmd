#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <algorithm>

#include "main.h"
#include "option.h"
#include "likelihood.h"
#include <boost/math/distributions/chi_squared.hpp>

// needed by alglib
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "ap.h"

struct fn_data {
    std::string         base;
    std::vector<double> errRateV;
};

// composite log likelihood: l_c(theta)
// mat theta is a column vector which has 4 elements for A, C, G, T, respectively.
//
double composite_LogLikelihood (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta )
{
    double l_c(0.0);

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i]/3;
        switch( base[i] ) {
            case 'A': l_c += log( (1-4*e) * theta[0] + e ); break;
            case 'C': l_c += log( (1-4*e) * theta[1] + e ); break;
            case 'G': l_c += log( (1-4*e) * theta[2] + e ); break;
            case 'T': l_c += log( (1-4*e) * theta[3] + e ); break;
            default: cerr << "unknown base in " << base << endl,exit(1);
        }
    }

    return l_c;
}

// -composite score function: -U_c(theta)
// return a column vector
//
alglib::real_1d_array composite_score (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta )
{
    alglib::real_1d_array U_c = "[0,0,0,0]";

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i]/3;

        switch( base[i] ) {
            case 'A': U_c[0] -= (1-4*e) / ( (1-4*e)*theta[0] + e );  break;
            case 'C': U_c[1] -= (1-4*e) / ( (1-4*e)*theta[1] + e );  break;
            case 'G': U_c[2] -= (1-4*e) / ( (1-4*e)*theta[2] + e );  break;
            case 'T': U_c[3] -= (1-4*e) / ( (1-4*e)*theta[3] + e );  break;
            default: cerr << "unknown base in " << base << endl, exit(1);
        }
    }

    return U_c;
}

// gradient optimization
void function1_grad (
        const alglib::real_1d_array  &x,
        double                       &func,
        alglib::real_1d_array        &grad,
        void                         *opt_data )
{
    fn_data* objfn_data = reinterpret_cast<fn_data*>(opt_data);
    const std::string &base = objfn_data->base;
    const std::vector<double> &errRateV = objfn_data->errRateV;

    func = -composite_LogLikelihood( base, errRateV, x);
    grad = composite_score( base, errRateV, x );
}

string initAlleleFreq (
        mCharUlong       &fr,
        const double     &depth,
        const char       &except_b )
{
    ostringstream s;
    s << '[';

    if ( depth == fr[except_b] ) {
        for ( auto b : "ACGT" ) {
            b == except_b ? s << 0.0 : s << 0.333333333;
            b == 'T' ? s << ']' : s << ',';
            if ( b == 'T' ) break;
        }
    }
    else {
        for ( auto b : "ACGT" ) {
            b == except_b ? s << 0.0 : s << fr[b] / (depth - fr[except_b]);
            b == 'T' ? s << ']' : s << ',';
            if ( b== 'T' ) break;
        }
    }
    return s.str();
}

string initAlleleFreq (
        mCharUlong       &fr,
        double           depth,
        const set<char>  &except_bs )
{
    ostringstream s;

    for ( auto b : "ACGT" ) {
        set<char>::const_iterator it = except_bs.find(b);
        if ( it != except_bs.end() ) {
            fr[b] = 0;
            depth -= fr[b];
        }

        if ( b== 'T' ) break;
    }

    if ( depth == 0 ) {
        s  << "[0.25,0.25,0.25,0.25]";
        return s.str();
    }

    s << '[';
    for ( auto b : "ACGT" ) {
        set<char>::const_iterator it = except_bs.find(b);
        it != except_bs.end() ? s << 0.0 : s << fr[b] / depth;  // depth here don't contain num of except_bs
        b == 'T' ? s << ']' : s << ',';
        if ( b== 'T' ) break;
    }

    return s.str();
}

string _upBoundary(const char except_b)
{
    ostringstream s;
    s << '[';

    for ( auto b : "ACGT" ) {
        b == except_b ? s << 0 : s << 1;
        b == 'T' ? s << ']' : s << ',';
        if ( b == 'T' ) break;
    }
    return s.str();
}

string _upBoundary(const set<char> except_bs)
{
    ostringstream s;
    s << '[';

    for ( auto b : "ACGT" ) {
        set<char>::const_iterator it = except_bs.find(b);
        it != except_bs.end() ? s << 0 : s << 1;
        b == 'T' ? s << ']' : s << ',';
        if ( b == 'T' ) break;
    }
    return s.str();
}

map<char, vector<double> > llh_genotype(const string &s, const string &q, const Option &opt)
{
//    mCharDouble ntP; // nt => pvalue
    map<char, vector<double> > ntPF; // nt => pvalue, fraction

    boost::math::chi_squared X2_dist(1);

    mCharUlong fr;
    for ( auto b : "ACGTN" ) {
        fr[b] = 0;
        if ( b == 'N' ) break;
    }

    string new_s(""), new_q("");
    for ( size_t i(0); i != s.size(); i++ ) {
        if ( lowQuality(q[i], opt) || s[i] == 'N' || s[i] == '*' )  continue; // fr['N'] == 0
        fr[ s[i] ]++;
        new_s += s[i];
        new_q += q[i];
    }

    double depth(new_s.size());
    vector<double> errV = quaToErrorRate(new_q, opt);

    if ( depth == 0 ) return ntPF;

    fn_data data;
    data.base = new_s;
    data.errRateV = errV;

    // var for alglib
    alglib::minbleicstate state;
    alglib::minbleicreport rep;
    double epsg(0.000001);
    double epsf(0.0);
    double epsx(0.0);
    alglib::ae_int_t maxits(0);

    // constraint: sum of frequency of 4 alleles == 1
    alglib::real_2d_array c = "[[1,1,1,1,1]]";  // sum of four allele == 1
    alglib::integer_1d_array ct = "[0]";    // equal
    alglib::real_1d_array bndl = "[0,0,0,0]";  // lower boundary

    // four allele maximize
    double cl_4(0.0);
    try {
        string AFstr = initAlleleFreq(fr, depth, 'N');
        alglib::real_1d_array alg_x = AFstr.c_str();
        alglib::real_1d_array bndu = "[1,1,1,1]";

        alglib::minbleiccreate(alg_x, state);
        alglib::minbleicsetlc(state, c, ct);
        alglib::minbleicsetbc(state, bndl, bndu);
        alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
        alglib::minbleicoptimize(state, function1_grad, NULL, &data );
        alglib::minbleicresults(state, alg_x, rep);

        if ( opt.debug ) {
            printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
            printf("%s\n", alg_x.tostring(20).c_str());
        }

        cl_4 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
        if ( opt.debug ) cout << "cl_4: " << setprecision(20) << cl_4 << endl;
    }
    catch ( alglib::ap_error &e ) {
        cerr << "catch error: " << e.msg << " at seq[" << new_s << "] qua[" << new_q << "]" << endl;
    }

    map<char, string> init_V;
    map<char, string> bndu_V;

    for ( auto b : "ACGT" ) {
        init_V[b] = initAlleleFreq(fr, depth, b);
        bndu_V[b] = _upBoundary(b);
        if ( b == 'T' ) break;
    }

    for ( mCharUlong::const_iterator it = fr.begin(); it != fr.end(); it++ )
    {
        if ( it->second < opt.minSupOnEachStrand || it->second/depth < opt.minFractionInFam ) continue;

        double cl_3(0.0);
        try {
            alglib::real_1d_array alg_x = init_V[ it->first ].c_str();
            alglib::real_1d_array bndu = bndu_V[ it->first ].c_str();

            alglib::minbleiccreate(alg_x, state);
            alglib::minbleicsetlc(state, c, ct);
            alglib::minbleicsetbc(state, bndl, bndu);
            alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
            alglib::minbleicoptimize(state, function1_grad, NULL, &data );
            alglib::minbleicresults(state, alg_x, rep);

            if ( opt.debug ) {
                printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
                printf("%s\n", alg_x.tostring(20).c_str());
            }

            cl_3 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
            if ( opt.debug ) cout << "cl_3: " << cl_3 << endl;
        }
        catch ( alglib::ap_error &e ) {
            cerr << "catch error: " << e.msg << " at seq[" << new_s << "] qua[" << new_q
                << "] for base[" << it->first << "]" << endl;
        }

        if ( cl_4 - cl_3 > opt.lhrGapCutoff ) {
//            ntP[ it->first ] = 1 - boost::math::cdf(X2_dist, 2*(cl_4 - cl_3) );
            ntPF[ it->first ].push_back( 1 - boost::math::cdf(X2_dist, 2*(cl_4 - cl_3)) );
        }
    }

    if ( ntPF.size() == 1 ) {
        ntPF[ ntPF.begin()->first ].push_back(1.0);
        return ntPF;
    }
    else if ( ntPF.size() > 1 ) {
        set<char> except_bs;
        for ( auto b : "ACGT" ) {
            map<char, vector<double> >::const_iterator it = ntPF.find(b);
            if ( it == ntPF.end() ) {  // not in ntPF
                except_bs.insert(b);
            }

            if ( b == 'T' ) break;
        }

        string AFstr = initAlleleFreq(fr, depth, except_bs);
        alglib::real_1d_array alg_x = AFstr.c_str();

        try {
            string upBnd = _upBoundary(except_bs);
            alglib::real_1d_array bndu = upBnd.c_str();

            alglib::minbleiccreate(alg_x, state);
            alglib::minbleicsetlc(state, c, ct);
            alglib::minbleicsetbc(state, bndl, bndu);
            alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
            alglib::minbleicoptimize(state, function1_grad, NULL, &data );
            alglib::minbleicresults(state, alg_x, rep);

            if ( opt.debug ) {
                printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
                printf("%s\n", alg_x.tostring(20).c_str());
            }

            cl_4 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
            if ( opt.debug ) cout << "cl_4: " << setprecision(20) << cl_4 << endl;
        }
        catch ( alglib::ap_error &e ) {
            cerr << "catch error: " << e.msg << " at seq[" << new_s << "] qua[" << new_q << "]" << endl;
        }

        string st = "ACGT";
        map<char, double> mBaseFrac;
        for ( int i(0); i != 4; i++ ) {
            mBaseFrac[ st[i] ] = alg_x[i];
        }

        for ( auto &p : ntPF ) {
            ntPF[ p.first ].push_back( mBaseFrac[p.first] );
        }

        return ntPF;
    }
    else if ( ntPF.size() > 4 ) {
        cerr << "ntPF contain unknown base" << endl, exit(1);
    }
    else {
        return ntPF;
    }
}

