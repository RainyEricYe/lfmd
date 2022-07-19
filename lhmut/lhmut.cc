/*
 */

#include "main.h"
#include "likelihood.h"
#include "option.h"

#include "iomanip"

int main( int argc, char **argv )
{
    Option opt;
    opt.parse(argc, argv);

    ifstream inf(opt.infileName.c_str());
    ofstream outf(opt.outfileName.c_str());

    if ( !inf.is_open() )  cerr << "open error: " << opt.infileName << endl, exit(1);
    if ( !outf.is_open() ) cerr << "open error: " << opt.outfileName << endl, exit(1);

    outf <<
        "Chrom\tTemplate\tPos\tDepths\tMuts\t"
        "Tcount\tCcount\tGcount\tAcount\t"
        "inscount\tdelcount\tNcount"
        << endl;

    string line;
    while ( getline(inf, line) ) {
        string chr, baseStr, quaStr;
        ulong pos(0), depth(0);
        char ref;

        istringstream itm(line);
        itm >> chr >> pos >> ref >> depth >> baseStr >> quaStr;

        replace(baseStr, "^", "", 1); // delete lable after ^, be careful of ^$
        replace(baseStr, "$", "");
        convertBase(baseStr, ref); // to upper & , . to ref

        mStrUlong mInsNum = fetchInDel(baseStr, '+');
        mStrUlong mDelNum = fetchInDel(baseStr, '-');

        mCharUlong mBaseNum = countBaseNum(baseStr);

        depth = baseStr.size();
        if ( depth != quaStr.size() ) {
            cerr << "parse error: length of bases and qualities are different "
                << chr << ' ' << pos << endl;
            continue;
        }

        ulong mutN(0);
        for ( auto &p : mBaseNum ) {
            if ( p.first == ref || p.first == '*' ) continue;  // not ref, ins, del
            mutN += p.second;
        }

        outf << chr << '\t' << ref << '\t' << pos << '\t'
             << depth << '\t' << mutN << '\t';

        string bs = "TCGA";
        for ( auto &b : bs ) {
            mCharUlong::const_iterator i = mBaseNum.find(b);
            outf << ( i == mBaseNum.end() ? 0 : i->second ) << '\t';
        }

        vector<pStrUlong> vInsNum = selectInDel(mInsNum);
        vector<pStrUlong> vDelNum = selectInDel(mDelNum);

        outf << ( vInsNum.empty() ? 0 : vInsNum[0].second ) << '\t';
        outf << ( vDelNum.empty() ? 0 : vDelNum[0].second ) << '\t';

        mCharUlong::const_iterator nit = mBaseNum.find('N');
        outf << ( nit == mBaseNum.end() ? 0 : nit->second );


        // detail info below

        for ( auto &p : vInsNum ) outf << "\tI+" << p.first << ':' << p.second;
        for ( auto &p : vDelNum ) outf << "\tD-" << p.first << ':' << p.second;

        map<char, vector<double> > ntPF = llh_genotype(baseStr, quaStr, opt);

        outf << std::scientific << std::setprecision(6);
        for ( auto &p : ntPF ) {
            outf << '\t' << p.first << ":" << p.second.at(0) << ":" << p.second.at(1);
        }

        outf << endl;
    }

    inf.close();
    outf.close();

    return 0;
}
