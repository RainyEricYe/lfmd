#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

//mCharDouble llh_genotype(const string &s, const string &q, const Option &opt);
map<char, vector<double> > llh_genotype(const string &s, const string &q, const Option &opt);

#endif // LIKELIHOOD_H_
