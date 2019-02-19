
#include "SVDQuestTripartitionScorer.hpp"

#include <phylonaut/CladeExtractor.hpp>
#include <phylonaut/ASTRALCladeExtractor.hpp>

#include <newick.hpp>
#include <phylonaut/wASTRAL.hpp>
#include <util/Logger.hpp>
#include <util/Timer.hpp>

#include <limits>
#include <fstream>
#include <cmath>
#if defined(_WIN32) || defined(_WIN64)
/* We are on Windows */
# define strtok_r strtok_s
#endif


void write_nex(string infile, string outfile, TaxonSet& ts) {
  ifstream inf(infile);


  int ntaxa, nchars;
  inf >> ntaxa >> nchars;

  ofstream of(outfile);

  INFO << ts.size() << endl;
  of << "#nexus\n";
  of << "Set AllowPunct=yes;" << endl;
  of << "BEGIN TAXA;"<< endl;
  of << "DIMENSIONS NTAX=" << ts.size()<< ";" << endl;

  of << "    TAXLABELS" << endl;

  for (int i = 0; i < ts.size(); i++) {
    of <<"    " << i << endl;
  }
  of << "    ;" << endl << "END;" << endl << endl;

  of << "BEGIN CHARACTERS;" << endl;
  of << "    DIMENSIONS NCHAR=" << nchars << ";" << endl;
  of << "    FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=.;" << endl;
  of << "    MATRIX" << endl <<  endl;

  string str;
  vector<string> sequences(ts.size());


  while (getline(inf, str)) {
    if (!str.size()) {
      continue;
    }
    string taxname = strtok(&(str[0]), " \t");
    cout << "taxon " << taxname << endl;
    if (!taxname.size()) {
      continue;
    }
    string s = strtok(NULL, " \t\n");
    sequences[ts[taxname]] += s;
  }

  for (int i = 0; i < ts.size(); i++) {
    of << i << "    " << sequences[i] << endl;
  }


  of << endl << ";" << endl << "END;";
  of.close();
  inf.close();
}

void parse_quartet(QuartetDict& qd, char* c) {
  char* saveptr;
  int i = 0;
  Taxon taxa[4];

  while(char* token=strtok_r(c, "|,", &saveptr)) {
    c = NULL;
    taxa[i] = atoi(token) - 1;
    i++;
  }

  qd.set(taxa[0], taxa[1], taxa[2], taxa[3], 1);
}

void SVDQuestTripartitionScorer::runPaup(Config& conf) {

  Timer::start("GetQuartetWeights");
  string nexfile = outname + ".svdquest.nex";
  string quartetfile = outname + ".svdquest.quartets";
  string paupfile = outname + ".pauptree";
  INFO << "Writing nexus file " << nexfile << endl;
  write_nex(alignmentfile, nexfile, ts);
  INFO << "Done" << endl;

  string command = "exe " + nexfile + ";\r\n svd evalQuartets=all \r\n qfile=" + quartetfile + "\r\n qformat=qmc \r\n ambigs=missing; \r\n savetrees file="+ paupfile + " format=Newick root=yes;\n";



  INFO << "Running PAUP* " << command << endl;
  FILE* paupstream;
#if defined(_WIN32) || defined(_WIN64)
  paupstream = _popen("paup -n", "w");
#else
  if (wine)
	  paupstream = popen("wine paup4c -n", "w");
  else
	  paupstream = popen("paup -n", "w");
#endif
  
  fputs(command.c_str(), paupstream);

  fclose(paupstream);

  INFO << "DONE" << endl;


  ifstream infile(quartetfile);

  qd = new QuartetDict(ts, "");

  Quartet q(ts);
  string s;
  double w;
  while(!infile.eof()) {
    getline(infile, s);
    if (s.size() == 0)
      continue;
    parse_quartet(*qd, &(s[0]));

  }

  Timer::stop("GetQuartetWeights");



    Timer::start("GetAdditionalClades");

    ifstream pauptree_stream(paupfile);
    string pauptree;
    pauptree_stream >> pauptree;
    pauptree_stream.close();
    pauptree = unmap_newick_names(pauptree, ts);
    ofstream paupstream_o(paupfile);
    paupstream_o << pauptree << endl;
    paupstream_o.close();


  if (!nostar) {

    ASTRALCladeExtractor ce(gtreefile, paupfile);


    unordered_set<Clade> newclades = ce.extract(ts);

    conf.add_clades(newclades.begin(), newclades.end());

    INFO << "Added " << newclades.size() << " clades" << endl;

    DEBUG << conf.get_clades().size() << endl;
  }
    Timer::stop("GetAdditionalClades");


}

void SVDQuestTripartitionScorer::setup(Config& conf, vector<Clade>& clades)
{

  if (doRunPaup) {
    runPaup(conf);
  }

  Timer::start("MakeBSDict");

  for (const Clade& clade : conf.get_clades()) {
    vector<Taxon> nonmembers;
    for(size_t i = 0; i < ts.size(); i++) {
      if (!clade.contains(i))
  	nonmembers.push_back(i);
    }

    map<pair<Taxon, Taxon>, double>& mp = W[clade.get_taxa()];

    if (clade.size() < 2) {
      continue;
    }

    for (size_t i = 0; i < nonmembers.size(); i++) {
      for (size_t j = i+1; j < nonmembers.size(); j++) {

  	double d = 0;
  	for (Taxon k : clade) {
  	  for (Taxon l : clade) {
  	    if (k > l)
  	      d += (*qd)(nonmembers[i],nonmembers[j],k,l);
  	  }
  	}

  	mp[make_pair(nonmembers[i], nonmembers[j])] = d;
  	mp[make_pair(nonmembers[j], nonmembers[i])] = d;
      }
    }
  }
  Timer::stop("MakeBSDict");





}

double SVDQuestTripartitionScorer::score(const Tripartition& t) {

  double val = 0;


  for (Taxon c : t.a2)
    for (Taxon d : t.rest)
      val +=  W[t.a1.get_taxa()][make_pair(c,d)];



  for (Taxon c : t.a1)
    for (Taxon d : t.rest)
      val +=  W[t.a2.get_taxa()][make_pair(c,d)];


  for (Taxon c : t.a2)
    for (Taxon d : t.a2)
      if (c > d)
	val += W[t.a1.get_taxa()][make_pair(c, d)];

  return val;
}
