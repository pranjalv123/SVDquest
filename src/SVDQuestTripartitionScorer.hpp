#ifndef SVDQUEST_TRIPARTITION_SCORER_HPP__
#define SVDQUEST_TRIPARTITION_SCORER_HPP__

#include "phylonaut/TripartitionScorer/TripartitionScorer.hpp"


class SVDQuestTripartitionScorer : public TripartitionScorer{
public:
  SVDQuestTripartitionScorer(string& alignment, string& outname, string& astralpath, string& gtreefile) : 
    alignmentfile(alignment), 
    outname(outname), 
    astralpath(astralpath), 
    gtreefile(gtreefile), 
    doRunPaup(true) {};
  
  SVDQuestTripartitionScorer(SVDQuestTripartitionScorer& other, string& gtreefile) : 
    gtreefile(gtreefile), 
    qd(other.qd), 
    doRunPaup(false) {};

  void runPaup(Config& conf);
  
  virtual void setup(Config& conf, vector<Clade>& clades);
  virtual double score(const Tripartition& t);

  bool wine;
  bool nostar;
private:
  unordered_map<clade_bitset, map<pair<Taxon, Taxon>, double> >  W;
  string alignmentfile, outname, astralpath, gtreefile;
  QuartetDict* qd;
  bool better(double newscore, double oldscore) {
    return newscore > oldscore;
  }
  bool doRunPaup;
};



#endif
