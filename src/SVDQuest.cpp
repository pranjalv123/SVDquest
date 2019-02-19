#include <phylonaut/wASTRAL.hpp>
#include <phylonaut/DefaultTaxonSetExtractor.hpp>
#include <phylonaut/ASTRALCladeExtractor.hpp>
#include "SVDQuestTripartitionScorer.hpp"

#include <util/Logger.hpp>
#include <util/Timer.hpp>
#include <phylonaut/SingleTreeAnalysis.hpp>
#include <phylonaut/ConsensusTreeAnalysis.hpp>
#include <phylonaut/CountTreesAnalysis.hpp>
#include <phylonaut/ScoreAnalysis.hpp>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <cstdlib>
#if !(defined(_WIN32) || defined(_WIN64))
	#include <libgen.h>
#endif

using namespace std;


int main(int argc, char** argv) {
  vector<string> wastral_args;

  string input;
  string output = "";
  string extra = "";
  string alignment;
  int debug = 0;
  bool useDP = 0;
  bool wine=true;
  vector<string> output_labels;

  bool nostar = false;
  string scoreTree = "";


  bool getScore=false;
  bool getSingle=true;
  bool getGreedy=true;
  bool getMajority=true;
  bool getStrict=true;
  bool getAll=false;
  bool getCount=true;


  Logger::disable("DEBUG");
  Logger::enable("INFO");
  Logger::enable("PROGRESS");

  Config conf;
  conf.matrix=0;
  conf.analyses.push_back(new ScoreAnalysis());
  output_labels.push_back("score");
  string astralpath;
  for(int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-i" || string(argv[i]) == "--input") {
      assert(argc > i+1);
      i++;
      input = string(argv[i]);
    }
    if (string(argv[i]) == "--astral") {
      assert(argc > i+1);
      i++;
      astralpath = string(argv[i]);
    }
    if (string(argv[i]) == "--nowine") {
      wine=false;
    }
    if (string(argv[i]) == "-a" || string(argv[i]) == "--alignment") {
      assert(argc > i+1);
      i++;
      alignment = string(argv[i]);
    }
    if (string(argv[i]) == "-o"  || string(argv[i]) == "--output") {
      assert(argc > i+1);
      i++;
      fclose(fopen(argv[i], "a"));
#if (defined(_WIN32) || defined(_WIN64))
	  char buf[1000];
	  _fullpath(buf, argv[i], 1000);
	  output = buf;

#else
      output = string(realpath(argv[i], NULL));
#endif
    }
    if (string(argv[i]) == "-e" || string(argv[i]) == "--extra") {
      assert(argc > i+1);
      i++;
#if (defined(_WIN32) || defined(_WIN64))
	  char buf[1000];
	  _fullpath(buf, argv[i], 1000);
	  extra = buf;

#else
	  extra = string(realpath(argv[i], NULL));
#endif
    }
    if (string(argv[i]) == "--debug") {
      Logger::enable("DEBUG");
      Logger::enable("INFO");
      Logger::enable("PROGRESS");
      DEBUG << "Debug enabled\n";
    }
    if (string(argv[i]) == "--score") {
      getScore = true;
      assert(argc > i+1);
      i++;
      scoreTree = string(argv[i]);

    }
    if (string(argv[i]) == "--noscore") {
      getScore=false;
    }

    if (string(argv[i]) == "--nosingle") {
      getSingle=false;
    }
    if (string(argv[i]) == "--nogreedy") {
      getSingle=false;
    }
    if (string(argv[i]) == "--nomajority") {
      getMajority=false;
    }
    if (string(argv[i]) == "--nostrict") {
      getStrict=false;
    }
    if (string(argv[i]) == "--nocount") {
      getCount=false;
    }

    if (string(argv[i]) == "--score") {
      getScore=true;
    }
    if (string(argv[i]) == "--single") {
      getSingle=true;
    }
    if (string(argv[i]) == "--greedy") {
      getSingle=true;
    }
    if (string(argv[i]) == "--majority") {
      getMajority=true;
    }
    if (string(argv[i]) == "--strict") {
      getStrict=true;
    }
    if (string(argv[i]) == "--count") {
      getCount=true;
    }

    if (string(argv[i]) == "--profile") {
      conf.profile="profile";
    }
    if (string(argv[i]) == "--no-precompute") {
      useDP=1;
    }
    if (string(argv[i]) == "--no-star") {
      nostar=1;
    }




  }
  
  if (input.size() == 0 || output.size() == 0) {
    cerr << "SVDQuest -i <source tree file> -o <output file> [-e <extra trees>]" << endl;
    exit(-1);
  }


  if (getSingle) {
    conf.analyses.push_back(new SingleTreeAnalysis());
    output_labels.push_back("single");
  }
  if (getGreedy) {
    conf.analyses.push_back(new ConsensusTreeAnalysis(0.0));
    output_labels.push_back("greedy");
  }
  if (getMajority) {
    conf.analyses.push_back(new ConsensusTreeAnalysis(0.5));
    output_labels.push_back("majority");
  }
  if (getStrict) {
    conf.analyses.push_back(new ConsensusTreeAnalysis(1.0));
    output_labels.push_back("strict");
  }
  if (getAll) {
    // conf.analyses.push_back();
    // output_labels.push_back("all");
  }
  if (getCount) {
    conf.analyses.push_back(new CountTreesAnalysis());
    output_labels.push_back("count");
  }
  

  string quartetFile = output + ".svdQuartets";

  ifstream infile(input);
  int nlines = 0;
  string buffer;

  while(getline(infile, buffer)) {
    if (buffer.size() > 3) {
      nlines += 1;
    }
  }

  conf.scorer = new SVDQuestTripartitionScorer(alignment, output, astralpath, input);
  dynamic_cast<SVDQuestTripartitionScorer*>(conf.scorer)->nostar = nostar;
  dynamic_cast<SVDQuestTripartitionScorer*>(conf.scorer)->wine = wine;
  DEBUG << conf.scorer->clades_size() << endl;

  conf.taxon_extractor = new DefaultTaxonSetExtractor(input);



  if (getScore) {
    INFO << "Only one tree provided\nScoring tree instead of finding optimal tree" << endl;
    conf.extractors.push_back(new ASTRALCladeExtractor(scoreTree, "", false, true));
  } else {

    if (extra != "")
      conf.extractors.push_back(new ASTRALCladeExtractor(input, extra));
    else
      conf.extractors.push_back(new ASTRALCladeExtractor(input));
  }
  vector<string> trees = wASTRAL(conf);

  for (int i = 0; i < trees.size(); i++) {
    cout << output_labels.at(i) << endl;
    cout << trees.at(i) << endl;
    if (output.size() ) {
      ofstream outfile(output + "." + output_labels.at(i));
      outfile << trees.at(i) << endl;
      outfile.close();
    }
  }
  ofstream outfile(output + ".timing");
  Timer::writeAll(outfile);
  outfile.close();
  Timer::reset();

  Config conf_score;
  conf_score.matrix=0;

  if (!getScore) {

    conf_score.analyses.push_back(new ScoreAnalysis());
    conf_score.matrix = 0;

    string pauptreefile = output + ".pauptree";
    conf_score.taxon_extractor = conf.taxon_extractor;
    conf_score.extractors.push_back(new ASTRALCladeExtractor(pauptreefile, "", false, true));
    conf_score.scorer = new SVDQuestTripartitionScorer(*dynamic_cast<SVDQuestTripartitionScorer*>(conf.scorer), pauptreefile);

    trees = wASTRAL(conf_score);

    ofstream outfile2(output + ".pauptree_score");
    outfile2 << trees.at(0) << endl;
    outfile2.close();
    INFO << "PAUP score: " << trees.at(0) << endl;
  }
}
