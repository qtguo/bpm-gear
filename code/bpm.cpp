#include "stdafx.h"
#include "log.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"



void init_random_seed() {
  srand(time(NULL) + getpid());
  // Randomize the seed for generating random numbers
  dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
}

void output_command(int argc, char* argv[]) {
  printf("The command line is: ");

  for (int i = 0; i < argc; i++) {
    printf("%s ", argv[i]);
  }

  printf("\n");
}

int main(int argc, char* argv[]) {
#ifdef VTUNE_PROF
  __itt_pause();
#endif

  output_command(argc, argv);

  TArgument Arg(argc, argv);

  if (Arg._probDist == PROB_DIST_ERROR) {
    log_info("The input probability distribution is not supported: %s", Arg._probDistStr.c_str());
    log_info("The supported probability distribution: weights, wc, uniform, skewed");
    return 0;
  }

  if (Arg._func == FUNC_ERROR) {
    log_info("The input func is not supported: %s", Arg._funcStr.c_str());
    log_info("The supported func: format, im");
  }

  init_random_seed();

  const std::string infilename = Arg._dir + "/" + Arg._graphname;
  if (Arg._func == FORMAT) {
    // Format the graph
    GraphBase::FormatGraph(infilename, Arg._probDist, Arg._wcVar, Arg._probEdge, Arg._skewType);
    return 0;
  }

  Timer mainTimer("main");
  // Load the reverse graph
  Graph graph;
  GraphBase::LoadGraph(graph, infilename);

  int probDist = GraphBase::LoadGraphProbDist(infilename);

  // Initialize a result object to record the results
  TResult tRes;
  TAlg tAlg(graph, tRes);
  tAlg.setCascadeModel(Arg._cascadeModel, (ProbDist)probDist);
  // tAlg.set_prob_dist((ProbDist)probDist);  // Set propagation model
  tAlg.set_sample_type(Arg._sampleOpt, Arg._vanilla);
  tAlg.set_greedy_type(Arg._greedyType, Arg._decayGreedy, Arg._decayRate);

  tAlg.setCost(Arg._cmc);
  if (Arg._func == PM) {
    log_info("epsilon: %f", Arg._eps);
    tAlg.dima(Arg._eps, 1.0 / graph.size());
    double totalCost = tAlg.getSeedTotalCost();
    log_info("profit for the selected seed: %f", tAlg.runRRsetVldt() - totalCost);
  }

  log_info("BATCH SIZE: %d", BATCH_SIZE);

  if (Arg._func == EVAL_GDY) {
    tAlg.setRRsetNum(Arg._rrnum);
    tAlg.evalGreedy();
  }

  if (Arg._func == IM) {
    log_info("eps: %f", Arg._eps);
    tAlg.subsimOnly(Arg._seedsize, Arg._eps, 1.0 / graph.size());
  }

  if (Arg._func == EVAL_RRSET) {
    tAlg.setRRsetNum(Arg._rrnum);
    tAlg.evalRRset();
  }
  tAlg.RefreshHypergraph();
  return 0;
}