#pragma once

class Argument {
 public:
  // Function parameter.
  // format: format graph
  // im: influence maximization
  std::string _funcStr = "im";
  FuncType _func = EVAL_GDY;

  // The number of nodes to be selected. Default is 50.
  int _seedsize = 50;

  // For the uniform setting, every edge has the same diffusion probability.
  float _probEdge = float(0.1);

  // Error threshold 1-1/e-epsilon.
  double _eps = 0.1;

  // Failure probability delta. Default is 1/#nodes.
  double _delta = -1.0;

  CostModelCtx _cmc;
  bool _sampleOpt = true;
  double _scaleFactor = 1.0;
  std::string _greedyTypeStr = GDY_STR_GEAR;

  int64_t _rrnum = DEFAULT_RRSET_NUM;
  // Graph name. Default is "facebook".
  std::string _graphname = "facebook";

  // Probability distribution
  // weights: graph data with weights
  // wc: wc setting
  // uniform: uniform setting
  // skewed: skewed distribution
  std::string _probDistStr = "wc";
  ProbDist _probDist = WC;
  std::string _skewType = "exp";
  CascadeModel _cascadeModel = IC;

  // Directory
  std::string _dir = "graphInfo";

  // Result folder
  std::string _resultFolder = "result";

  // File name of the result
  std::string _outFileName;

  // wc variant
  double _wcVar = 1.0;

  // sample RR set with the vanilla method
  bool _vanilla = false;

  // use hist algorithm
  int _greedyType = GDY_GEAR;
  bool _decayGreedy = true;
  double _decayRate = 0.9;

  Argument(int argc, char* argv[]) {
    std::string param, value;

    for (int ind = 1; ind < argc; ind++) {
      if (argv[ind][0] != '-') break;

      std::stringstream sstr(argv[ind]);
      getline(sstr, param, '=');
      getline(sstr, value, '=');

      if (!param.compare("-func"))
        _funcStr = value;
      else if (!param.compare("-seedsize"))
        _seedsize = stoi(value);
      else if (!param.compare("-eps"))
        _eps = stod(value);
      else if (!param.compare("-delta"))
        _delta = stod(value);
      else if (!param.compare("-gname"))
        _graphname = value;
      else if (!param.compare("-dir"))
        _dir = value;
      else if (!param.compare("-outpath"))
        _resultFolder = value;
      else if (!param.compare("-pdist"))
        _probDistStr = value;
      else if (!param.compare("-model"))
        _cascadeModel = (CascadeModel)(value == "lt");
      else if (!param.compare("-pedge"))
        _probEdge = stof(value);
      else if (!param.compare("-wcvariant"))
        _wcVar = stod(value);
      else if (!param.compare("-skew"))
        _skewType = value;
      else if (!param.compare("-vanilla"))
        _vanilla = (value == "1");
      else if (!param.compare("-sampleOpt"))
        _sampleOpt = (value == "1");
      else if (!param.compare("-decayGreedy"))
        _decayGreedy = (value == "1");
      else if (!param.compare("-decayRate"))
        _decayRate = stof(value);

      else if (!param.compare("-greedyType"))
        _greedyTypeStr = value;
      // cost model
      else if (!param.compare("-randomCost"))
        _cmc.randomCost = (value == "1");
      else if (!param.compare("-basicCost"))
        _cmc.basicCost = stof(value);
      else if (!param.compare("-beta"))
        _cmc.beta = stof(value);
      // budget allocation
      else if (!param.compare("-budget"))
        _cmc.budget = stof(value);
      // RR sets number
      else if (!param.compare("-rrnum"))
        _rrnum = stoll(value);
    }

    if (_wcVar <= 0) {
      // wrong input
      _wcVar = 1.0;
    }

    decode_func_type();
    decode_greedy_type();
    decode_prob_dist();
  }

  void build_outfilename(int seedSize, ProbDist dist, Graph& graph) {
    std::string distStr;

    if (dist == WEIGHTS) {
      _probDistStr = "weights";
    } else if (dist == WC) {
      _probDistStr = "wc";
    } else if (dist == UNIFORM) {
      _probDistStr = "uniform";

      for (int i = 0; i < graph.size(); i++) {
        if (graph[i].size() > 0) {
          _probEdge = graph[i][0].second;
          break;
        }
      }
    } else {
      _probDistStr = "skewed";
    }

    _outFileName = TIO::BuildOutFileName(_graphname, "subsim", seedSize, _probDistStr, _probEdge);

    return;
  }

  void decode_prob_dist() {
    if (_probDistStr == "wc") {
      _probDist = WC;
    } else if (_probDistStr == "uniform") {
      _probDist = UNIFORM;
    } else if (_probDistStr == "skewed") {
      _probDist = SKEWED;
    } else if (_probDistStr == "weights") {
      _probDist = WEIGHTS;
    } else {
      _probDist = PROB_DIST_ERROR;
    }
  }

  void decode_func_type() {
    if (_funcStr == "format") {
      _func = FORMAT;
    } else if (_funcStr == "pm") {
      _func = PM;
    } else if (_funcStr == "evalGreedy") {
      _func = EVAL_GDY;
    } else if (_funcStr == "im") {
      _func = IM;
    } else if (_funcStr == "evalRRset") {
      _func = EVAL_RRSET;
    } else {
      _func = FUNC_ERROR;
    }
  }

  void decode_greedy_type() {
    if (_greedyTypeStr == GDY_STR_SIMPLE) {
      _greedyType = GDY_SIMPLE;
    } else if (_greedyTypeStr == GDY_STR_DISTORTED) {
      _greedyType = GDY_DISTORTED;
    } else if (_greedyTypeStr == GDY_STR_GEAR) {
      _greedyType = GDY_GEAR;
    } else if (_greedyTypeStr == GDY_STR_SAMPLE) {
      _greedyType = GDY_SAMPLE;
    } else {
      _greedyType = GDY_ERROR;
    }
  }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;
