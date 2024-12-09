#pragma once

class Alg {
 private:
  /// _numV: number of nodes in the graph.
  uint32_t _numV;
  /// _numE: number of edges in the graph.
  size_t _numE;
  /// _numRRsets: number of RR sets.
  size_t _numRRsets = 0;
  /// Upper bound in the last round for __mode=1.
  double _boundLast = DBL_MAX;
  /// The minimum upper bound among all rounds for __model=2.
  double _boundMin = DBL_MAX;
  /// Two hyper-graphs, one is used for selecting seeds and the other is used for validating
  /// influence.
  THyperGraph _hyperGraph, _hyperGraphVldt;
  const Graph& _graph;

  iVector<offset_t> offsetCsr;
  iVector<nid_t> nbrCsr;
  iVector<ew_t> wCsr;

  /// Result object.
  TResult& _res;
  /// Seed set.
  Nodelist _vecSeed;
  ProbDist _probDist = WC;

  double _baseNumRRsets = 0.0;

  std::vector<rrid_t> _coverage;
  std::vector<bool> _isSeed;
  std::vector<std::vector<nid_t>> _vecThreshold;
  double _decayRate = 0.9;

  std::vector<uint32_t> _vecOutDegree;
  std::vector<uint32_t> _vecVldtInf;
  std::vector<std::pair<uint32_t, uint32_t>> _vecOrderedOutDeg;

  CascadeModel _cmodel;

  bool _sampleOpt = true;
  int _greedyType = GDY_GEAR;
  bool _decayGreedy = true;

  Graph _fGraph;
  std::vector<bool> _bfsBool;
  std::vector<dsfmt_t> _dsfmtSeedVec;

  std::vector<double> _vecCost;
  double _maxCost;
  double _minCost;
  double _totalBudget;
  double _basicCost;

  // variables for greedy strategy
  double _currCost = 0;
  double _lastCoverage = 0;
  double _sumUtils = 0.0;

  // variable for sample greedy
  bool sampleGreedyInit = false;

  // variable for gear
  double _maxEstRatio;
  std::vector<std::vector<uint32_t>> _greedyCandidateVec;
  int _currVecPos = -1;

  // variable for distorted greedy
  std::vector<uint32_t> _candidateVec;
  uint32_t _distortedTotalNum = 0;
  uint32_t _distortedIterId = 0;
  bool _distortedInit = false;

  /// Maximum coverage by lazy updating.
  double MaxCoverVanilla(const int targetSize);

  /// Maximum coverage by maintaining the top-k marginal coverage.
  double MaxCoverTopK(const int targetSize);
  /// Maximum coverage.
  double MaxCover(const int targetSize);

 public:
  Alg(const Graph& graph, TResult& tRes)
      : _graph(graph), _hyperGraph(graph), _hyperGraphVldt(graph), _res(tRes) {
    _numV = _hyperGraph.get_nodes();
    _numE = _hyperGraph.get_edges();
    _coverage.resize(_numV, 0);
    _isSeed.resize(_numV, false);

    log_info("graph info: n = %u, m=%lu", _numV, _numE);
#if 0
    createForwardGraph(graph);
#endif
    _dsfmtSeedVec.resize(THD_NUM);
    for (int tid = 0; tid < THD_NUM; tid++) {
      dsfmt_init_gen_rand(&_dsfmtSeedVec[tid], (uint32_t)rand());
    }
  }

  ~Alg() {}

  void convertToCsr() {
    offsetCsr.re_allocate(_numV);
    nbrCsr.re_allocate(_numE);
    wCsr.re_allocate(_numE);
    int64_t offset = 0;

    for (nid_t i = 0; i < _numV; i++) {
      nid_t nbrSize = _graph[i].size();
      offsetCsr.push_back({(int64_t)offset, (int64_t)(offset + nbrSize)});
      if (i == 0) {
        log_trace("offsetCsr: %ld, %ld", offsetCsr[0].first, offsetCsr[0].second);
      }

      assert(nbrSize >= 0);
      for (nid_t nbr = 0; nbr < nbrSize; nbr++) {
        nbrCsr[offset] = _graph[i][nbr].first;
        wCsr[offset] = _graph[i][nbr].second;

        offset++;
      }
      //   offset += nbrSize;
    }
    printf("m=%lu, offset: %ld\n", _numE, offset);
  }

  void setCascadeModel(CascadeModel model, const ProbDist dist) {
    log_info("cascade model: %s", model == IC ? "IC" : "LT");
    _cmodel = model;
    _hyperGraph.set_cascade_model(model);
    _hyperGraphVldt.set_cascade_model(model);

    _probDist = dist;
    _hyperGraph.set_prob_dist(dist);
    _hyperGraphVldt.set_prob_dist(dist);
  }

  void set_sample_type(bool sampleOpt, const bool isVanilla) {
    _sampleOpt = sampleOpt;
    if (isVanilla) {
      std::cout << "Vanilla sampling method is used" << std::endl;
    }

    _hyperGraph.set_vanilla_sample(isVanilla, _sampleOpt);
    _hyperGraphVldt.set_vanilla_sample(isVanilla, _sampleOpt);

    if (_sampleOpt) {
      convertToCsr();
      log_trace("offsetCsr Addr: %p, nbr addr: %p, wCsr addr: %p", offsetCsr.m_data, nbrCsr.m_data,
                wCsr.m_data);
      _hyperGraph.ctx.init(_cmodel, _probDist, _numV);
      _hyperGraphVldt.ctx.init(_cmodel, _probDist, _numV);
      _hyperGraph.ctx.setGraphAdj(offsetCsr.m_data, nbrCsr.m_data, wCsr.m_data);
      _hyperGraphVldt.ctx.setGraphAdj(offsetCsr.m_data, nbrCsr.m_data, wCsr.m_data);

      log_trace("offsetCsr Addr: %p, nbr addr: %p, wCsr addr: %p", _hyperGraph.ctx.offsetBase,
                _hyperGraph.ctx.adjBase, _hyperGraph.ctx.wBase);
    }
  }

  void set_greedy_type(int type, bool decay, double rate) {
    _greedyType = type;
    _decayGreedy = decay;
    _decayRate = rate;
  }

  void setCost(CostModelCtx& cmc) {
    // assert(_fGraph.size() != 0);
    _vecCost.resize(_numV);
    _basicCost = cmc.basicCost;

    _maxCost = cmc.basicCost - 1;
    _minCost = FLT_MAX;
    _totalBudget = 0;
    auto& rgraph = _hyperGraph._graph;
    size_t max_deg = 0;
    if (cmc.basicCost != 0.0) {
      for (uint32_t i = 0; i < _numV; i++) {
        double cost = cmc.basicCost + cmc.beta * rgraph[i].size();
        _vecCost[i] = cost;
        _totalBudget += cost;
        _maxCost = std::max(cost, _maxCost);
        _minCost = std::min(cost, _minCost);
      }
    }
    // log_info("max degree: %ld", max_deg);
    if (cmc.budget > 1.0) {
      _totalBudget = cmc.budget;
    } else {
      _totalBudget = cmc.budget * _totalBudget;
    }

    log_info("total budget: %f, max Cost: %f, basic cost: %f", _totalBudget, _maxCost, _basicCost);
  }

  void setRRsetNum(uint64_t rrnum) {
    _numRRsets = rrnum;
    log_info("rr set number: %lu", _numRRsets);
  }

  void createForwardGraph(const Graph& rGraph) {
    assert(_fGraph.size() == 0);
    _fGraph.resize(_numV);
    for (uint32_t t = 0; t < _numV; t++) {
      for (auto& s : rGraph[t]) {
        _fGraph[s.first].push_back(std::make_pair(t, s.second));
      }
    }

    _bfsBool.resize(_numV, false);
  }

  uint32_t oneSimulation(const std::vector<uint32_t>& seeds) {
    uint32_t numVisited = 0;
    std::vector<uint32_t> vecBfs;
    vecBfs.reserve(0.1 * _numV);
    for (auto seed : seeds) {
      _bfsBool[seed] = true;
      vecBfs.push_back(seed);
      numVisited++;
    }
    uint32_t currId = 0;

    while (currId < numVisited) {
      const uint32_t frontier = vecBfs[currId++];
      for (auto& elem : _fGraph[frontier]) {
        //   for (auto& elem : _graph[frontier]) {
        const uint32_t nbr = elem.first;
        if (_bfsBool[nbr]) continue;
        const auto randDouble = dsfmt_gv_genrand_open_close();

        if (randDouble > elem.second) continue;
        vecBfs.push_back(nbr);
        numVisited++;
        _bfsBool[nbr] = true;
      }
    }

    for (uint32_t i = 0; i < numVisited; i++) {
      _bfsBool[vecBfs[i]] = false;
    }

    return numVisited;
  }

  double calUtility(uint32_t influence) { return influence; }

  double runMonteCarloIC(uint32_t iterNum) {
    double utility = 0.0;

    for (uint32_t i = 0; i < iterNum; i++) {
      uint32_t inf = oneSimulation(_vecSeed);
      utility += calUtility(inf);
    }

    utility = utility / iterNum;

    double totalCost = 0.0;
    for (auto seed : _vecSeed) {
      totalCost += _vecCost[seed];
    }
    log_info("total cost: %f", totalCost);
    log_info("profit: %f", utility);

    return utility - totalCost;
  }

  void runMonteCarloOneThread(int tid, uint32_t sid, uint32_t eid, std::vector<double>& vecInf) {
    if (_cmodel == IC) {
      return runMonteCarloICOneThread(tid, sid, eid, vecInf);
    } else {
      return runMonteCarloLTOneThread(tid, sid, eid, vecInf);
    }
  }

  void runMonteCarloICOneThread(int tid, uint32_t sid, uint32_t eid, std::vector<double>& vecInf) {
    std::vector<uint32_t> vecBfs;
    vecBfs.reserve(0.1 * _numV);
    std::vector<bool> bfsBoolThread(_numV, false);
    // log_info("tid: %d, start: %u, end: %u", tid, sid, eid);

    const std::vector<uint32_t> seeds = _vecSeed;
    dsfmt_t* dsfmt = &_dsfmtSeedVec[tid];

    for (uint32_t sim = sid; sim < eid; sim++) {
      uint32_t numVisited = 0;
      vecBfs.resize(0);
      for (auto seed : seeds) {
        bfsBoolThread[seed] = true;
        vecBfs.push_back(seed);
        numVisited++;
      }
      uint32_t currId = 0;

      while (currId < numVisited) {
        const uint32_t frontier = vecBfs[currId++];
        for (auto& elem : _graph[frontier]) {
          const uint32_t nbr = elem.first;
          if (bfsBoolThread[nbr]) continue;
          const auto randDouble = dsfmt_genrand_open_close(dsfmt);

          if (randDouble > elem.second) continue;
          vecBfs.push_back(nbr);
          numVisited++;
          bfsBoolThread[nbr] = true;
        }
      }

      for (uint32_t i = 0; i < numVisited; i++) {
        bfsBoolThread[vecBfs[i]] = false;
      }
      vecInf[sim] = calUtility(numVisited);
    }
    return;
  }

  // weighted cascade
  void runMonteCarloLTOneThread(int tid, uint32_t sid, uint32_t eid, std::vector<double>& vecInf) {
    std::vector<int> Q;
    std::vector<bool> activated(_numV, false);
    std::vector<int> sampleEdge(_numV, -1);
    // log_info("tid: %d, start: %u, end: %u", tid, sid, eid);

    const std::vector<uint32_t> seeds = _vecSeed;
    dsfmt_t* dsfmt = &_dsfmtSeedVec[tid];

    for (uint32_t sim = sid; sim < eid; sim++) {
      uint32_t numVisited = 0;
      Q.resize(0);
      for (auto seed : seeds) {
        activated[seed] = true;
        Q.push_back(seed);
        const auto innbrs = _graph[seed].size();
        if (innbrs != 0) {
          sampleEdge[seed] = _graph[seed][dsfmt_genrand_uint32(dsfmt) % innbrs].first;
        }
        numVisited++;
      }
      uint32_t currId = 0;

      while (currId < numVisited) {
        const uint32_t frontier = Q[currId++];
        for (auto& elem : _fGraph[frontier]) {
          const uint32_t out_nbr = elem.first;
          if (activated[out_nbr]) continue;
          if (sampleEdge[out_nbr] < 0) {
            // in-nbrs must not be empty
            sampleEdge[out_nbr] =
                _graph[out_nbr][dsfmt_genrand_uint32(dsfmt) % _graph[out_nbr].size()].first;
          }
          if (sampleEdge[out_nbr] == frontier) {
            Q.push_back(out_nbr);
            numVisited++;
            activated[out_nbr] = true;
          }
        }
      }

      std::fill(activated.begin(), activated.end(), false);
      std::fill(sampleEdge.begin(), sampleEdge.end(), -1);
      vecInf[sim] = calUtility(numVisited);
      //   vecInf[sim] = numVisited;
      // log_info("numVisited: %u", numVisited);
    }
    return;
  }

  double getSeedTotalCost() {
    double totalCost = 0.0;
    for (auto seed : _vecSeed) {
      totalCost += _vecCost[seed];
    }
    return totalCost;
  }

  double runMonteCarloMultiThreads(uint32_t iterNum) {
    std::vector<std::future<void>> futures(THD_NUM);
    uint32_t avg_thread_num = iterNum / THD_NUM + 1;
    std::vector<double> vecInf(iterNum, 0);
    for (int tid = 0; tid < THD_NUM; tid++) {
      uint32_t startThd = tid * avg_thread_num;
      uint32_t endThd = startThd + avg_thread_num;
      endThd = std::min(iterNum, endThd);
      futures[tid] = std::async(std::launch::async, &Alg::runMonteCarloOneThread, this, tid,
                                startThd, endThd, std::ref(vecInf));
    }

    std::for_each(futures.begin(), futures.end(), std::mem_fn(&std::future<void>::wait));

    double util = std::accumulate(vecInf.begin(), vecInf.end(), 0.0) / ((double)iterNum);
    log_info("total exact utils: %f", util);

    return util;
  }

  double runRRsetVldt() {
    std::unordered_set<uint32_t> seedSet(_vecSeed.begin(), _vecSeed.end());
    double util = _hyperGraphVldt.EvalSeedSetInf(seedSet, 1024 * 10000);
    // log_info("use vanilla RRset get utility: %f", util);
    return util;
  }

  /// Set cascade model.
  void set_prob_dist(const ProbDist weight);
  void set_vanilla_sample(const bool isVanilla);

  void RefreshHypergraph() {
    _hyperGraph.RefreshHypergraph();
    _hyperGraphVldt.RefreshHypergraph();
  }
  /// Evaluate influence spread for the seed set constructed
  double EfficInfVldtAlg();
  /// Evaluate influence spread for a given seed set
  double EfficInfVldtAlg(const Nodelist vecSeed);

  double CalUtility(uint32_t inf) { assert(inf <= _numV); }

  double estimateRRSize();

  bool solveQuadratic(double a, double b, double c, double& largerRoot) {
    double discriminant = b * b - 4 * a * c;
    double sqrtDiscriminant = sqrt(abs(discriminant));

    if (discriminant < 0) return false;

    double root_1 = (-b + sqrtDiscriminant) / (2 * a);
    double root_2 = (-b - sqrtDiscriminant) / (2 * a);

    largerRoot = (root_1 > root_2) ? root_1 : root_2;
    log_info("root1: %f, root2: %f, largerRoot: %f", root_1, root_2, largerRoot);
    return true;
  }

  double dima(const double epsilon, const double delta) {
    Timer timerDima("DIMA");

    double gamma = 0.5 * (1.0 - 1.0 / exp(1));
    log_info("gamma: %f", gamma);

    const rrid_t numRbase = std::ceil((10000 * log(_totalBudget)) / BATCH_SIZE) * BATCH_SIZE;

    double time1 = 0, time2 = 0, time3 = 0;
    int idx = 1;
    rrid_t prevSize = 0;
    while (true) {
      const auto numR = numRbase << (idx - 1);
      std::cout << ": " << idx << " RR set: " << numR << std::endl;
      idx++;

      timerDima.get_operation_time();
      if (_sampleOpt) {
        _hyperGraph.BuildRRsetsWithPrefect(numR);
        _hyperGraphVldt.BuildRRsetsWithPrefect(numR);
      } else {
        _hyperGraph.BuildRRsets(numR);      // R1
                                            //   _hyperGraph.createInvertedList(prevSize, numR);
        _hyperGraphVldt.BuildRRsets(numR);  // R2
        //   _hyperGraphVldt.createInvertedList(prevSize, numR);
      }

      _numRRsets = _hyperGraph.get_RR_sets_size();
      log_info("complete build RRsets: %ld", _numRRsets);
      time1 += timerDima.get_operation_time();
      if (_decayGreedy) {
        log_info("decay");
        nodeSelectionThreshold();
      } else {
        log_info("no decay");
        nodeSelectionWithRRset();
      }
      time2 += timerDima.get_operation_time();

      log_info("seed size: %lu", _vecSeed.size());
      double utilR1 = _sumUtils / _numRRsets * _numV;
      double utilR2 = _hyperGraphVldt.CalculateInf(_vecSeed);
      time3 += timerDima.get_operation_time();

      log_info("R1 util: %f, R2 util: %f, cost: %f", utilR1, utilR2, _currCost);

      double t = (utilR1 - _currCost) / (utilR2 - _currCost);
      double rhs_epsilon_1 =
          1.0 * utilR2 / log(6.0 * idx * idx / delta) * ((double)_numRRsets / _numV);

      log_info("c-1: %f", rhs_epsilon_1);
      double eps_1 = 0;
      if (!solveQuadratic(rhs_epsilon_1 - 1, -3, -2, eps_1)) {
        log_info("epsilon 1 has not real roots or less than zero: %f", eps_1);
        continue;
      }

      double rhs_epsilon_2 = 1.0 * (utilR2 - (1.0 - eps_1) * _currCost) /
                             log(6.0 * idx * idx / delta) * ((double)_numRRsets / _numV);
      double eps_2 = 0;
      if (!solveQuadratic(rhs_epsilon_2, -2, -2, eps_2)) {
        log_info("epsilon 2 has not real roots or less than zero: %f", eps_2);
        continue;
      }

      log_info("t: %f, condition value: %f, %f", t, (t - 1) / t + eps_1 / gamma + eps_2,
               eps_1 / gamma + eps_2);
      if ((t - 1) / t + eps_1 / gamma + eps_2 < epsilon && eps_1 / gamma + eps_2 <= epsilon) {
        log_info("find satisfying seed set. total time: %f", timerDima.get_total_time());
        log_info("generate RRset time: %f, selection times: %f, evaluation: %f", time1, time2,
                 time3);
        return 0;
      }
    }

    return 0.0;
  }

  bool decayCostAwareGreedy(std::vector<double>& vecUtil, std::vector<bool>& isSeed, int& seed) {
    double maxRatio = -FLT_MAX;
    uint32_t maxSeed = _numV;
    double leftBudget = _totalBudget - _currCost;
    const double tmp1 = log(DECAY_FACTOR);

    if (leftBudget < _minCost) {
      // if (_totalBudget < _currCost + _vecCost[maxSeed]) {
      log_info("exceeding budget.");
      return false;
    }

    if (_currVecPos == -1) {
      _maxEstRatio = -FLT_MAX;
      _greedyCandidateVec.clear();
      _greedyCandidateVec.resize(DECAY_LIST_LEN + 1);

      for (uint32_t i = 0; i < _numV; i++) {
        if (_vecCost[i] > _totalBudget) continue;
        double ratio = (vecUtil[i] / (double)_numRRsets) / _vecCost[i];
        if (ratio > _maxEstRatio) {
          _maxEstRatio = ratio;
          maxSeed = i;
        }
      }

      const double tmp1 = log(DECAY_FACTOR);
      for (uint32_t i = 0; i < _numV; i++) {
        if (_vecCost[i] > _totalBudget) continue;
        double ratio = (vecUtil[i] / (double)_numRRsets) / _vecCost[i];
        if (ratio < 1e-8) {
          _greedyCandidateVec[DECAY_LIST_LEN].push_back(i);
          continue;
        }
        int64_t pos = std::ceil(log(ratio / _maxEstRatio) / tmp1);
        pos = (pos <= DECAY_LIST_LEN) ? pos : DECAY_LIST_LEN;
        _greedyCandidateVec[pos].push_back(i);
      }
      _currVecPos = 0;
      log_info("number of nodes at pos 0: %lu", _greedyCandidateVec[0].size());
    }

    while (_currVecPos < DECAY_LIST_LEN) {
      if (_greedyCandidateVec[_currVecPos].size() == 0) {
        _currVecPos++;
        continue;
      }

      bool foundNext = false;
      double threshold = _maxEstRatio * std::pow(DECAY_FACTOR, _currVecPos);
      // log_info("threshold: %f, currVecPos: %d, leftBudget: %f", threshold, _currVecPos,
      // leftBudget);
      auto& currVec = _greedyCandidateVec[_currVecPos];
      int lastElem = currVec.size();
      for (int i = lastElem - 1; i >= 0; i--) {
        auto node = currVec[i];
        if (_vecCost[node] > leftBudget) continue;
        double ratio = (vecUtil[node] / (double)_numRRsets) / _vecCost[node];
        if (ratio >= threshold) {
          foundNext = true;
          maxSeed = node;
          currVec.resize(i);
          break;
        } else {
          if (ratio < 1e-8) {
            _greedyCandidateVec[DECAY_LIST_LEN].push_back(node);
            continue;
          }
          int64_t pos = std::ceil(log(ratio / _maxEstRatio) / tmp1);
          pos = (pos <= DECAY_LIST_LEN) ? pos : DECAY_LIST_LEN;
          _greedyCandidateVec[pos].push_back(node);
        }
      }

      if (foundNext) {
        break;
      } else {
        _currVecPos++;
      }
    }

    if (_currVecPos == DECAY_LIST_LEN) {
      bool foundNext = false;
      auto& currVec = _greedyCandidateVec[DECAY_LIST_LEN];
      int lastElem = currVec.size();
      for (int i = lastElem - 1; i >= 0; i--) {
        auto node = currVec[i];
        if (_vecCost[node] > leftBudget) continue;
        maxSeed = node;
        foundNext = true;
        currVec.resize(i);
      }

      if (!foundNext) {
        maxSeed = _numV;
        currVec.resize(0);
      }
    }

    if (maxSeed == _numV) {
      // if (_totalBudget < _currCost + _vecCost[maxSeed]) {
      log_info("exceeding budget.");
      return false;
    }
    if (vecUtil[maxSeed] / (double)_numRRsets - _vecCost[maxSeed] < 0.0f) {
      log_info("reaching negative marginal gain.");
      return false;
    }

    seed = maxSeed;
    return true;
  }

  bool costAwareGreedy(std::vector<rrid_t>& vecUtil, std::vector<bool>& isSeed, int& seed) {
    double maxRatio = -FLT_MAX;
    uint32_t maxSeed = _numV;
    double leftBudget = _totalBudget - _currCost;
    for (uint32_t i = 0; i < _numV; i++) {
      if (isSeed[i]) continue;
      if (_vecCost[i] > leftBudget) continue;
      double ratio = ((double)_numV * vecUtil[i] / (double)_numRRsets) / _vecCost[i];
      if (ratio > maxRatio) {
        maxRatio = ratio;
        maxSeed = i;
      }
    }

    if (maxSeed == _numV) {
      // if (_totalBudget < _currCost + _vecCost[maxSeed]) {
      log_info("exceeding budget.");
      return false;
    }
    if ((double)_numV * vecUtil[maxSeed] / (double)_numRRsets - _vecCost[maxSeed] < 0.0f) {
      log_info("reaching negative marginal gain.");
      return false;
    }

    // log_info("max coverage: %ld", vecUtil[maxSeed]);
    seed = maxSeed;
    return true;
  }

  bool sampleGreedy(std::vector<rrid_t>& vecUtil, std::vector<bool>& isSeed, int& seed) {
    double maxRatio = -FLT_MAX;
    uint32_t maxSeed = _numV;
    double leftBudget = _totalBudget - _currCost;

    if (!sampleGreedyInit) {
      const double prob = 0.4142;
      for (uint32_t i = 0; i < _numV; i++) {
        double r = dsfmt_gv_genrand_open_close();
        if (r > prob) {
          isSeed[i] = true;
        }
      }
      sampleGreedyInit = true;
      log_info("sample greedy init.");
    }

    for (uint32_t i = 0; i < _numV; i++) {
      if (isSeed[i]) continue;
      if (_vecCost[i] > leftBudget) continue;
      double profit = ((double)_numV * vecUtil[i] / (double)_numRRsets) - _vecCost[i];
      double ratio = profit / _vecCost[i];
      if (ratio > maxRatio) {
        maxRatio = ratio;
        maxSeed = i;
      }
    }

    if (maxSeed == _numV || maxRatio <= 0) {
      // if (_totalBudget < _currCost + _vecCost[maxSeed]) {
      log_info("exceeding budget.");
      return false;
    }

    seed = maxSeed;
    return true;
  }

  bool simpleGreedy(std::vector<rrid_t>& vecUtil, std::vector<bool>& isSeed, int& seed) {
    double maxGain = -FLT_MAX;
    uint32_t maxSeed = _numV;
    double leftBudget = _totalBudget - _currCost;
    for (uint32_t i = 0; i < _numV; i++) {
      if (isSeed[i]) continue;
      if (_vecCost[i] > leftBudget) continue;
      double gain = ((double)_numV * vecUtil[i] / (double)_numRRsets) - _vecCost[i];
      if (gain > maxGain) {
        maxGain = gain;
        maxSeed = i;
      }
    }

    if (maxSeed == _numV) {
      // if (_totalBudget < _currCost + _vecCost[maxSeed]) {
      log_info("exceeding budget.");
      return false;
    }

    if (maxGain < 0.0f) {
      log_info("Greedy reaches negative marginal gain.");
      return false;
    }

    seed = maxSeed;
    return true;
  }

  bool distortedGreedy(std::vector<rrid_t>& vecUtil, std::vector<bool>& isSeed, int& seed) {
    double maxGain = -FLT_MAX;
    uint32_t maxSeed = _numV;
    double leftBudget = _totalBudget - _currCost;
    uint32_t seedNum = 0;  // include the nodes larger than leftBudget

    if (!_distortedInit) {
      _distortedInit = true;
      for (int i = 0; i < _numV; i++) {
        if (_vecCost[i] > _totalBudget) {
          _distortedIterId++;
          isSeed[i] = true;
        }
      }
    }

    if (_distortedIterId == _numV) {
      log_info("left budget is too small");
      return false;
    }

    // int maxCandidate = leftBudget / _minCost;
    const double base = 1 - 1.0 / (double)_numV;

    double distortedFactor = std::pow(base, _numV - _distortedIterId);
    log_info("distorted id: %u, factor: %f", _distortedIterId, distortedFactor);
    for (nid_t node = 0; node < _numV; node++) {
      if (isSeed[node]) continue;

      double gain =
          distortedFactor * ((double)_numV * vecUtil[node] / (double)_numRRsets) - _vecCost[node];
      if (gain > maxGain) {
        maxGain = gain;
        maxSeed = node;
      }
    }

    if (maxSeed == _numV) {
      log_info("budget is used up");
      return false;
    }

    if (maxGain < 0.0f) {
      log_info("profit is negative");
      return false;
    }

    _distortedIterId++;
    leftBudget -= _vecCost[maxSeed];
    for (nid_t node = 0; node < _numV; node++) {
      if (isSeed[node]) continue;
      if (_vecCost[node] > leftBudget) {
        _distortedIterId++;
        isSeed[node] = true;
      }
    }

    seed = maxSeed;
    return true;
  }

  bool costAwareVanillaGreedy(std::vector<double>& vecUtil, std::vector<bool>& isSeed, int& seed) {
    double maxRatio = -FLT_MAX;
    uint32_t maxSeed = _numV;
    double leftBudget = _totalBudget - _currCost;
    for (uint32_t i = 0; i < _numV; i++) {
      if (isSeed[i]) continue;
      if (_vecCost[i] > leftBudget) continue;
      double ratio = (double)_numV * (double)vecUtil[i] / (double)_numRRsets / _vecCost[i];
      if (maxRatio < ratio) {
        maxRatio = ratio;
        maxSeed = i;
      }
    }

    if (maxSeed == _numV) {
      // if (_totalBudget < _currCost + _vecCost[maxSeed]) {
      log_info("exceeding budget.");
      return false;
    }

    if ((double)_numV * (double)vecUtil[maxSeed] / (double)_numRRsets - _vecCost[maxSeed] < 0.0f) {
      log_info("reaching negative marginal gain");
      return false;
    }

    seed = maxSeed;
    return true;
  }

  bool greedySelectNode(std::vector<rrid_t>& vecUtil, std::vector<bool>& isSeed, int& seed) {
    bool cont = false;
    switch (_greedyType) {
      case GDY_GEAR:
        cont = costAwareGreedy(vecUtil, isSeed, seed);
        break;
      case GDY_SIMPLE:
        cont = simpleGreedy(vecUtil, isSeed, seed);
        break;
      case GDY_DISTORTED:
        cont = distortedGreedy(vecUtil, isSeed, seed);
        break;

      case GDY_SAMPLE:
        cont = sampleGreedy(vecUtil, isSeed, seed);
        break;

      default:
        break;
    }

    return cont;
  }

  rrid_t getCoverageAndMaxDeg(const HyperGraph& hg, std::vector<rrid_t>& coverage) {
    int64_t totalCoverage = 0;
    rrid_t maxDeg = 0;

    if (_sampleOpt) {
      for (nid_t i = 0; i < _numV; i++) {
        coverage[i] = hg._myInvertedList[i].m_num;
        maxDeg = std::max(maxDeg, coverage[i]);
        totalCoverage += coverage[i];
      }
    } else {
      for (nid_t i = 0; i < _numV; i++) {
        coverage[i] = hg._FRsets[i].size();
        maxDeg = std::max(maxDeg, coverage[i]);
        totalCoverage += coverage[i];
      }
    }
    // log_info("average rr size: %f", (double)totalCoverage / _numRRsets);
    // log_info("max ratio: %f", maxDeg);
    return maxDeg;
  }

  double getCoverageAndMaxRatio(const HyperGraph& hg, std::vector<rrid_t>& coverage) {
    int64_t totalCoverage = 0;
    double maxRatio = -1.0;

    if (_sampleOpt) {
      for (nid_t i = 0; i < _numV; i++) {
        coverage[i] = hg._myInvertedList[i].m_num;
        maxRatio = std::max(maxRatio, coverage[i] / _vecCost[i]);
        totalCoverage += coverage[i];
      }
    } else {
      for (nid_t i = 0; i < _numV; i++) {
        coverage[i] = hg._FRsets[i].size();
        maxRatio = std::max(maxRatio, coverage[i] / _vecCost[i]);
        totalCoverage += coverage[i];
      }
    }
    log_info("average rr size: %f", (double)totalCoverage / _numRRsets);
    // log_info("max ratio: %f", maxRatio);
    return maxRatio * _numV / _numRRsets;
  }

  void updateCoverageInfo(HyperGraph& hg, std::vector<bool>& edgeMark,
                          std::vector<rrid_t>& coverage, const nid_t seedNode) {
    rrid_t originalDeg;
    rrid_t* startAddr;
    if (_sampleOpt) {
      originalDeg = hg._myInvertedList[seedNode].m_num;
      startAddr = hg._myInvertedList[seedNode].m_data;
      for (rrid_t i = 0; i < originalDeg; i++) {
        const rrid_t edgeIdx = startAddr[i];
        if (edgeMark[edgeIdx]) continue;

        edgeMark[edgeIdx] = true;

        // const auto& rraddr = hg.ctx.RRs[edgeIdx].addr;
        // const nid_t s = hg.ctx.RRs[edgeIdx].size;

        const auto& rraddr = hg.ctx.optRRs[edgeIdx]->m_data;
        const nid_t s = hg.ctx.optRRs[edgeIdx]->m_num;

        for (nid_t curr = 0; curr < s; curr++) {
          const auto& nodeIdx = rraddr[curr];
          if (coverage[nodeIdx] == 0) continue;  // This node is seed, skip
          coverage[nodeIdx]--;
        }
      }

    } else {
      originalDeg = hg._FRsets[seedNode].size();
      startAddr = &(hg._FRsets[seedNode][0]);

      for (rrid_t i = 0; i < originalDeg; i++) {
        const rrid_t edgeIdx = startAddr[i];
        if (edgeMark[edgeIdx]) continue;

        edgeMark[edgeIdx] = true;

        for (auto nodeIdx : hg._RRsets[edgeIdx]) {
          if (coverage[nodeIdx] == 0) continue;  // This node is seed, skip
          coverage[nodeIdx]--;
        }
      }
    }
  }

  void nodeSelectionWithRRset() {
    getCoverageAndMaxRatio(_hyperGraph, _coverage);
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();
    _currCost = 0;
    _sumUtils = 0;

    while (true) {
      int nextSeed = -1;
      if (!greedySelectNode(_coverage, _isSeed, nextSeed)) {
        log_info("size of seed set : %lu, total cost: %f", _vecSeed.size(), getSeedTotalCost());

        std::fill(_isSeed.begin(), _isSeed.end(), false);
        return;
      }
      _currCost += _vecCost[nextSeed];

      _vecSeed.push_back(nextSeed);
      _isSeed[nextSeed] = true;
      _sumUtils += _coverage[nextSeed];

      updateCoverageInfo(_hyperGraph, edgeMark, _coverage, nextSeed);

      assert(_coverage[nextSeed] == 0);
    }
  }

  void nodeSelectionThreshold() {
    switch (_greedyType) {
      case GDY_GEAR:
        nodeSelectionCAThreshold();
        break;
      case GDY_SIMPLE:
        log_info("it is not implemented");
        exit(1);
        break;
      case GDY_DISTORTED:
        log_info("it is not implemented");
        break;

      case GDY_SAMPLE:
        nodeSelectionSampleThreshold();
        break;

      default:
        break;
    }
  }

  void nodeSelectionCAThreshold() {
    const double maxRatio = getCoverageAndMaxRatio(_hyperGraph, _coverage);
    _vecThreshold.clear();
    double minRatio = (double)_numV / ((double)_numRRsets * _maxCost);
    const nid_t maxIndex = std::ceil(std::log(minRatio / maxRatio) / std::log(_decayRate));
    _vecThreshold.resize(maxIndex + 1);

    const double c = (double)_numV / (double)_numRRsets;
    for (nid_t i = 0; i < _numV; i++) {
      if (_vecCost[i] > _totalBudget) continue;

      if (_coverage[i] == 0) {
        _vecThreshold[maxIndex].push_back(i);
      } else {
        double currRatio = _coverage[i] * c / _vecCost[i];
        nid_t index = std::ceil(std::log(currRatio / maxRatio) / std::log(_decayRate));
        _vecThreshold[index].push_back(i);
      }
    }

    log_info("size of vecThreshold[0]: %lu", _vecThreshold[0].size());

    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();
    _currCost = 0;
    _sumUtils = 0;
    double leftBudget = _totalBudget;

    // decayRate^power
    for (nid_t p = 0; p < maxIndex; p++) {
      auto& vecRatio = _vecThreshold[p];

      const nid_t len = vecRatio.size();
      for (nid_t idx = len; idx--;) {
        auto node = vecRatio[idx];
        if (_vecCost[node] >= leftBudget) continue;

        if (_coverage[node] == 0) continue;

        double currRatio = _coverage[node] * c / _vecCost[node];
        nid_t currPower = std::ceil(std::log(currRatio / maxRatio) / std::log(_decayRate));
        if (currPower > p) {
          _vecThreshold[currPower].push_back(node);
          continue;
        }
        // log_info("ratio: %f, curr deg: %u", _coverage[node] * c/_vecCost[node], currPower);

        if (currRatio - 1 <= 0) {
          log_info("exceeding budget.");
          goto finished;
        }

        _currCost += _vecCost[node];
        leftBudget = _totalBudget - _currCost;
        _vecSeed.push_back(node);
        _isSeed[node] = true;
        _sumUtils += _coverage[node];

        updateCoverageInfo(_hyperGraph, edgeMark, _coverage, node);
        assert(_coverage[node] == 0);

        if (leftBudget < _basicCost) {
          log_info("budget left is smaller than the mininum cost");
          goto finished;
        }
      }
      vecRatio.clear();
    }
    log_info("reach coverage 0, basicCost: %f", _basicCost);

  finished:
    for (auto node : _vecSeed) {
      _isSeed[node] = false;
    }
    return;
  }

  void nodeSelectionSampleThreshold() {
    const double maxRatio = getCoverageAndMaxRatio(_hyperGraph, _coverage);
    _vecThreshold.clear();
    double minRatio = (double)_numV / ((double)_numRRsets * _maxCost);
    const nid_t maxIndex = std::ceil(std::log(minRatio / maxRatio) / std::log(_decayRate));
    _vecThreshold.resize(maxIndex + 1);

    if (!sampleGreedyInit) {
      const double prob = 0.4142;
      for (uint32_t i = 0; i < _numV; i++) {
        double r = dsfmt_gv_genrand_open_close();
        if (r > prob) {
          _isSeed[i] = true;
        }
      }
      sampleGreedyInit = true;
      log_info("sample greedy init.");
    }

    const double c = (double)_numV / (double)_numRRsets;
    for (nid_t i = 0; i < _numV; i++) {
      if (_isSeed[i] || _vecCost[i] > _totalBudget) continue;

      if (_coverage[i] == 0) {
        _vecThreshold[maxIndex].push_back(i);
      } else {
        double currRatio = _coverage[i] * c / _vecCost[i];
        nid_t index = std::ceil(std::log(currRatio / maxRatio) / std::log(_decayRate));
        _vecThreshold[index].push_back(i);
      }
    }

    log_info("size of vecThreshold[0]: %lu", _vecThreshold[0].size());

    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();
    _currCost = 0;
    _sumUtils = 0;
    double leftBudget = _totalBudget;

    // decayRate^power
    for (nid_t p = 0; p < maxIndex; p++) {
      auto& vecRatio = _vecThreshold[p];

      const nid_t len = vecRatio.size();
      for (nid_t idx = len; idx--;) {
        auto node = vecRatio[idx];
        if (_isSeed[node] || _vecCost[node] >= leftBudget) continue;

        if (_coverage[node] == 0) continue;

        double currRatio = _coverage[node] * c / _vecCost[node];
        nid_t currPower = std::ceil(std::log(currRatio / maxRatio) / std::log(_decayRate));
        if (currPower > p) {
          _vecThreshold[currPower].push_back(node);
          continue;
        }
        // log_info("ratio: %f, curr deg: %u", _coverage[node] * c/_vecCost[node], currPower);

        if (currRatio - 1 <= 0) {
          log_info("exceeding budget.");
          goto finished;
        }

        _currCost += _vecCost[node];
        leftBudget = _totalBudget - _currCost;
        _vecSeed.push_back(node);
        _isSeed[node] = true;
        _sumUtils += _coverage[node];

        updateCoverageInfo(_hyperGraph, edgeMark, _coverage, node);
        assert(_coverage[node] == 0);

        if (leftBudget < _basicCost) {
          log_info("budget left is smaller than the mininum cost");
          goto finished;
        }
      }
      vecRatio.clear();
    }
    log_info("reach coverage 0, basicCost: %f", _basicCost);

  finished:
    for (auto node : _vecSeed) {
      _isSeed[node] = false;
    }
    return;
  }

  double evalGreedyWithVanillaRRset() {
    Timer timerDIMA("DIMA");
#ifdef VTUNE_PROF
    __itt_resume();
#endif

    double totalRRsetTime = 0;
    _hyperGraph.BuildRRsets(_numRRsets);

    double t1 = timerDIMA.get_operation_time();
    log_info("time for construct RRsets: %f", t1);
    totalRRsetTime += t1;

#ifdef VTUNE_PROF
    __itt_pause();
#endif

    nodeSelectionWithRRset();
    log_info("selection time: %f", timerDIMA.get_operation_time());
    log_info("total running time: %f", timerDIMA.get_total_time());
    return 0;
  }

  void evalGreedy() {
    log_info("sample optimized: %s, greedyType: %d", _sampleOpt ? "yes" : "no", _greedyType);

    evalGreedyWithVanillaRRset();

    double totalCost = getSeedTotalCost();
    //   log_info("forward profit for the selected seed: %f",
    //    runMonteCarloMultiThreads(DEFAULT_FWD_NUM) - totalCost);
    log_info("profit for the selected seed: %f", runRRsetVldt() - totalCost);
  }

  void evalRRset() {
    Timer timerEvalRR("DIMA");

#ifdef VTUNE_PROF
    __itt_resume();
#endif

    _hyperGraph.BuildRRsets(_numRRsets);

#ifdef VTUNE_PROF
    __itt_pause();
#endif

    double usedTime = timerEvalRR.get_operation_time();
    log_info("build RR set time: %f, each one: %f", usedTime, usedTime / _numRRsets);
    _hyperGraph.getRRsetAvgSize(_numRRsets);
  }

  double subsimOnly(const int targetSize, const double epsilon, const double delta);
};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;