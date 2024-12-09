#ifndef _CACHEOPT_H_
#define _CACHEOPT_H_

// #include "SFMT/dSFMT/dSFMT.h"

#define PREFETCH_HINT _MM_HINT_T0
// #define _mm_prefetch(x, y)

enum SUBSIM_STATUS { FRONTIER_WORKING, FRONTIER_FINISHED, LT_WORKING, RRSET_COMPLETED };
enum CASCADE_MODEL { MODEL_IC, MODEL_LT };

#ifndef BATCH_SIZE
#define BATCH_SIZE 16
#endif

// rr context
struct spread {
  offset_t frontierOffset;
  nid_t frontier;
  nid_t frontierNodeId;
  nid_t nbrSize;

  // record the process of jumping action
  nid_t currNbr;
  nid_t currNbrId;
  ew_t jumpingProb;
  ew_t landingProb;
  // optimized WC
  iVector<nid_t>* rrPtr;
  iVector<nid_t> wcTmp;
  const nid_t* ltPos;
  nid_t ltTarget;
};

struct batchElem {
  spread spd;
  bArray visitedFlag;
};

struct batchCtx {
  // the number of RR sets
  const offset_t* offsetBase;
  const ew_t* wBase;
  const nid_t* adjBase;
  nid_t graphSize;

  rrid_t numCompleted = 0;
  rrid_t numWorking = 0;
  int8_t spdModel;
  int8_t spdCas;

  batchElem elems[BATCH_SIZE];
  bid_t elemIsIdle[BATCH_SIZE];
  bid_t elemStatus[BATCH_SIZE];

  dsfmt_t sfmt;
  std::vector<iVector<nid_t>*> optRRs;

  void init(int8_t cascade, int8_t model, nid_t numGraph) {
    spdCas = cascade;
    spdModel = model;
    graphSize = numGraph;

    for (bid_t i = 0; i < BATCH_SIZE; i++) {
      elems[i].visitedFlag.makeSpace(graphSize);
      elemIsIdle[i] = 1;
      elemStatus[i] = RRSET_COMPLETED;
    }

    srand(time(NULL));
    dsfmt_init_gen_rand(&sfmt, rand());
  }

  void setGraphAdj(offset_t* offset, nid_t* adj, ew_t* weights) {
    offsetBase = offset;
    adjBase = adj;
    wBase = weights;
  }

  void batchICOpt_Execute(nid_t* targets, rrid_t numRRsets) {
    const nid_t maxLen = graphSize + 1;

    for (bid_t i = 0; i < BATCH_SIZE; i++) {
      auto& elem = elems[i];

      nid_t target = targets[numWorking];
      elem.spd.rrPtr = new iVector<nid_t>();
      elem.spd.rrPtr->push_back(target);
      elem.spd.frontier = 0;
      elem.visitedFlag.setTrue(target);
      elemIsIdle[i] = 0;
      elemStatus[i] = FRONTIER_WORKING;
      numWorking++;
    }
    while (numCompleted < numRRsets) {
      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        log_trace("step 1: i=%u, status: %u, idle: %u", i, elemStatus[i], elemIsIdle[i]);
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        auto& spd = elem.spd;
        spd.currNbr = 0;

        spd.frontierNodeId = elem.spd.rrPtr->m_data[spd.frontier];
        log_trace("i: %u, frontier: %u", i, spd.frontierNodeId);
        _mm_prefetch((void*)(offsetBase + spd.frontierNodeId), PREFETCH_HINT);
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        log_trace("step 2.1: i=%d, status: %d", i, elemStatus[i]);
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        auto& spd = elem.spd;
        spd.frontierOffset = offsetBase[spd.frontierNodeId];
        spd.nbrSize = spd.frontierOffset.second - spd.frontierOffset.first;
        log_trace("i: %d, first: %ld, second: %ld, size: %d", i, spd.frontierOffset.first,
                  spd.frontierOffset.second, spd.nbrSize);

        if (spd.nbrSize == 0) {
          elemStatus[i] = FRONTIER_FINISHED;
        }
        const ew_t* probAddr = wBase + spd.frontierOffset.first + spd.currNbr;
        _mm_prefetch((void*)(probAddr), PREFETCH_HINT);
      }

      bool moreWork = false;
      while (1) {
        moreWork = false;
        for (bid_t i = 0; i < BATCH_SIZE; i++) {
          log_trace("step 2.2: i=%d, status: %d", i, elemStatus[i]);
          auto& elem = elems[i];

          if (elemIsIdle[i]) continue;

          if (elemStatus[i] == FRONTIER_WORKING) {
            moreWork = true;

            auto& spd = elem.spd;
            spd.jumpingProb = wBase[spd.frontierOffset.first + spd.currNbr];

            if (spd.jumpingProb < 1.0) {
              double log_p = std::log(1.0 - spd.jumpingProb);
              nid_t incr = truncSkipLen(std::log(dsfmt_genrand_open_close(&sfmt)) / log_p, maxLen);
              log_trace("increment: %d, prob: %f", incr, spd.jumpingProb);

              spd.currNbr += incr;
            } else {
            }

            int64_t os = spd.frontierOffset.first + spd.currNbr;
            if (os >= spd.frontierOffset.second) {
              elemStatus[i] = FRONTIER_FINISHED;
              continue;
            }

            _mm_prefetch((void*)(adjBase + os), PREFETCH_HINT);

            _mm_prefetch((void*)(wBase + os), PREFETCH_HINT);
          }
        }

        if (!moreWork) break;

        moreWork = false;
        for (bid_t i = 0; i < BATCH_SIZE; i++) {
          log_trace("step 2.3: i=%d, status: %d", i, elemStatus[i]);
          auto& elem = elems[i];
          if (elemIsIdle[i]) continue;

          if (elemStatus[i] == FRONTIER_WORKING) {
            moreWork = true;

            auto& spd = elem.spd;
            spd.currNbrId = adjBase[spd.frontierOffset.first + spd.currNbr];
            log_trace("currNbr: %d", spd.currNbr);

            spd.landingProb = wBase[spd.frontierOffset.first + spd.currNbr];

            _mm_prefetch((void*)(elem.visitedFlag.getData(spd.currNbrId)), PREFETCH_HINT);
          }
        }

        for (bid_t i = 0; i < BATCH_SIZE; i++) {
          log_trace("step 2.4: i=%d, status: %d", i, elemStatus[i]);
          auto& elem = elems[i];
          if (elemIsIdle[i]) continue;

          if (elemStatus[i] == FRONTIER_WORKING) {
            moreWork = true;

            auto& spd = elem.spd;
            spd.currNbrId = adjBase[spd.frontierOffset.first + spd.currNbr];

            double prob = dsfmt_genrand_open_close(&sfmt);
            bool accept = true;
            if (prob > spd.landingProb / spd.jumpingProb) {
              accept = false;
            }

            if (accept && !elem.visitedFlag.getBool(spd.currNbrId)) {
              elem.spd.rrPtr->push_back(spd.currNbrId);
              elem.visitedFlag.setTrue(spd.currNbrId);
            }

            spd.currNbr++;
            if (spd.currNbr >= spd.nbrSize) {
              elemStatus[i] = FRONTIER_FINISHED;
            } else {
              spd.jumpingProb = wBase[spd.frontierOffset.first + spd.currNbr];
            }
          }
        }
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        log_trace("step 3: i=%d, status: %d", i, elemStatus[i]);
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        assert(elemStatus[i] == FRONTIER_FINISHED);

        auto& spd = elem.spd;

        spd.frontier++;
        spd.currNbr = 0;

        if (spd.frontier < spd.rrPtr->m_num) {
          elemStatus[i] = FRONTIER_WORKING;
          continue;
        } else {
          log_trace("RR set complete: i=%d, numCompleted: %ld, numWorking: %ld", i, numCompleted,
                    numWorking);
          numCompleted++;
        }

        // restore the environment
        for (nid_t i = 0; i < spd.rrPtr->m_num; i++) {
          elem.visitedFlag.setFalse(spd.rrPtr->m_data[i]);
        }
        optRRs.push_back(spd.rrPtr);
        spd.rrPtr = new iVector<nid_t>();

        // new sampling
        if (numWorking < numRRsets) {
          elemStatus[i] = FRONTIER_WORKING;
          elem.spd.frontier = 0;
          int t = targets[numWorking++];
          elem.spd.rrPtr->push_back(t);
          elem.visitedFlag.setTrue(t);
        } else {
          elemIsIdle[i] = true;
        }
      }
    }
  }

  void batchWCOpt_Execute(nid_t* targets, rrid_t numRRsets) {
    const nid_t maxJumpLen = graphSize + 1;
    numCompleted = 0;

    numWorking = 0;
    for (bid_t i = 0; i < BATCH_SIZE; i++) {
      auto& elem = elems[i];
      elem.spd.rrPtr = new iVector<nid_t>();

      nid_t target = targets[numWorking];
      elem.spd.rrPtr->push_back(target);
      elem.spd.frontier = 0;
      elem.visitedFlag.setTrue(target);

      elemIsIdle[i] = 0;
      elemStatus[i] = FRONTIER_WORKING;
      numWorking++;
    }

    while (numCompleted < numRRsets) {
      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        log_trace("step 1: i=%d, status: %d, idle: %d", i, elemStatus[i], elemIsIdle[i]);
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        auto& spd = elem.spd;

        spd.frontierNodeId = elem.spd.rrPtr->m_data[spd.frontier];
        log_trace("elem batch id: %d, frontier: %u, frontierNodeId: %u", i, spd.frontier,
                  spd.frontierNodeId);
        _mm_prefetch((void*)(offsetBase + spd.frontierNodeId), PREFETCH_HINT);
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        auto& spd = elem.spd;
        spd.frontierOffset = offsetBase[spd.frontierNodeId];
        nid_t nbrsNum = spd.frontierOffset.second - spd.frontierOffset.first;
        log_trace("i: %d, first: %ld, second: %ld, size: %d", i, spd.frontierOffset.first,
                  spd.frontierOffset.second, nbrsNum);

        elemStatus[i] = FRONTIER_WORKING;

        spd.wcTmp.clean();
        if (nbrsNum > 1) {
          double jumpingProb = (1.0 / (double)nbrsNum);
          const double log_p = std::log(1 - jumpingProb);
          nid_t startPos =
              truncSkipLen(std::log(dsfmt_genrand_open_close(&sfmt)) / log_p, maxJumpLen);

          while (startPos < nbrsNum) {
            spd.wcTmp.push_back(startPos);
            double randNum = dsfmt_genrand_open_close(&sfmt);
            nid_t incr = truncSkipLen(std::log(randNum) / log_p, maxJumpLen);
            log_trace("incr = %d, %f, %f", incr, randNum, log_p);
            startPos += incr;
          }

          if (spd.wcTmp.m_num == 0) {
            elemStatus[i] = FRONTIER_FINISHED;
          } else {
            for (auto j = 0; j < spd.wcTmp.m_num; j++) {
              _mm_prefetch((void*)(adjBase + spd.frontierOffset.first + spd.wcTmp[j]),
                           PREFETCH_HINT);
            }
          }

        } else if (nbrsNum == 1) {
          spd.wcTmp.push_back(0);
          _mm_prefetch((void*)(adjBase + spd.frontierOffset.first), PREFETCH_HINT);
        } else {
          elemStatus[i] = FRONTIER_FINISHED;
        }
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        if (elemStatus[i] == FRONTIER_WORKING) {
          auto& spd = elem.spd;
          const nid_t* addr = adjBase + spd.frontierOffset.first;
          for (auto j = 0; j < spd.wcTmp.m_num; j++) {
            const nid_t currNbr = spd.wcTmp[j];
            log_trace("step 3: j = %d, nbr=%d", j, currNbr);
            const nid_t target = addr[currNbr];
            spd.wcTmp[j] = target;
            _mm_prefetch((void*)(elem.visitedFlag.getData(target)), PREFETCH_HINT);
          }
        }
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        if (elemStatus[i] == FRONTIER_WORKING) {
          auto& spd = elem.spd;
          nid_t cursor = 0;
          for (auto j = 0; j < spd.wcTmp.m_num; j++) {
            const nid_t currNbrId = spd.wcTmp[j];
            if (!elem.visitedFlag.getBool(currNbrId)) {
              spd.rrPtr->push_back(currNbrId);
              elem.visitedFlag.setTrue(currNbrId);
            }
          }
          elemStatus[i] = FRONTIER_FINISHED;
        }
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];
        if (elemIsIdle[i]) continue;

        assert(elemStatus[i] == FRONTIER_FINISHED);

        auto& spd = elem.spd;

        spd.frontier++;
        if (spd.frontier < spd.rrPtr->m_num) {
          elemStatus[i] = FRONTIER_WORKING;
          continue;
        } else {
          log_trace("RR set complete: i=%d, numCompleted: %ld, numWorking: %ld", i, numCompleted,
                    numWorking);
          numCompleted++;
        }

        {  // restore the environment
          for (auto i = 0; i < spd.rrPtr->m_num; i++) {
            log_trace("set False: %ld", elem.spd.rrPtr->m_data[i]);
            elem.visitedFlag.setFalse(elem.spd.rrPtr->m_data[i]);
          }

          optRRs.push_back(elem.spd.rrPtr);
          elem.spd.rrPtr = new iVector<nid_t>();
        }

        if (numWorking < numRRsets) {
          elemStatus[i] = FRONTIER_WORKING;
          elem.spd.frontier = 0;
          int t = targets[numWorking++];
          elem.spd.rrPtr->push_back(t);
          elem.visitedFlag.setTrue(t);
        } else {
          elemIsIdle[i] = true;
        }
      }
    }

    int64_t totalNum = optRRs.size();
    log_info("total rrset: %ld", totalNum);
  }

  void batchLT_Execute(nid_t* targets, rrid_t numRRsets) {
    const nid_t maxJumpLen = graphSize + 1;

    numWorking = 0;
    numCompleted = 0;
    for (bid_t i = 0; i < BATCH_SIZE; i++) {
      auto& elem = elems[i];
      elem.spd.rrPtr = new iVector<nid_t>();
      nid_t target = targets[numWorking++];
      log_trace("working id: %ld, target: %ld", numWorking, target);
      elem.spd.rrPtr->push_back(target);
      elem.spd.ltTarget = target;
      elem.visitedFlag.setTrue(target);

      elemIsIdle[i] = false;
      elemStatus[i] = LT_WORKING;
    }

    while (numCompleted < numRRsets) {
      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        log_trace("step 1: i=%d, status: %d, idle: %d", i, elemStatus[i], elemIsIdle[i]);
        auto& elem = elems[i];
        if (elemIsIdle[i] || elemStatus[i] != LT_WORKING) continue;

        auto& spd = elem.spd;
        log_trace("elem batch id: %d, frontier: %u, frontierNodeId: %u", i, spd.rrPtr->m_num,
                  spd.ltTarget);
        _mm_prefetch((void*)(offsetBase + spd.ltTarget), PREFETCH_HINT);
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];
        if (elemIsIdle[i] || elemStatus[i] != LT_WORKING) continue;

        auto& spd = elem.spd;
        spd.frontierOffset = offsetBase[spd.ltTarget];
        nid_t nbrsNum = spd.frontierOffset.second - spd.frontierOffset.first;
        log_trace("i: %d, first: %ld, second: %ld, size: %d", i, spd.frontierOffset.first,
                  spd.frontierOffset.second, nbrsNum);

        if (nbrsNum != 0) {
          const int64_t nbr = dsfmt_genrand_uint32(&sfmt) % nbrsNum;
          spd.ltPos = (adjBase + spd.frontierOffset.first + nbr);
          _mm_prefetch((void*)(spd.ltPos), PREFETCH_HINT);
        } else {
          elemStatus[i] = RRSET_COMPLETED;
        }
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];
        if (elemIsIdle[i] || elemStatus[i] != LT_WORKING) continue;

        auto& spd = elem.spd;
        spd.ltTarget = *(spd.ltPos);

        _mm_prefetch((void*)(elem.visitedFlag.getData(spd.ltTarget)), PREFETCH_HINT);
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];
        if (elemIsIdle[i] || elemStatus[i] != LT_WORKING) continue;

        auto& spd = elem.spd;
        if (!elem.visitedFlag.getBool(spd.ltTarget)) {
          spd.rrPtr->push_back(spd.ltTarget);
          elem.visitedFlag.setTrue(spd.ltTarget);
        } else {
          elemStatus[i] = RRSET_COMPLETED;
        }
      }

      for (bid_t i = 0; i < BATCH_SIZE; i++) {
        auto& elem = elems[i];

        // restore environment
        if (elemIsIdle[i] || elemStatus[i] == LT_WORKING) continue;

        auto& spd = elem.spd;
        for (auto i = 0; i < spd.rrPtr->m_num; i++) {
          log_trace("set False: %ld", elem.spd.rrPtr->m_data[i]);
          elem.visitedFlag.setFalse(elem.spd.rrPtr->m_data[i]);
        }

        optRRs.push_back(elem.spd.rrPtr);
        elem.spd.rrPtr = new iVector<nid_t>();

        numCompleted++;
        if (numWorking < numRRsets) {
          elemStatus[i] = LT_WORKING;
          elem.spd.frontier = 0;
          int t = targets[numWorking++];
          log_trace("working id: %ld, target: %ld", numWorking, t);
          elem.spd.rrPtr->push_back(t);
          elem.spd.ltTarget = t;
          elem.visitedFlag.setTrue(t);
        } else {
          elemIsIdle[i] = true;
        }
      }
    }

    int64_t totalNum = optRRs.size();
    log_info("total rrset: %ld", totalNum);
  }

  void geneIC_RRsets(rrid_t prevSize, rrid_t numRRsets, InvList_t* invList) {
    rrid_t leftNumRRsets = numRRsets - prevSize;
    assert(leftNumRRsets >= 0);
    log_info("graph Size: %u", graphSize);
    nid_t* targets = geneTargets((rrid_t)leftNumRRsets, graphSize);

    if (spdCas == IC && spdModel == WC) {
      log_info("cascade: IC, prob dist: WC");
      batchWCOpt_Execute(targets, leftNumRRsets);
    } else if (spdCas == LT && spdModel == WC) {
      log_info("cascade: LT, prob dist: WC");
      batchLT_Execute(targets, leftNumRRsets);
    } else {
      batchICOpt_Execute(targets, leftNumRRsets);
    }

    delete[] targets;
#ifdef BUILD_INV_LIST
    createInvListWithPrefecth(invList, prevSize);
#endif
  }

  nid_t* geneTargets(size_t numRRsets, nid_t graphSize) {
    nid_t* targets = new nid_t[numRRsets];
    for (rrid_t i = 0; i < numRRsets; i++) {
      targets[i] = dsfmt_genrand_uint32(&sfmt) % graphSize;
      ;
    }

    return targets;
  }

  void createInvListWithPrefecth(InvList_t* invList, rrid_t startOffset = 0) {
    Timer timeInvertedIndex;

    size_t cnt = 0;
    // log_info("start: %ld, RR number: %ld", startOffset, optRRs.size());
    size_t totalLen = optRRs.size();
    for (auto i = startOffset; i < totalLen; i++) {
      const auto& meta = *optRRs[i];
      const nid_t* rrset = meta.m_data;
      const size_t hyperidx = i;

      if (i < totalLen - 1) {
        const auto nextMetaPtr = optRRs[i + 1];
        _mm_prefetch((void*)nextMetaPtr, PREFETCH_HINT);
      }

      const size_t len = meta.m_num;
      cnt += len;
      size_t curr = 0;
      {
        while (BATCH_SIZE >= 16 && curr + 16 <= len) {
          const uint32_t id0 = rrset[curr];
          const uint32_t id1 = rrset[curr + 1];
          const uint32_t id2 = rrset[curr + 2];
          const uint32_t id3 = rrset[curr + 3];
          const uint32_t id4 = rrset[curr + 4];
          const uint32_t id5 = rrset[curr + 5];
          const uint32_t id6 = rrset[curr + 6];
          const uint32_t id7 = rrset[curr + 7];
          const uint32_t id8 = rrset[curr + 8];
          const uint32_t id9 = rrset[curr + 9];
          const uint32_t id10 = rrset[curr + 10];
          const uint32_t id11 = rrset[curr + 11];
          const uint32_t id12 = rrset[curr + 12];
          const uint32_t id13 = rrset[curr + 13];
          const uint32_t id14 = rrset[curr + 14];
          const uint32_t id15 = rrset[curr + 15];

          _mm_prefetch((void*)(&invList[id0]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id1]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id2]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id3]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id4]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id5]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id6]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id7]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id8]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id9]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id10]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id11]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id12]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id13]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id14]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id15]), PREFETCH_HINT);

          void* addr0 = invList[id0].get_next_position();
          _mm_prefetch(addr0, PREFETCH_HINT);

          void* addr1 = invList[id1].get_next_position();
          _mm_prefetch(addr1, PREFETCH_HINT);

          void* addr2 = invList[id2].get_next_position();
          _mm_prefetch(addr2, PREFETCH_HINT);

          void* addr3 = invList[id3].get_next_position();
          _mm_prefetch(addr3, PREFETCH_HINT);

          void* addr4 = invList[id4].get_next_position();
          _mm_prefetch(addr4, PREFETCH_HINT);

          void* addr5 = invList[id5].get_next_position();
          _mm_prefetch(addr5, PREFETCH_HINT);

          void* addr6 = invList[id6].get_next_position();
          _mm_prefetch(addr6, PREFETCH_HINT);

          void* addr7 = invList[id7].get_next_position();
          _mm_prefetch(addr7, PREFETCH_HINT);

          void* addr8 = invList[id8].get_next_position();
          _mm_prefetch(addr8, PREFETCH_HINT);

          void* addr9 = invList[id9].get_next_position();
          _mm_prefetch(addr9, PREFETCH_HINT);

          void* addr10 = invList[id10].get_next_position();
          _mm_prefetch(addr10, PREFETCH_HINT);

          void* addr11 = invList[id11].get_next_position();
          _mm_prefetch(addr11, PREFETCH_HINT);

          void* addr12 = invList[id12].get_next_position();
          _mm_prefetch(addr12, PREFETCH_HINT);

          void* addr13 = invList[id13].get_next_position();
          _mm_prefetch(addr13, PREFETCH_HINT);

          void* addr14 = invList[id14].get_next_position();
          _mm_prefetch(addr14, PREFETCH_HINT);

          void* addr15 = invList[id15].get_next_position();
          _mm_prefetch(addr15, PREFETCH_HINT);

          *(size_t*)addr0 = hyperidx;
          *(size_t*)addr1 = hyperidx;
          *(size_t*)addr2 = hyperidx;
          *(size_t*)addr3 = hyperidx;
          *(size_t*)addr4 = hyperidx;
          *(size_t*)addr5 = hyperidx;
          *(size_t*)addr6 = hyperidx;
          *(size_t*)addr7 = hyperidx;
          *(size_t*)addr8 = hyperidx;
          *(size_t*)addr9 = hyperidx;
          *(size_t*)addr10 = hyperidx;
          *(size_t*)addr11 = hyperidx;
          *(size_t*)addr12 = hyperidx;
          *(size_t*)addr13 = hyperidx;
          *(size_t*)addr14 = hyperidx;
          *(size_t*)addr15 = hyperidx;

          invList[id0].increase_num();
          invList[id1].increase_num();
          invList[id2].increase_num();
          invList[id3].increase_num();
          invList[id4].increase_num();
          invList[id5].increase_num();
          invList[id6].increase_num();
          invList[id7].increase_num();
          invList[id8].increase_num();
          invList[id9].increase_num();
          invList[id10].increase_num();
          invList[id11].increase_num();
          invList[id12].increase_num();
          invList[id13].increase_num();
          invList[id14].increase_num();
          invList[id15].increase_num();

          curr += 16;
        }
        while (BATCH_SIZE >= 8 && curr + 8 <= len) {
          const uint32_t id0 = rrset[curr];
          const uint32_t id1 = rrset[curr + 1];
          const uint32_t id2 = rrset[curr + 2];
          const uint32_t id3 = rrset[curr + 3];
          const uint32_t id4 = rrset[curr + 4];
          const uint32_t id5 = rrset[curr + 5];
          const uint32_t id6 = rrset[curr + 6];
          const uint32_t id7 = rrset[curr + 7];

          _mm_prefetch((void*)(&invList[id0]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id1]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id2]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id3]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id4]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id5]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id6]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id7]), PREFETCH_HINT);

          void* addr0 = invList[id0].get_next_position();
          _mm_prefetch(addr0, PREFETCH_HINT);

          void* addr1 = invList[id1].get_next_position();
          _mm_prefetch(addr1, PREFETCH_HINT);

          void* addr2 = invList[id2].get_next_position();
          _mm_prefetch(addr2, PREFETCH_HINT);

          void* addr3 = invList[id3].get_next_position();
          _mm_prefetch(addr3, PREFETCH_HINT);

          void* addr4 = invList[id4].get_next_position();
          _mm_prefetch(addr4, PREFETCH_HINT);

          void* addr5 = invList[id5].get_next_position();
          _mm_prefetch(addr5, PREFETCH_HINT);

          void* addr6 = invList[id6].get_next_position();
          _mm_prefetch(addr6, PREFETCH_HINT);

          void* addr7 = invList[id7].get_next_position();
          _mm_prefetch(addr7, PREFETCH_HINT);

          *(size_t*)addr0 = hyperidx;
          *(size_t*)addr1 = hyperidx;
          *(size_t*)addr2 = hyperidx;
          *(size_t*)addr3 = hyperidx;
          *(size_t*)addr4 = hyperidx;
          *(size_t*)addr5 = hyperidx;
          *(size_t*)addr6 = hyperidx;
          *(size_t*)addr7 = hyperidx;

          invList[id0].increase_num();
          invList[id1].increase_num();
          invList[id2].increase_num();
          invList[id3].increase_num();
          invList[id4].increase_num();
          invList[id5].increase_num();
          invList[id6].increase_num();
          invList[id7].increase_num();

          curr += 8;
        }

        while (BATCH_SIZE >= 4 && curr + 4 <= len) {
          const uint32_t id0 = rrset[curr];
          const uint32_t id1 = rrset[curr + 1];
          const uint32_t id2 = rrset[curr + 2];
          const uint32_t id3 = rrset[curr + 3];

          _mm_prefetch((void*)(&invList[id0]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id1]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id2]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id3]), PREFETCH_HINT);

          void* addr0 = invList[id0].get_next_position();
          _mm_prefetch(addr0, PREFETCH_HINT);

          void* addr1 = invList[id1].get_next_position();
          _mm_prefetch(addr1, PREFETCH_HINT);

          void* addr2 = invList[id2].get_next_position();
          _mm_prefetch(addr2, PREFETCH_HINT);

          void* addr3 = invList[id3].get_next_position();
          _mm_prefetch(addr3, PREFETCH_HINT);

          *(size_t*)addr0 = hyperidx;
          *(size_t*)addr1 = hyperidx;
          *(size_t*)addr2 = hyperidx;
          *(size_t*)addr3 = hyperidx;

          invList[id0].increase_num();
          invList[id1].increase_num();
          invList[id2].increase_num();
          invList[id3].increase_num();
          curr += 4;
        }

        while (BATCH_SIZE >= 2 && curr + 2 <= len) {
          const uint32_t id0 = rrset[curr];
          const uint32_t id1 = rrset[curr + 1];

          _mm_prefetch((void*)(&invList[id0]), PREFETCH_HINT);
          _mm_prefetch((void*)(&invList[id1]), PREFETCH_HINT);

          void* addr0 = invList[id0].get_next_position();
          _mm_prefetch(addr0, PREFETCH_HINT);

          void* addr1 = invList[id1].get_next_position();
          _mm_prefetch(addr1, PREFETCH_HINT);

          *(size_t*)addr0 = hyperidx;
          *(size_t*)addr1 = hyperidx;

          invList[id0].increase_num();
          invList[id1].increase_num();

          curr += 2;
        }

        for (; curr < len; curr++) {
          const uint32_t id = rrset[curr];
          invList[id].push_back(hyperidx);
        }
      }
    }
  }
};
#endif
