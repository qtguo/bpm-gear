#pragma once

typedef std::vector<nid_t> Nodelist;
/// Edge structure, neighbor id and the edge weight
typedef std::pair<nid_t, float> Edge;
/// Edgelist structure from one source/target node
typedef std::vector<Edge> Edgelist;
/// Graph structure
typedef std::vector<Edgelist> Graph;
/// One forward reachable set
typedef std::vector<rrid_t> FRset;
/// A set of forward reachable sets
typedef std::vector<FRset> FRsets;
double decay_factor = 1.0;
/// One reverse reachable set
typedef std::vector<nid_t> RRset;
/// A set of reverse reachable sets
typedef std::vector<RRset> RRsets;

// typedef char int8_t;
// typedef unsigned char uint8_t;
// typedef long long int64_t;
// typedef unsigned long long uint64_t;

typedef std::pair<int64_t, int64_t> offset_t;

bool optflag;
enum CascadeModel { IC = 0, LT = 1 };
enum ProbDist { WEIGHTS, UNIFORM, WC, SKEWED, PROB_DIST_ERROR };
enum FuncType { FORMAT, PM, EVAL_GDY, EVAL_RRSET, IM, FUNC_ERROR };
enum GreedyType { GDY_GEAR = 0, GDY_SIMPLE = 1, GDY_DISTORTED = 2, GDY_SAMPLE = 3, GDY_ERROR };

typedef struct CostModel {
  bool randomCost = false;
  double basicCost = 1;
  double budget = 100;
  double beta = 0.005;
} CostModelCtx;

/// Node element with id and a property value
typedef struct NodeElement {
  int id;
  double value;
} NodeEleType;

/// Smaller operation for node element
struct smaller {
  bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const {
    return (Ele1.value < Ele2.value);
  }
};

/// Greater operation for node element
struct greater {
  bool operator()(const NodeEleType& Ele1, const NodeEleType& Ele2) const {
    return (Ele1.value > Ele2.value);
  }
};
