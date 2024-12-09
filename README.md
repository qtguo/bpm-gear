# bpm

# Description
The technical report and implementation code for the paper "Efficient Algorithms for Budgeted Profit Maximization with Theoretical Guarantees". The codes are implemented based on the codes of SUBSIM: https://github.com/qtguo/subsim.git

# Guideline for running code
1. compile code: 
```
   make
```

2. Format the dataset
```
    bash download.sh
   ./bpm -func=format -gname=pokec
```

2. Run Gear on Dataset pokec
```shell
# Running Gear Algorithm
# Budget 100
./bpm -func=evalGreedy -gname=pokec -rrnum=10240000 -budget=100 -greedyType=gear -model=ic -sampleOpt=1

```
