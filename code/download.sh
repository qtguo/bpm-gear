set -x
mkdir graphInfo
cd graphInfo
wget https://snap.stanford.edu/data/soc-pokec-relationships.txt.gz
gzip -d soc-pokec-relationships.txt.gz 
mv soc-pokec-relationships.txt pokec
