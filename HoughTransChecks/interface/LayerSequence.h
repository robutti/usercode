#ifndef HoughTransChecks_LayerSequence_H
#define HoughTransChecks_LayerSequence_H
#include <vector>
#include <map>

using namespace std;

class LayerSequence {
 public:
  enum {nLayers = 27};

 private:
  vector<int>* nextLayers_[nLayers];
  map<unsigned int, int> layerId_;
  vector<bool> layerPattern_;  // 2^(n. of layers) elements
  unsigned int nConsLayers_;
  bool FindNext_(int pattern, int lastLayer, unsigned int nLayersFound = 1);

 public:
  LayerSequence(int nConsLayers = 3);
  ~LayerSequence();
  bool testPattern(int pattern) const {return layerPattern_[pattern];}
  bool testPattern(vector<unsigned int> layerVector);
  template <typename T> bool testPattern(map<unsigned int, T> layerMap);
};

template <typename T>
bool LayerSequence::testPattern(map<unsigned int, T> layerMap)
{
  int pattern = 0;
  for (typename map<unsigned int, T>::iterator itLyr = layerMap.begin(); itLyr != layerMap.end(); itLyr++)
    pattern |= (1 << layerId_[(*itLyr).first]);
  return testPattern(pattern);
}

#endif
