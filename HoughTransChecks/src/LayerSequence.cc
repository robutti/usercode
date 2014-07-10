#include <iostream>

#include "../interface/LayerSequence.h"

LayerSequence::LayerSequence(int nConsLayers) :
  nConsLayers_(nConsLayers)
{
  layerPattern_.reserve(1<<nLayers);

  // Initialize allowed next layers
  for (int iLyr = 0; iLyr < nLayers; iLyr++) {
    nextLayers_[iLyr] = new vector<int>;
  }
  nextLayers_[0]->push_back(1); nextLayers_[0]->push_back(3);
  nextLayers_[1]->push_back(2); nextLayers_[1]->push_back(3);
  nextLayers_[2]->push_back(3); nextLayers_[2]->push_back(5);
  nextLayers_[3]->push_back(4); nextLayers_[3]->push_back(5);
  nextLayers_[4]->push_back(5); nextLayers_[4]->push_back(18); nextLayers_[4]->push_back(19); nextLayers_[4]->push_back(20);
  nextLayers_[5]->push_back(6); nextLayers_[5]->push_back(9);
  nextLayers_[6]->push_back(9); nextLayers_[6]->push_back(12);
  nextLayers_[9]->push_back(10); nextLayers_[9]->push_back(12); nextLayers_[9]->push_back(18); nextLayers_[9]->push_back(19); nextLayers_[9]->push_back(20);
  nextLayers_[10]->push_back(11); nextLayers_[10]->push_back(12); nextLayers_[10]->push_back(18); nextLayers_[10]->push_back(19); nextLayers_[10]->push_back(20); nextLayers_[10]->push_back(21);
  nextLayers_[11]->push_back(12); nextLayers_[11]->push_back(18); nextLayers_[11]->push_back(19); nextLayers_[11]->push_back(20); nextLayers_[11]->push_back(21); nextLayers_[11]->push_back(22); nextLayers_[11]->push_back(23); nextLayers_[11]->push_back(24); nextLayers_[11]->push_back(25); nextLayers_[11]->push_back(26);
  nextLayers_[12]->push_back(13); nextLayers_[12]->push_back(18);
  nextLayers_[13]->push_back(18);
  nextLayers_[18]->push_back(19); nextLayers_[18]->push_back(21); nextLayers_[18]->push_back(22); nextLayers_[18]->push_back(23); nextLayers_[18]->push_back(24); nextLayers_[18]->push_back(25); nextLayers_[18]->push_back(26);
  nextLayers_[19]->push_back(20); nextLayers_[19]->push_back(21); nextLayers_[19]->push_back(22); nextLayers_[19]->push_back(23); nextLayers_[19]->push_back(24); nextLayers_[19]->push_back(25); nextLayers_[19]->push_back(26);
  nextLayers_[20]->push_back(21); nextLayers_[20]->push_back(22); nextLayers_[20]->push_back(23); nextLayers_[20]->push_back(24); nextLayers_[20]->push_back(25); nextLayers_[20]->push_back(26);
  nextLayers_[21]->push_back(22); nextLayers_[21]->push_back(23); nextLayers_[21]->push_back(24); nextLayers_[21]->push_back(25); nextLayers_[21]->push_back(26);
  nextLayers_[22]->push_back(23); nextLayers_[22]->push_back(25); nextLayers_[22]->push_back(26);
  nextLayers_[24]->push_back(25);
  nextLayers_[25]->push_back(26);

  // Initialize layer ID map
  layerId_[17] = 0; layerId_[18] = 1; layerId_[19] = 2;
  layerId_[33] = 3; layerId_[34] = 4;
  layerId_[49] = 5; layerId_[50] = 6; layerId_[51] = 7; layerId_[52] = 8;
  layerId_[65] = 9; layerId_[66] = 10; layerId_[67] = 11;
  layerId_[81] = 12; layerId_[82] = 13; layerId_[83] = 14; layerId_[84] = 15; layerId_[85] = 16; layerId_[86] = 17;
  layerId_[97] = 18; layerId_[98] = 19; layerId_[99] = 20; layerId_[100] = 21; layerId_[101] = 22; layerId_[102] = 23; layerId_[103] = 24; layerId_[104] = 25; layerId_[105] = 26;

  // Loop over layer patterns
  for (int iPat = 0; iPat < (1<<nLayers); iPat++) {
    layerPattern_[iPat] = false;
    for (int iLyr = 0; iLyr < nLayers - nConsLayers + 1; iLyr++)
      if ((iPat & (1 << iLyr)) && FindNext_(iPat, iLyr)) {
        layerPattern_[iPat] = true;
        break;
      }
  }
}

//

LayerSequence::~LayerSequence()
{
  for (int iLyr = 0; iLyr < nLayers; iLyr++)
    delete nextLayers_[iLyr];
}

//

bool LayerSequence::FindNext_(int pattern, int lastLayer, unsigned int nLayersFound)
{
  if (pattern == 0)
    return false;

  // Look for allowed next layer
  for (vector<int>::iterator itLyr = nextLayers_[lastLayer]->begin(); itLyr != nextLayers_[lastLayer]->end(); itLyr++) {
    if (pattern & (1 << *itLyr)) {
      if (nLayersFound == nConsLayers_ - 1)
        return true;
      else if (FindNext_(pattern, *itLyr, nLayersFound + 1))
        return true;
    }
  }
  return false;         
}

//

bool LayerSequence::testPattern(vector<unsigned int> layerVector)
{
  int pattern = 0;
  for (vector<unsigned int>::iterator itLyr = layerVector.begin(); itLyr != layerVector.end(); itLyr++)
    pattern |= (1 << layerId_[*itLyr]);
  return testPattern(pattern);
}

//

