#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <sstream>
#include <vector>
#include <map>

#include "../interface/LayerSequence.h"

using namespace std;

int main(int argc, char* argv[])
{
  unsigned int layerId[27] = {17, 18, 19, 33, 34, 49, 50, 51, 52, 65, 66, 67, 81, 82, 83, 84, 85, 86, 97, 98, 99, 100, 101, 102, 103, 104, 105};

  LayerSequence lyrSeq;

  srand(time(NULL));
  for (int iPatt = 0; iPatt < (1 << 27); iPatt++) {
    if (!(rand()%100000)) {
      stringstream binPatt;
      vector<unsigned int> layerVector;
      map<unsigned int, bool> layerMap;
      for (int iBit = 26; iBit >= 0; iBit--) {
        unsigned int bit = (iPatt & (1 << iBit)) >> iBit;
        if (bit) {
          layerVector.push_back(layerId[iBit]);
          layerMap[layerId[iBit]] = true;
        }
        binPatt << bit;
      }
      cout << "Pattern " << binPatt.str() << ": ";
      cout << (lyrSeq.testPattern(iPatt) ? "OK" : "not OK") << "; "
           << (lyrSeq.testPattern(layerVector) ? "OK" : "not OK") << "; "
           << (lyrSeq.testPattern(layerMap) ? "OK" : "not OK");
      if ((lyrSeq.testPattern(iPatt) != lyrSeq.testPattern(layerVector)) ||
          (lyrSeq.testPattern(iPatt) != lyrSeq.testPattern(layerMap)))
        cout << " *** ERROR! ***";
      cout << endl;
    }
  }
  
  return 0;
}
