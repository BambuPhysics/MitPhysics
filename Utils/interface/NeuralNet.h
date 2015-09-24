#ifndef NEURALNET
#define NEURALNET 1

#include <vector>
#include "TMath.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

class NeuralNet
{
public:
  NeuralNet(unsigned int in, unsigned int out);
  ~NeuralNet();
 
  void AddLayer(unsigned int in, unsigned int out, float **W, float *b, bool isFinal = false);
  void AddBranchAddress(float *input, float mean, float stdev, const char *name="");
  void AllocateMemory();

  bool CheckIntegrity () const;
  float* Evaluate () const;

protected:
  unsigned int nIn;
  unsigned int nOut;
  std::vector<unsigned int> hiddenLayerSizes;
  std::vector<float*> inputs;
  std::vector<const char*> inputNames;
  std::vector<float**> Ws;
  std::vector<float*> bs;
  std::vector<float*> layers;
  std::vector<float> mus;
  std::vector<float> sigmas;

  mutable  bool integrityChecked;

};

#endif
