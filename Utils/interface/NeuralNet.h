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
 
  void AddLayer(unsigned int in, unsigned int out, double **W, double *b, bool isFinal = false);
  void AddBranchAddress(float *input, double mean, double stdev, const char *name="");
  void AddBranchAddress(float *input);
  void AddMuSigma(double mean, double stdev);
  void AllocateMemory();

  bool CheckIntegrity () const;
  double* Evaluate () const;

protected:
  unsigned int nIn;
  unsigned int nOut;
  std::vector<unsigned int> hiddenLayerSizes;
  std::vector<float*> inputs;
  std::vector<const char*> inputNames;
  std::vector<double**> Ws;
  std::vector<double*> bs;
  std::vector<double*> layers;
  std::vector<double> mus;
  std::vector<double> sigmas;

  mutable  bool integrityChecked;

};

#endif
