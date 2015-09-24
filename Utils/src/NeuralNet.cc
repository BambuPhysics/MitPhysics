#include "MitPhysics/Utils/interface/NeuralNet.h"


NeuralNet::NeuralNet(unsigned int in, unsigned int out):
  integrityChecked(false)
{
  nIn = in;
  nOut = out;
}

NeuralNet::~NeuralNet() { 
  for (float *layer : layers)
    delete[] layer;
  for (float *b : bs) {
    delete[] b;
  }
  for (unsigned int j=0; j!=nIn; ++j) {
    delete[]Ws[0][j];
  }
  for (unsigned int i=1; i!=Ws.size(); ++i) {
    for (unsigned int j=0; j!=hiddenLayerSizes[i-1]; ++j) {
      delete[]Ws[i][j];
    }
  }
  for (float **W : Ws) {
    delete[] W;
  }
}

void NeuralNet::AddLayer(unsigned int in, unsigned int out, float **W, float *b, bool isFinal) {
  unsigned int previousLayerSize = (hiddenLayerSizes.size()==0) ? nIn : hiddenLayerSizes.back();
  assert(previousLayerSize==in);
  if (isFinal)
    assert(nOut==out);
  // layer is valid
  if (!isFinal)
    hiddenLayerSizes.push_back(out);
  Ws.push_back(W);
  bs.push_back(b);
  integrityChecked = false;
}

void NeuralNet::AddBranchAddress(float *input, float mean, float stdev, const char *name/*=""*/) {
  inputs.push_back(input);
  mus.push_back(mean);
  sigmas.push_back(stdev);
  inputNames.push_back(name);
  integrityChecked = false;
}

void NeuralNet::AllocateMemory() {
  // pre-allocate memory for linear operations
  // faster than using TMatrix, which allocates memory for each operation
  float *tmpLayer;
  for (unsigned int iL=0; iL!=hiddenLayerSizes.size(); ++iL) {
    fprintf(stderr,"allocating layer of size %i\n",(int)hiddenLayerSizes[iL]);
    tmpLayer = new float[hiddenLayerSizes[iL]];
    layers.push_back(tmpLayer);
  }
  fprintf(stderr,"allocating layer of size %i\n",(int)nOut);
  tmpLayer = new float[nOut];
  layers.push_back(tmpLayer);
  integrityChecked = false;
}

bool NeuralNet::CheckIntegrity() const {
  // TODO: check input variables are ordered correctly
  if (inputs.size()!=nIn) {
    fprintf(stderr,"Network has inconsistent structure: %i!=%i\n",(int)inputs.size(),(int)nIn);
    return false;
  }
  unsigned int nW = Ws.size();
  if (hiddenLayerSizes.size()+1!=nW || nW!=bs.size() || nW!=layers.size()) {
    fprintf(stderr,"Network has inconsistent structure: %i %i %i %i\n",(int)hiddenLayerSizes.size(),(int)nW,(int)bs.size(),(int)layers.size());
    return false;
  }
  for (unsigned int i=0; i!=Ws.size(); ++i) {
    if (Ws[i]==NULL || bs[i]==NULL || layers[i]==NULL) {
      fprintf(stderr,"Network layer %i is not correctly provided: %p %p %p\n",(int)i,Ws[i],bs[i],layers[i]);
      return false;
    }
  }
  integrityChecked = true;
  return true;
}

float *NeuralNet::Evaluate() const {
	if (!integrityChecked)
  	assert(CheckIntegrity());
  for (unsigned int iOut=0; iOut!=hiddenLayerSizes[0]; ++iOut) {
    float outVal = bs[0][iOut];
    for (unsigned int iIn=0; iIn!=nIn; ++iIn) {
      float tmpIn = ( *(inputs[iIn]) - mus[iIn] ) / sigmas[iIn];
      outVal += tmpIn * Ws[0][iIn][iOut];
    }
    layers[0][iOut] = TMath::TanH(outVal);
  }
  for (unsigned int iLayer=1; iLayer!=hiddenLayerSizes.size(); ++iLayer) {
   for (unsigned int iOut=0; iOut!=hiddenLayerSizes[iLayer]; ++iOut) {
      float outVal = bs[iLayer][iOut];
      for (unsigned int iIn=0; iIn!=hiddenLayerSizes[iLayer-1]; ++iIn) {
        outVal += layers[iLayer-1][iIn] * Ws[iLayer][iIn][iOut];
        
      }
      layers[iLayer][iOut] = TMath::TanH(outVal);
    } 
  }
  float Z = 0;
  float *outLayer = layers.back();
  unsigned int nextToLast = hiddenLayerSizes.size();
  float *prevLayer = layers[nextToLast-1];
  for (unsigned int iOut=0; iOut!=nOut; ++iOut) {
    float outVal = bs[nextToLast][iOut];
    for (unsigned int iIn=0; iIn!=hiddenLayerSizes.back(); ++iIn) {
      outVal += prevLayer[iIn] * Ws[nextToLast][iIn][iOut];
    }
    outLayer[iOut] = TMath::Exp(outVal);
    Z += outLayer[iOut];
  }
  for (unsigned int iOut=0; iOut!=nOut; ++iOut) {
  	outLayer[iOut] /= Z;
  }
  return outLayer;
}

