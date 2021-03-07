#include "Precompilied.h"
#include "Gauss-LegendreNodes.h"

static const std::vector<real> nullVec = {};

static const std::vector<real> gauss1DNodes2 = { -0.5773502691896257645092,
                                                  0.5773502691896257645092 };

static const std::vector<real> gauss1DWeights2 = { 1,
                                                   1 };

static const std::vector<real> gauss1DNodes3 = { -0.7745966692414833770359,
                                                  0,
                                                  0.7745966692414833770359 };

static const std::vector<real> gauss1DWeights3 = { 0.5555555555555555555556,
                                                   0.8888888888888888888889,
                                                   0.5555555555555555555556 };

static const std::vector<real> gauss1DNodes4 = { -0.861136311594052575224,
                                                 -0.3399810435848562648027,
                                                  0.3399810435848562648027,
                                                  0.861136311594052575224 };

static const std::vector<real> gauss1DWeights4 = { 0.3478548451374538573731,
                                                   0.6521451548625461426269,
                                                   0.6521451548625461426269,
                                                   0.3478548451374538573731 };

static const std::vector<real> gauss1DNodes5 = { -0.9061798459386639927976,
                                                 -0.5384693101056830910363,
                                                  0,
                                                  0.5384693101056830910363,
                                                  0.9061798459386639927976 };

static const std::vector<real> gauss1DWeights5 = { 0.2369268850561890875143,
                                                   0.4786286704993664680413,
                                                   0.5688888888888888888889,
                                                   0.4786286704993664680413,
                                                   0.2369268850561890875143 };

const std::vector<real>& gauss1DNodesRef(const int numNodes)
{
  ASSERT(numNodes > 1, "Number of nodes must be greater than 1");

  switch (numNodes)
  {
  case 2:
    return gauss1DNodes2;
  case 3:
    return gauss1DNodes3;
  case 4:
    return gauss1DNodes4;
  case 5:
    return gauss1DNodes5;
  default:
    LOG("Number of nodes not supported", LogLevel::Error);
    return nullVec;
  }
}

const std::vector<real>& gauss1DWeightsRef(const int numNodes)
{
  ASSERT(numNodes > 1, "Number of nodes must be greater than 1");

  switch (numNodes)
  {
  case 2:
    return gauss1DWeights2;
  case 3:
    return gauss1DWeights3;
  case 4:
    return gauss1DWeights4;
  case 5:
    return gauss1DWeights5;
  default:
    LOG("Number of nodes not supported", LogLevel::Error);
    return nullVec;
  }
}

std::vector<real> gauss1DNodesLocal(const Mesh1D& mesh, const int elementIndex, const int numNodes)
{
  const real& a = mesh(elementIndex, EdgeType::Left).x;
  const real& b = mesh(elementIndex, EdgeType::Right).x;
  const std::vector<real>& t = gauss1DNodesRef(numNodes);

  std::vector<real> localNodes = std::vector<real>(numNodes);
  for (int i = 0; i < numNodes; ++i)
    localNodes[i] = (b - a) * t[i] / 2 + (a + b) / 2;
  return localNodes;
}

std::vector<real> gauss1DWeightsLocal(const Mesh1D& mesh, const int elementIndex, const int numNodes)
{
  const real& a = mesh(elementIndex, EdgeType::Left).x;
  const real& b = mesh(elementIndex, EdgeType::Right).x;
  const std::vector<real>& w = gauss1DWeightsRef(numNodes);

  std::vector<real> localWeights = std::vector<real>(numNodes);
  for (int i = 0; i < numNodes; ++i)
    localWeights[i] = (b - a) * w[i] / 2;
  return localWeights;
}
