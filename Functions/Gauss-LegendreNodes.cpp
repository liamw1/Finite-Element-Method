#include "Precompilied.h"
#include "Gauss-LegendreNodes.h"

static const std::vector<real> nullVec = {};
static const std::vector<std::array<real, 2>> nullVec2D = {};

static const std::vector<real> gauss1DNodes2 = { -0.5773502691896257645092,
                                                  0.5773502691896257645092 };

static const std::vector<real> gauss1DWeights2 = { 1.0,
                                                   1.0 };

static const std::vector<real> gauss1DNodes3 = { -0.7745966692414833770359,
                                                  0.0,
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
                                                  0.0,
                                                  0.5384693101056830910363,
                                                  0.9061798459386639927976 };

static const std::vector<real> gauss1DWeights5 = { 0.2369268850561890875143,
                                                   0.4786286704993664680413,
                                                   0.5688888888888888888889,
                                                   0.4786286704993664680413,
                                                   0.2369268850561890875143 };

static const std::vector<real> gauss1DNodes6 = { -0.9324695142031520278123,
                                                 -0.661209386466264513661,
                                                 -0.2386191860831969086305,
                                                  0.238619186083196908631,
                                                  0.661209386466264513661,
                                                  0.9324695142031520278123 };

static const std::vector<real> gauss1DWeights6 = { 0.1713244923791703450403,
                                                   0.3607615730481386075698,
                                                   0.4679139345726910473899,
                                                   0.46791393457269104739,
                                                   0.3607615730481386075698,
                                                   0.1713244923791703450403 };

static const std::vector<real> gauss1DNodes7 = { -0.9491079123427585245262,
                                                 -0.7415311855993944398639,
                                                 -0.4058451513773971669066,
                                                  0.0,
                                                  0.4058451513773971669066,
                                                  0.7415311855993944398639,
                                                  0.9491079123427585245262 };

static const std::vector<real> gauss1DWeights7 = { 0.1294849661688696932706,
                                                   0.2797053914892766679015,
                                                   0.38183005050511894495,
                                                   0.417959183673469387755,
                                                   0.38183005050511894495,
                                                   0.279705391489276667901,
                                                   0.129484966168869693271 };



static const std::vector<std::array<real, 2>> gauss2DNodes7 = { { 1.012865073234563e-01, 1.012865073234563e-01 },
                                                                { 7.974269853530872e-01, 1.012865073234563e-01 },
                                                                { 1.012865073234563e-01, 7.974269853530872e-01},
                                                                { 4.701420641051151e-01, 4.701420641051151e-01 },
                                                                { 5.971587178976981e-02, 4.701420641051151e-01 },
                                                                { 4.701420641051151e-01, 5.971587178976981e-02 },
                                                                { 3.333333333333333e-01, 3.333333333333333e-01 } };

static const std::vector<real> gauss2DWeights7 = { 6.296959027241358e-02,
                                                   6.296959027241358e-02,
                                                   6.296959027241358e-02,
                                                   6.619707639425308e-02,
                                                   6.619707639425308e-02,
                                                   6.619707639425308e-02,
                                                   1.125000000000000e-01 };


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
  case 6:
    return gauss1DNodes6;
  case 7:
    return gauss1DNodes7;
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
  case 6:
    return gauss1DWeights6;
  case 7:
    return gauss1DWeights7;
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



const std::vector<std::array<real, 2>>& gauss2DNodesRef(const int numNodes)
{
  ASSERT(numNodes > 1, "Number of nodes must be greater than 1");

  switch (numNodes)
  {
  case 7:
    return gauss2DNodes7;
  default:
    LOG("Number of nodes not supported", LogLevel::Error);
    return nullVec2D;
  }
}

const std::vector<real>& gauss2DWeightsRef(const int numNodes)
{
  ASSERT(numNodes > 1, "Number of nodes must be greater than 1");

  switch (numNodes)
  {
  case 7:
    return gauss2DWeights7;
  default:
    LOG("Number of nodes not supported", LogLevel::Error);
    return nullVec;
  }
}

std::vector<std::array<real, 2>> gauss2DNodesLocal(const Mesh2D& mesh, const int elementIndex, const int numNodes)
{
  const int& K = elementIndex;
  const MeshNode2D& A1 = mesh(K, 0);
  const MeshNode2D& A2 = mesh(K, 1);
  const MeshNode2D& A3 = mesh(K, 2);
  const std::vector<std::array<real, 2>>& t = gauss2DNodesRef(numNodes);

  // Create transformation matrix from reference domain to K
  Matrix B = Matrix(2);
  B[0][0] = A2.x - A1.x;  B[0][1] = A3.x - A1.x;
  B[1][0] = A2.y - A1.y;  B[1][1] = A3.y - A1.y;

  std::vector<std::array<real, 2>> localNodes = std::vector<std::array<real, 2>>(numNodes);
  for (int i = 0; i < numNodes; ++i)
  {
    localNodes[i][0] = B[0][0] * t[i][0] + B[0][1] * t[i][1] + A1.x;
    localNodes[i][1] = B[1][0] * t[i][0] + B[1][1] * t[i][1] + A1.y;
  }
  return localNodes;
}

std::vector<real> gauss2DWeightsLocal(const Mesh2D& mesh, const int elementIndex, const int numNodes)
{
  const int& K = elementIndex;
  const MeshNode2D& A1 = mesh(K, 0);
  const MeshNode2D& A2 = mesh(K, 1);
  const MeshNode2D& A3 = mesh(K, 2);
  const std::vector<real>& w = gauss2DWeightsRef(numNodes);

  // Create transformation matrix from reference domain to K
  Matrix B = Matrix(2);
  B[0][0] = A2.x - A1.x;
  B[0][1] = A3.x - A1.x;
  B[1][0] = A2.y - A1.y;
  B[1][1] = A3.y - A1.y;

  // Calculate area of K
  Matrix A = Matrix(3);
  A[0][0] = A1.x;  A[0][1] = A1.y;  A[0][2] = 1.0;
  A[1][0] = A2.x;  A[1][1] = A2.y;  A[1][2] = 1.0;
  A[2][0] = A3.x;  A[2][1] = A3.y;  A[2][2] = 1.0;
  const real area = 0.5 * abs(A.determinant());
  
  // Debug
  ASSERT(area != 0.0, "Area matrix is singular");

  std::vector<real> localWeights = std::vector<real>(numNodes);
  for (int i = 0; i < numNodes; ++i)
    localWeights[i] = 2 * area * w[i];
  return localWeights;
}

std::vector<std::array<real, 2>> gaussEdgeNodesLocal(const Mesh2D& mesh, const int edgeIndex, const int numNodes)
{
  const MeshNode2D& A1 = mesh.meshNodes[mesh.edgeArray[edgeIndex][0]];
  const MeshNode2D& A2 = mesh.meshNodes[mesh.edgeArray[edgeIndex][1]];
  const std::vector<real>& t = gauss1DNodesRef(numNodes);

  std::vector<std::array<real,2>> localNodes = std::vector<std::array<real,2>>(numNodes);
  for (int i = 0; i < numNodes; ++i)
  {
    localNodes[i][0] = (A2.x - A1.x) * t[i] / 2 + (A1.x + A2.x) / 2;
    localNodes[i][1] = (A2.y - A1.y) * t[i] / 2 + (A1.y + A2.y) / 2;
  }
  return localNodes;
}

std::vector<real> gaussEdgeWeightsLocal(const Mesh2D& mesh, const int edgeIndex, const int numNodes)
{
  const MeshNode2D& A1 = mesh.meshNodes[mesh.edgeArray[edgeIndex][0]];
  const MeshNode2D& A2 = mesh.meshNodes[mesh.edgeArray[edgeIndex][1]];
  const std::vector<real>& w = gauss1DWeightsRef(numNodes);

  std::vector<real> localWeights = std::vector<real>(numNodes);
  for (int i = 0; i < numNodes; ++i)
  {
    const real distance = sqrtl((A2.x - A1.x) * (A2.x - A1.x) + (A2.y - A1.y) * (A2.y - A1.y));
    localWeights[i] = distance * w[i] / 2;
  }
  return localWeights;
}
