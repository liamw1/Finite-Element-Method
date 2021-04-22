#include "Precompilied.h"
#include "UnstructuredMesh2D.h"

static int initializeSize(const std::string meshFile)
{
  std::ifstream file(meshFile);

  int n = 0;
  std::string line;
  while (getline(file, line))
  {
    // Count each line between ">StartElements" and ">EndElements"
    if (line == ">StartElements")
    {
      getline(file, line);
      while (line != ">EndElements")
      {
        if (file.eof())
          LOG("List of elements has no ending marker", LogLevel::Error);

        ++n;
        getline(file, line);
      }
    }
  }
  
  file.close();
  return n;
}

static int initializeNumNodes(const std::string meshFile)
{
  std::ifstream file(meshFile);

  int n = 0;
  std::string line;
  while (getline(file, line))
  {
    // Count each line between ">StartNodes" and ">EndNodes"
    if (line == ">StartNodes")
    {
      getline(file, line);
      while (line != ">EndNodes")
      {
        if (file.eof())
          LOG("List of nodes has no ending marker", LogLevel::Error);

        ++n;
        getline(file, line);
      }
    }
  }

  file.close();
  return n;
}

static int initializeNumEdges(const std::string meshFile)
{
  return initializeSize(meshFile) + initializeNumNodes(meshFile) - 1;
}

static std::vector<double> split(const std::string& str, const char delim)
{
  std::vector<double> numbers;
  std::stringstream stream(str);
  std::string item;

  while (getline(stream, item, delim))
    numbers.push_back(std::stod(item));
  return numbers;
}

UnstructuredMesh2D::UnstructuredMesh2D(const std::string meshFile)
  : Mesh2D(initializeSize(meshFile), initializeNumNodes(meshFile), initializeNumEdges(meshFile))
{
  // Debug
  ASSERT(size > 0, "Invalid mesh size: A mesh must have a least one element!");

  // Initialize data structures
  meshNodes = new MeshNode2D[numNodes];
  connectivityMatrix = Array2D<int>(size, 3);
  edgeArray = Array2D<int>(numEdges, 2);
  edgeTypeMatrix = Array2D<int>(numNodes, numNodes);
  edgeMatrix = Array2D<int>(numNodes, numNodes);

  // Set node coordinates and connectivity matrix
  std::ifstream file(meshFile);
  std::string line;
  int nodeIndex = 0;
  int elementIndex = 0;
  while (getline(file, line))
  {
    // Read each line after ">StartNodes" and ending at ">EndNodes"
    if (line == ">StartNodes")
    {
      getline(file, line);
      while (line != ">EndNodes")
      {
        if (file.eof())
          LOG("List of nodes has no ending marker", LogLevel::Error);

        std::vector<double> coords = split(line, ' ');
#pragma warning(suppress: 6386)
        meshNodes[nodeIndex].x = coords[0];
        meshNodes[nodeIndex].y = coords[1];

        getline(file, line);
        ++nodeIndex;
      }
    }

    // Read each line after ">StartElements" and ending at ">EndElements"
    if (line == ">StartElements")
    {
      getline(file, line);
      while (line != ">EndElements")
      {
        if (file.eof())
          LOG("List of elements has no ending marker", LogLevel::Error);

        std::vector<double> element = split(line, ' ');
        for (int i = 0; i < 3; ++i)
          connectivityMatrix[elementIndex][i] = (int)element[i];

        getline(file, line);
        ++elementIndex;
      }
    }
  }
  file.close();

  // Form edge type matrix
  for (int i = 0; i < numNodes; ++i)
    for (int j = 0; j < numNodes; ++j)
      edgeTypeMatrix[i][j] = 0;
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
      {
        // Grabbing indices of j-th and k-th vertex in element
        const int& vj = connectivityMatrix[i][j];
        const int& vk = connectivityMatrix[i][k];

        ++edgeTypeMatrix[vj][vk];
      }

  // Form edge array
  int E = 0;
  for (int i = 0; i < numNodes; ++i)
    for (int j = i + 1; j < numNodes; ++j)
    {
      if (edgeTypeMatrix[i][j] > 0)
      {
        edgeArray[E][0] = i;
        edgeArray[E][1] = j;
        ++E;
      }
    }

  // Form edge matrix
  for (int i = 0; i < numNodes; ++i)
    for (int j = 0; j < numNodes; ++j)
      edgeMatrix[i][j] = 0;
  for (int i = 0; i < numEdges; ++i)
  {
    edgeMatrix[edgeArray[i][0]][edgeArray[i][1]] = i + 1;
    edgeMatrix[edgeArray[i][1]][edgeArray[i][0]] = i + 1;
  }
}

UnstructuredMesh2D::UnstructuredMesh2D(UnstructuredMesh2D&& other) noexcept
  : Mesh2D(other.size, other.numNodes, other.numEdges)
{
}

MeshNode2D UnstructuredMesh2D::operator()(const int elementIndex, const int nodeIndex) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < size, "Element index must be less than the number of elements");
  ASSERT(nodeIndex >= 0, "Node index must be non-negative");
  ASSERT(nodeIndex < 3, "Node index must be less than 3 on a triangular mesh");

  return meshNodes[connectivityMatrix[elementIndex][nodeIndex]];
}