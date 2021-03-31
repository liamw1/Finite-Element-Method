#include "Mesh2D.h"

Mesh2D::Mesh2D(const int n, const int nodes, const int edges)
  : size(n), numNodes(nodes), numEdges(edges)
{
}

Mesh2D::Mesh2D(Mesh2D&& other) noexcept
  : size(other.size),
    numNodes(other.numNodes),
    numEdges(other.numEdges),
    connectivityMatrix(std::move(other.connectivityMatrix)),
    edgeArray(std::move(other.edgeArray)),
    edgeTypeMatrix(std::move(other.edgeTypeMatrix)),
    edgeMatrix(std::move(other.edgeMatrix))
{
  meshNodes = other.meshNodes;
  other.meshNodes = nullptr;
}

Mesh2D::~Mesh2D()
{
}
