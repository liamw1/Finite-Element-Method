#pragma once
#include "Precompilied.h"
#include "Mesh1D.h"
#include "Mesh2D.h"
#include "Matrix.h"

const std::vector<real>& gauss1DNodesRef(const int numNodes);

const std::vector<real>& gauss1DWeightsRef(const int numNodes);

std::vector<real> gauss1DNodesLocal(const Mesh1D& mesh, const int elementIndex, const int numNodes);

std::vector<real> gauss1DWeightsLocal(const Mesh1D& mesh, const int elementIndex, const int numNodes);



const std::vector<std::array<real, 2>>& gauss2DNodesRef(const int numNodes);

const std::vector<real>& gauss2DWeightsRef(const int numNodes);

std::vector<std::array<real, 2>> gauss2DNodesLocal(const Mesh2D& mesh, const int elementIndex, const int numNodes);

std::vector<real> gauss2DWeightsLocal(const Mesh2D& mesh, const int elementIndex, const int numNodes);