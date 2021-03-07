#pragma once
#include "Precompilied.h"
#include "Mesh1D.h"

const std::vector<real>& gauss1DNodesRef(const int numNodes);

const std::vector<real>& gauss1DWeightsRef(const int numNodes);

std::vector<real> gauss1DNodesLocal(const Mesh1D& mesh, const int elementIndex, const int numNodes);

std::vector<real> gauss1DWeightsLocal(const Mesh1D& mesh, const int elementIndex, const int numNodes);