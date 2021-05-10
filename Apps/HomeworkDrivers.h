#pragma once
#include "Precompilied.h"
#include "Utilities/Print.h"

// Meshes
#include "Meshing/1D/UniformMesh1D.h"
#include "Meshing/2D/UniformRectangularMesh2D.h"
#include "Meshing/2D/UnstructuredMesh2D.h"

// Equation systems
#include "EquationSystems/1D/Elliptic1DABCF.h"
#include "EquationSystems/2D/Elliptic2DABCF.h"
#include "EquationSystems/2D/StokesFluid.h"

#include "ErrorAnalysis/ErrorAnalysis.h"

void Interpolation_Driver();
void Hwk4_C1_Driver();
void Hwk4_C2_Driver();
void Hwk4_C3_Driver();
void Hwk6_B2_Driver();
void Hwk8_C1_Driver();
void Hwk9_C1_Driver();
void Hwk9_C2_Driver();
void Hwk9_C3_Driver();
void Hwk10_C1_Driver();
void Hwk10_C2_Driver();