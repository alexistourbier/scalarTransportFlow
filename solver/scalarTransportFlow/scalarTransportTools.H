
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "explicitInterfaceDiffFlux.H"
#include "zoneDistribute.H"
#include "processorPolyPatch.H"

namespace Foam
{

// Interpolation in the bulk of the phase i.e. if the cell is not an interface cell
void bulkInterpolate
(
    const volScalarField& C,
    const boolList& interfaceCells,
    const surfaceScalarField& Cf_2ndOrder,
    surfaceScalarField& Cf
);

// Computes the diffusion coefficient if the cell is not an interface cell
void bulkDiffusivity
(
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const dimensionedScalar& Diffusivity,
    surfaceScalarField& D
);


// Interpolation of a scalar in an interface cell, consistently with isoAdvector
void interfaceInterpolate
(
    const volScalarField& C,
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const surfaceScalarField& alphaf,
    const surfaceScalarField& alphaPhi,
    surfaceScalarField& Cf
);
void phaseAverageInterpolate
(
    const volScalarField& C,
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const surfaceScalarField& alphaf,
    const surfaceScalarField& alphaPhi,
    surfaceScalarField& Cf
);
// function that does the final computation used in interfaceInterpolate
scalar interpolateOnFaceI
(
    const scalar& alphaf,
    const scalar& COwn,
    const scalar& CNei,
    const scalar& alphaOwn,
    const scalar& alphaNei,
    const scalar& alphaPhiFaceI
);
scalar interpolateFaceConcentrationOnFaceI
(
    const scalar& alphaf,
    const scalar& COwn,
    const scalar& CNei,
    const scalar& alphaOwn,
    const scalar& alphaNei,
    const scalar& alphaPhiFaceI
);
// function used to cut the diffusivity on faces of an interface cell
void correctDiffusivity
(
    const volScalarField& alpha,
    const volScalarField& alphaNew,
    const surfaceScalarField& alphaf,
    const volScalarField& C1,
    const volScalarField& C2,
    const dimensionedScalar& D,
    surfaceScalarField& Deff,
    const boolList& interfaceCells,
    const dimensionedScalar& H
);

// compute the diffusion coefficient for the mass transfer through the interface
void interfaceDiffusivity
(
    const dimensionedScalar& D,
    volScalarField& Di,
    const boolList& interfaceCells
);

// limiters to ensure the face values stay between 0 and 1
void limitFaceValues
(
    surfaceScalarField& Cf
);



scalar diffusivityOfFaceI
(
    const scalar& alphaOwn,
    const scalar& C1Own,
    const scalar& C2Own,
    const scalar& alphaNei,
    const scalar& C1Nei,
    const scalar& C2Nei,
    const scalar& alphaOwnNew,
    const scalar& alphaNeiNew,
    const scalar& alphaf,
    const dimensionedScalar& H,
    const dimensionedScalar& D
);

// synchronisation of processor patches when the value has been updated by interfaceInterpolate on one processor
void syncProcPatches
(
    surfaceScalarField& Cf,
    const boolList& interfaceCells
);

// computation of the interfacial concentration from the liquid concentration field
// volScalarField interfaceLiquidConcentration
// (
//     const volScalarField& Cliq,
//     const volScalarField& Cgas,
//     const scalar& H,
//     const boolList& interfaceCells
// );

// computation of the interfacial concentration from the gas concentration field
volScalarField interfaceConcentration
(
    const volScalarField& Cliq,
    const volScalarField& Cgas,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dimensionedScalar& H,
    const boolList& interfaceCells
);

volScalarField degreeOfFilling
(
    const volScalarField& Cliq,
    const volScalarField& Cgas,
    const dimensionedScalar& H
);

label getCellLabel
(
    const label& localIndex,
    const volScalarField& C
);

void getNeighbours
(
    const label cellI,
    const vector faceCentre,
    const vector faceNormal,
    const volVectorField& centre,
    const volScalarField& C,
    const volScalarField& alpha,
    const Map<vector>& mapCC,
    const Map<scalar>& mapC,
    const Map<scalar>& mapAlpha,
    zoneDistribute& exchangeFields,
    DynamicField<scalar>& CField,
    DynamicField<scalar>& alphaField,
    DynamicField<vector>& distField,
    DynamicField<label>& indicies
);

volScalarField explicitInterfaceDiffFlux
(
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const volScalarField& CLVi,
    const volScalarField& CGVi,
    const volScalarField& CLVLi,
    const volScalarField& CGVGi,
    const boolList& interfaceCells_0,
    const boolList& interfaceCells,
    const volVectorField& interfaceCenters_0,
    const volVectorField& interfaceCenters,
    const volVectorField& interfaceNormals_0,
    const volVectorField& interfaceNormals,
    const volScalarField& D,
    const dimensionedScalar& H
);
}