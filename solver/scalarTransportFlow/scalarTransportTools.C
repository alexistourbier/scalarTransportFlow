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
    Set of functions to:
        - interpolate cell centered values onto cell faces consistently with isoAdvector
        - correct diffusion coefficient in a two-field framework
\*---------------------------------------------------------------------------*/

#include "scalarTransportTools.H"

void Foam::bulkInterpolate
(
    const volScalarField& C,
    const boolList& interfaceCells,
    const surfaceScalarField& Cf_2ndOrder,
    surfaceScalarField& Cf
)
{
    const fvMesh& mesh = C.mesh();

    forAll(C, cellI)
    {
        if (interfaceCells[cellI] == 0)
        {
            labelList faces = mesh.cells()[cellI];
            forAll(faces, faceI)
            {
                if (mesh.isInternalFace(faces[faceI]))
                {
                    Cf[faces[faceI]] = Cf_2ndOrder[faces[faceI]];
                }
            }
        }
    }
}

void Foam::bulkDiffusivity
(
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const dimensionedScalar& Diffusivity,
    surfaceScalarField& D
)
{
    const scalar alphaTol = 1e-014;
    const fvMesh& mesh = alpha.mesh();
    const cellList& cellFaces = mesh.cells(); 
    forAll(alpha, cellI)
    {
        
        if (interfaceCells[cellI] == 0)
        {
            const cell& faces = cellFaces[cellI];
            
    
            forAll(faces, faceI)
            {
                if (mesh.isInternalFace(faces[faceI]))
                {
                    if (alpha[cellI] > (1 - alphaTol))
                    {
                        D[faces[faceI]] = Diffusivity.value();
                    }
                    else
                    {
                        D[faces[faceI]] = 0;
                    }
                }
            }
        }
    }
}

void Foam::interfaceInterpolate
(
    const volScalarField& C,
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const surfaceScalarField& alphaf,
    const surfaceScalarField& alphaPhi,
    surfaceScalarField& Cf
)
{
    const fvMesh& mesh = Cf.mesh();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(interfaceCells, cellI)
    {
        labelList interfaces = mesh.cells()[cellI];
        if (interfaceCells[cellI] == 1)
        {
            forAll(interfaces, faceI)
            {
                if (mesh.isInternalFace(interfaces[faceI]))
                {
                    const scalar COwn = C[owner[interfaces[faceI]]];
                    const scalar CNei = C[neighbour[interfaces[faceI]]];
                    const scalar alphaOwn = alpha[owner[interfaces[faceI]]];
                    const scalar alphaNei = alpha[neighbour[interfaces[faceI]]];
                    const scalar alphafaceI = alphaf[interfaces[faceI]];
                    const scalar alphaPhiFaceI = alphaPhi[interfaces[faceI]];
                    Cf.primitiveFieldRef()[interfaces[faceI]] = interpolateOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                }
                else
                {
                    const label patchI = boundaryMesh.whichPatch(interfaces[faceI]);        // index of the patch to which faceI belongs to
                    const label bfaceI = boundaryMesh[patchI].whichFace(interfaces[faceI]); // - mesh.nInternalFaces(); // local faceI index
                    if (Cf.boundaryField()[patchI].coupled())
                    {
                        const scalar COwn = C.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar CNei = C.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphaOwn = alpha.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar alphaNei = alpha.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphafaceI = alphaf.boundaryField()[patchI][bfaceI];
                        const scalar alphaPhiFaceI = alphaPhi.boundaryField()[patchI][bfaceI];
                        Cf.boundaryFieldRef()[patchI][bfaceI] = interpolateOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                    }
                    // else // to be implemented
                    // {
                    //     Cf.boundaryFieldRef()[patch][bfaceI] = alphaf.boundaryField()[patch][bfaceI]*C.primitiveField()[faceCells[bFaceI]]/alpha.primitiveField()[cellI];
                    // }
                }
            }
        }
    }
}

void Foam::phaseAverageInterpolate
(
    const volScalarField& C,
    const volScalarField& alpha,
    const boolList& interfaceCells,
    const surfaceScalarField& alphaf,
    const surfaceScalarField& alphaPhi,
    surfaceScalarField& Cf
)
{
    const fvMesh& mesh = Cf.mesh();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    forAll(interfaceCells, cellI)
    {
        labelList interfaces = mesh.cells()[cellI];
        if (interfaceCells[cellI] == 1)
        {
            forAll(interfaces, faceI)
            {
                if (mesh.isInternalFace(interfaces[faceI]))
                {
                    const scalar COwn = C[owner[interfaces[faceI]]];
                    const scalar CNei = C[neighbour[interfaces[faceI]]];
                    const scalar alphaOwn = alpha[owner[interfaces[faceI]]];
                    const scalar alphaNei = alpha[neighbour[interfaces[faceI]]];
                    const scalar alphafaceI = alphaf[interfaces[faceI]];
                    const scalar alphaPhiFaceI = alphaPhi[interfaces[faceI]];
                    Cf.primitiveFieldRef()[interfaces[faceI]] = interpolateFaceConcentrationOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                }
                else
                {
                    const label patchI = boundaryMesh.whichPatch(interfaces[faceI]);        // index of the patch to which faceI belongs to
                    const label bfaceI = boundaryMesh[patchI].whichFace(interfaces[faceI]); // - mesh.nInternalFaces(); // local faceI index
                    if (Cf.boundaryField()[patchI].coupled())
                    {
                        const scalar COwn = C.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar CNei = C.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphaOwn = alpha.boundaryField()[patchI].patchInternalField()()[bfaceI];
                        const scalar alphaNei = alpha.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                        const scalar alphafaceI = alphaf.boundaryField()[patchI][bfaceI];
                        const scalar alphaPhiFaceI = alphaPhi.boundaryField()[patchI][bfaceI];
                        Cf.boundaryFieldRef()[patchI][bfaceI] = interpolateFaceConcentrationOnFaceI(alphafaceI, COwn, CNei, alphaOwn, alphaNei, alphaPhiFaceI);
                    }
                    // else // to be implemented
                    // {
                    //     Cf.boundaryFieldRef()[patch][bfaceI] = alphaf.boundaryField()[patch][bfaceI]*C.primitiveField()[faceCells[bFaceI]]/alpha.primitiveField()[cellI];
                    // }
                }
            }
        }
    }
}


Foam::scalar Foam::interpolateOnFaceI
(
    const scalar& alphaf,
    const scalar& COwn,
    const scalar& CNei,
    const scalar& alphaOwn,
    const scalar& alphaNei,
    const scalar& alphaPhiFaceI
)
{
    const scalar alphaTol = 1e-012;
    scalar Cf = 0;
    if (alphaf > alphaTol)
    {
        if (alphaPhiFaceI > 0 && alphaOwn > 0)
        {
            Cf = alphaf * COwn / alphaOwn;
        }
        if (alphaPhiFaceI < 0 && alphaNei > 0)
        {
            Cf = alphaf * CNei / alphaNei;
        }
    }
    else
    {
        Cf = 0;
    }
    return Cf;
}

Foam::scalar Foam::interpolateFaceConcentrationOnFaceI
(
    const scalar& alphaf,
    const scalar& COwn,
    const scalar& CNei,
    const scalar& alphaOwn,
    const scalar& alphaNei,
    const scalar& alphaPhiFaceI
)
{
    const scalar alphaTol = 1e-012;
    scalar Cf = 0;
    if (alphaf > alphaTol)
    {
        if (alphaPhiFaceI > 0 && alphaOwn > 0)
        {
            Cf = alphaf * COwn;
        }
        if (alphaPhiFaceI < 0 && alphaNei > 0)
        {
            Cf = alphaf * CNei;
        }
    }
    else
    {
        Cf = 0;
    }
    return Cf;
}


void Foam::limitFaceValues
(
    surfaceScalarField& Cf
)
{
    forAll(Cf, faceI)
    {
        if (Cf[faceI] < 0)
        {
            Cf[faceI] = 0;
        }    
    }
}

void Foam::correctDiffusivity
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
)
{
    const fvMesh& mesh = alpha.mesh();
    const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    // cut off the bulk diffusivity on every face of an interface cell
    forAll(mesh.faces(), faceI)
    {
        if (mesh.isInternalFace(faceI))
        {
            const label& own = owner[faceI];
            const label& nei = neighbour[faceI];
            const scalar alphaOwn = alpha[own];
            const scalar C1Own = C1[own];
            const scalar C2Own = C2[own];
            const scalar alphaNei = alpha[nei];
            const scalar C1Nei = C1[nei];
            const scalar C2Nei = C2[nei];
            const scalar alphaOwnNew = alphaNew[own];
            const scalar alphaNeiNew = alphaNew[nei];
            const scalar alphafFaceI = alphaf[faceI];
            Deff[faceI] = diffusivityOfFaceI(alphaOwn, C1Own, C2Own, alphaNei, C1Nei, C2Nei, alphaOwnNew, alphaNeiNew, alphafFaceI, H, D);
        }
        else
        {
            const label patchI = boundaryMesh.whichPatch(faceI);        // index of the patch to which faceI belongs to
            const label bfaceI = boundaryMesh[patchI].whichFace(faceI); // - mesh.nInternalFaces(); // local faceI index
            if (Deff.boundaryField()[patchI].coupled())
            {
                const scalar alphaOwn = alpha.boundaryField()[patchI].patchInternalField()()[bfaceI];
                const scalar C1Own = C1.boundaryField()[patchI].patchInternalField()()[bfaceI];
                const scalar C2Own = C2.boundaryField()[patchI].patchInternalField()()[bfaceI];
                const scalar alphaOwnNew = alphaNew.boundaryField()[patchI].patchInternalField()()[bfaceI];
                const scalar alphaNei = alpha.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                const scalar C1Nei = C1.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                const scalar C2Nei = C2.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                const scalar alphaNeiNew = alphaNew.boundaryField()[patchI].patchNeighbourField()()[bfaceI];
                const scalar alphafFaceI = alphaf.boundaryField()[patchI][bfaceI];
                Deff.boundaryFieldRef()[patchI][bfaceI] = diffusivityOfFaceI(alphaOwn, C1Own, C2Own, alphaNei, C1Nei, C2Nei, alphaOwnNew, alphaNeiNew, alphafFaceI, H, D.value());
            }
        }
    }
}

void Foam::interfaceDiffusivity
(
    const dimensionedScalar& D,
    volScalarField& Di,
    const boolList& interfaceCells
)
{
    // set the diffusivity at the interface, only for the mass source term
    forAll(Di, cellI)
    {
        if (interfaceCells[cellI] == 1)
        {
            Di[cellI] = D.value();
        }
        else
        {
            Di[cellI] = 0;
        }
    }
}

Foam::scalar Foam::diffusivityOfFaceI
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
)
{
    scalar Deff = 0;
    if (H.value() == 0)
    {
        if (alphaOwn < 1 || alphaNei < 1)
        {
            Deff = 0;
        }
        if (alphaOwnNew < 1 || alphaNeiNew < 1)
        {
            Deff = 0;
        }
    }
    else
    {
        if ((alphaOwn) > 0 && (alphaNei) > 0)
        {
            Deff = alphaf * D.value();
        }
    }
    return Deff;
}

void Foam::syncProcPatches
(
    surfaceScalarField& Cf,
    const boolList& interfaceCells
)
{
    // references to the mesh and patches
    const fvMesh& mesh = Cf.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const fvPatchList& patches = mesh.boundary();
    // list containing the patches indexes
    DynamicList<label> procPatchLabels(pbm.size());
    // list for the communication of the surfaceScalarField
    List<DynamicList<label>> surfaceCellFacesOnProcPatches(0);
    List<DynamicList<scalar>> pFaceValues(0);
    // list used to verify if the patch needs information from the other processor
    List<DynamicList<bool>> requiresNeighbData(0);
    surfaceCellFacesOnProcPatches.resize(pbm.size());
    pFaceValues.resize(pbm.size());
    requiresNeighbData.resize(pbm.size());
    // Append all processor patch labels to the list
    forAll(patches, patchI)
    {
        const fvPatch& thisPatch = patches[patchI];
        const polyPatch& thisPolyPatch = thisPatch.patch();
        if (isA<processorPolyPatch>(thisPolyPatch) && thisPolyPatch.size() > 0)
        {
            procPatchLabels.append(patchI);
            forAll(thisPolyPatch, faceI)
            {
                surfaceCellFacesOnProcPatches[patchI].append(thisPolyPatch.start() + faceI);
                const label cellI = thisPatch.faceCells()[faceI];
                if (interfaceCells[cellI] == 0)
                {
                    requiresNeighbData[patchI].append(1);
                }
                else
                {
                    requiresNeighbData[patchI].append(0);
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels, i)
        {
            const label patchI = procPatchLabels[i];
            const processorPolyPatch &procPatch =
                refCast<const processorPolyPatch>(pbm[patchI]);
            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            forAll(pbm[patchI], faceI)
            {
                pFaceValues[patchI].append(Cf.boundaryField()[patchI][faceI]);
            }
            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches[patchI];
            const List<scalar>& faceValues =
                pFaceValues[patchI];
            toNbr << surfCellFacesOnProcPatch << faceValues;
        }

        pBufs.finishedSends();

        // Receive and combine
        forAll(procPatchLabels, patchLabelI)
        {
            const label patchI = procPatchLabels[patchLabelI];
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pbm[patchI]);
            List<label> faceIDs;
            List<scalar> nbrCfs;
            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            fromNeighb >> faceIDs >> nbrCfs;
            scalarField& localFaceValues = Cf.boundaryFieldRef()[patchI];
            forAll(faceIDs, i)
            {
                if (requiresNeighbData[patchI][i] == 1)
                {
                    localFaceValues[i] = nbrCfs[i];
                }
            }
        }
    }
    surfaceCellFacesOnProcPatches.clear();
    pFaceValues.clear();
    requiresNeighbData.clear();
}

Foam::volScalarField Foam::degreeOfFilling
(
    const volScalarField& Cliq,
    const volScalarField& Cgas,
    const dimensionedScalar& H
)
{
    volScalarField degreeOfFilling = Cliq * 0;
    forAll (Cliq, cellI)
    {
        scalar C1 = Cliq[cellI];
        scalar C2 = Cgas[cellI];
        if (C2 > 0)
        {
            degreeOfFilling[cellI] = C1 / (H.value() * C2);
        }
        else
        {
            degreeOfFilling[cellI] = 0;
        }
    }
    return degreeOfFilling;
}

Foam::label Foam::getCellLabel
(
    const label& localIndex,
    const volScalarField& C
)
{
    const fvMesh& mesh = C.mesh();
    if (localIndex < mesh.nCells())
    {
        return localIndex; // cellI
    }
    else // cell next to it
    {
        const label faceI = localIndex + mesh.nInternalFaces() - mesh.nCells();

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        // Boundary face. Find out which face of which patch
        const label patchI = pbm.whichPatch(faceI);
        if (patchI < 0 || patchI >= pbm.size())
        {
            FatalErrorInFunction << "Cannot find patch for face " << faceI
                                 << abort(FatalError);
        }
        const polyPatch& pp = pbm[patchI];
        const label patchFaceI = pp.whichFace(faceI);
        const label cellI = pp.faceCells()[patchFaceI];
        return cellI;
    }
}

void Foam::getNeighbours
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
)
{
    const vector n = faceNormal / mag(faceNormal);
    const labelListList& stencil = exchangeFields.getStencil();
    for (label i = 1; i < stencil[cellI].size(); ++i)
    {
        const label gblIdx = stencil[cellI][i];
        const vector neiCC = exchangeFields.getValue(centre, mapCC, gblIdx);
        const vector dist = neiCC - faceCentre;
        scalar cosAngle = (dist / mag(dist)) & n;
        if (cosAngle > 0.25) // roughly 75 deg
        {
            scalar neiC = exchangeFields.getValue(C, mapC, gblIdx);
            scalar neiAlpha = exchangeFields.getValue(alpha, mapAlpha, gblIdx);
            CField.append(neiC);
            alphaField.append(neiAlpha);
            distField.append(dist);
            indicies.append(gblIdx);
        }
    }
}

