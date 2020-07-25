/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "velocityABLWallModFunctionFvPatchField.H"
#include "volFields.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

velocityABLWallModFunctionFvPatchField::velocityABLWallModFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    printOn_(false)
{}


velocityABLWallModFunctionFvPatchField::velocityABLWallModFunctionFvPatchField
(
    const velocityABLWallModFunctionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    printOn_(ptf.printOn_)
{}


velocityABLWallModFunctionFvPatchField::velocityABLWallModFunctionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    printOn_(dict.lookupOrDefault<bool>("print", false))
{}


velocityABLWallModFunctionFvPatchField::velocityABLWallModFunctionFvPatchField
(
    const velocityABLWallModFunctionFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    printOn_(ptf.printOn_)
{}


velocityABLWallModFunctionFvPatchField::velocityABLWallModFunctionFvPatchField
(
    const velocityABLWallModFunctionFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    printOn_(ptf.printOn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityABLWallModFunctionFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

//  Set up access to the internal velocity field.
    const volVectorField& UCell = db().objectRegistry::lookupObject<volVectorField>("U");

//  Interpolate velocity on the cell faces.
    surfaceVectorField UFace = fvc::interpolate(UCell);

//  Set up access to the mesh
    const fvMesh& mesh = patch().boundaryMesh().mesh();

//  Get face areas (individual and global sum) (note that "gSum" is used as
//  opposed to "sum" because "gSum" is parallel-aware --- it gathers sums from
//  each processor to which this patch belongs and makes the global sum)
    const scalarField area = patch().magSf();
    scalar areaTotal = gSum(area);

//  Get perpendicular distance from cell center to boundary.  In other words,
//  the height of the z_1/2 grid level.
    // const scalarField z12 = 1.0/patch().deltaCoeffs();
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

//  Get face normal vectors
    vectorField normal = patch().nf();

//  Get resolved U vector at boundary face (UBoundaryFace refers to patch().lookupPatchField...)
    const fvPatchVectorField& UBoundaryFace = patch().lookupPatchField<volVectorField, vector>("U");

//  Get the components of the resolved velocity at z_1/2, U(z_1/2), that are
//  parallel to the boundary face, U_||(z_1/2), by using:
//
//        U_||(z_1/2) = U(z_1/2) - (U(z_1/2) dot |S_n|)*|S_n|
//
//  where S_n is the surface normal vector.
    vectorField UParallel12 = UBoundaryFace.patchInternalField();
    forAll(UBoundaryFace, faceI)
    {
	UParallel12[faceI] = UParallel12[faceI] - ((UParallel12[faceI] & normal[faceI]) * normal[faceI]);
    }

//  Now, the meat of the boundary condition: specify surface velocity!!
//  For the momentum equation with no molecular viscosity, it doesn't matter what
//  the velocity parallel to the surface is.  All that matters is that velocity
//  normal to the surface is zero.  The specified surface total stress creates
//  the surface drag to create a velocity profile.  However, we need to create
//  an appropriate surface normal velocity gradient at the z_1/2 level so that
//  the SGS model has a meaningful gradient to be used in its production term(s).
//  Therefore, we find the surface normal gradient at the z_1 level and "copy"
//  that to the z_1/2 level.  We do so by specifying a surface parallel velocity
//  that creates a surface normal gradient at z_1/2 equal to that at z_1.

//  First, initialize the velocity at the surface and level z_1 and the surface
//  normal gradient at level z_1/2.  ULocal is a reference to the object of this
//  class, the boundary velocity field itself.
    vectorField ULocal(patch().size(),vector::zero);
    vectorField snGradU1(patch().size(),vector::zero);
    vectorField U1(patch().size(),vector::zero);

    forAll(patch(), faceI)
    {
         //  Find the global cell indices for cells adjacent to this patch (the
         //  cells with centers at the z_1/2 level).
         label cellI = patch().faceCells()[faceI];

         //  Find the global face index for the face opposite the patch face (the
         //  z_1 level face).
         label oppFaceI = mesh.cells()[cellI].opposingFaceLabel(faceI+patch().patch().start(),mesh.faces());

         if(oppFaceI < 0)
         {
            FatalErrorInFunction
            << "Cell " << cellI  << " has " <<  mesh.cells()[cellI].nFaces()
            << " faces. Not a prism or hexagon"
            << exit(FatalError);
         }

         label oppPatchID = 0;
         if (oppFaceI >= UFace.size())
         {
             oppPatchID = mesh.boundaryMesh().whichPatch(oppFaceI);
             oppFaceI -= mesh.boundary()[oppPatchID].patch().start();
             U1[faceI] = UFace.boundaryField()[oppPatchID][oppFaceI];
             snGradU1[faceI] = (U1[faceI] - UParallel12[faceI])*deltaCoeffs[faceI];
         }
         else
         {
             U1[faceI] = UFace[oppFaceI];
             snGradU1[faceI] = (U1[faceI] - UParallel12[faceI])*deltaCoeffs[faceI];
         }

         ULocal[faceI] =  2.0*UParallel12[faceI] - U1[faceI];
         ULocal[faceI] -= ((ULocal[faceI] & normal[faceI]) * normal[faceI]);
    }

    vector ULocalAvg = gSum(ULocal * area) / areaTotal;
    vector U1Avg = gSum(U1 * area) / areaTotal;
    vector snGradU1Avg = gSum(snGradU1 * area) / areaTotal;

    if (printOn_)
    {
         Info << "<U_1> = " << U1Avg << tab << "<U_s> = " << ULocalAvg << tab << "<dU/dn> = " << snGradU1Avg << endl;
    }

    this->operator==(ULocal);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void velocityABLWallModFunctionFvPatchField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("print")     << printOn_   << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
   fvPatchVectorField,
   velocityABLWallModFunctionFvPatchField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
