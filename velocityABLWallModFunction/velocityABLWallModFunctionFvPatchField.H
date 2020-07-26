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

Class
    Foam::velocityABLWallModFunctionFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Wall boundary condition for velocity at the ground in atmospheric boundary
    layer incompressible flow in which there is wall roughness.  The velocity
    at the surface is specified such that the velocity gradient normal to the
    surface is equal to the normal velocity gradient at the cell face opposite
    to the boundary face, but there is no wall normal flow.  The purpose of this
    is to supply a proper vertical gradient for use in the turbulence model
    equation.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        print        | print output to screen  | no          | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
        lower
        {
           type            velocityABLWallModFunction;
           print           true;                       //-print output (true/false)
        }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    velocityABLWallModFunctionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef velocityABLWallModFunctionFvPatchField_H
#define velocityABLWallModFunctionFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class velocityABLWallModFunctionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class velocityABLWallModFunctionFvPatchField
:
    public fixedValueFvPatchVectorField
{

protected:

    // Protected data

        //- Print information
        bool printOn_;

        // ratio of friction velocity to von-Karman constant
        scalar  uStarByKappa_;

        // field to store opposing face indices
        labelField oppFaceIDs_;

public:

    //- Runtime type information
    TypeName("velocityABLWallModFunction");


    // Constructors

        //- Construct from patch and internal field
        velocityABLWallModFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        velocityABLWallModFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given velocityABLWallModFunctionFvPatchField
        //  onto a new patch
        velocityABLWallModFunctionFvPatchField
        (
            const velocityABLWallModFunctionFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        velocityABLWallModFunctionFvPatchField
        (
            const velocityABLWallModFunctionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new velocityABLWallModFunctionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        velocityABLWallModFunctionFvPatchField
        (
            const velocityABLWallModFunctionFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new velocityABLWallModFunctionFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            // find the opposite face ids of the wall patch faces
            void getOppositeFaceIDs();

            // return a const reference to the protected oppFaceIDs_
            const labelField& oppFaceIDs()
            {
                return oppFaceIDs_;
            }

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
