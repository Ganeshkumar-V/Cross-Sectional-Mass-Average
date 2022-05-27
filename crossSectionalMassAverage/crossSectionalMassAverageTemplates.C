/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2022 Ganeshkumar V, IIT Gandhinagar, India.
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::crossSectionalMassAverage::calcAverage()
{
  typedef GeometricField<vector, fvPatchField, volMesh> volVectorFieldx;
  typedef GeometricField<scalar, fvPatchField, volMesh> volScalarFieldx;

    if (foundObject<volVectorFieldx>(fieldName_, false))
    {
	const volVectorFieldx U = lookupObject<volVectorFieldx>(fieldName_);
	const volScalarField rho = lookupObject<volScalarField>("thermo:rho.air");
	volScalarField U_mag = mag(U);

	// Finding the number of CVs in the radial direction
	label ncells = 0;
	for (int i = 1; i < mesh_.V().size(); i++)
	{
	   if (
	       (mesh_.C()[i].z() > mesh_.C()[ncells].z()) &&
	       ((mesh_.C()[i].z() - mesh_.C()[ncells].z()) > 1e-8)
	      ) break;
	   ncells = ncells + 1;
	}
	ncells = ncells + 1;

	//  Find Average Velocity
	//volScalarField* U_avg = new volScalarField(0);
	volScalarField U_avg = U_mag;
	double sum_rhoUV = 0;
	double sum_rhoV = 0;
	double sum_avg = 0;
	int ndiv = mesh_.V().size()/ncells;
	for (int i = 0; i < ndiv; i++)
	{
	    int index = i*ncells;
	    // for all cells in radial direction
	    for (int n = 0; n < ncells; n++)
	    {
		int celli = index + n;
		sum_rhoUV = sum_rhoUV + U_mag.internalField()[celli]*mesh_.V()[celli]
			*rho.internalField()[celli];
		sum_rhoV = sum_rhoV + mesh_.V()[celli]*rho.internalField()[celli];
	    }
	    sum_avg = sum_rhoUV/sum_rhoV;
	    for (int n = 0; n < ncells; n++)
	    {
		U_avg[index+n] = sum_avg;
	    }
	    sum_rhoUV = 0;
	    sum_rhoV = 0;
	}
	return store
        (
            resultName_,
            Foam::mag(U_avg),
	    mesh_.changing() && mesh_.cache(resultName_)
        );
    }
    else if (foundObject<volScalarFieldx>(fieldName_, false))
    {
	const volScalarField U_mag = lookupObject<volScalarFieldx>(fieldName_);
        const volScalarField rho = lookupObject<volScalarField>("thermo:rho.air");
	// Finding the number of CVs in the radial direction
	label ncells = 0;
	for (int i = 1; i < mesh_.V().size(); i++)
	{
	   if (
	       (mesh_.C()[i].z() > mesh_.C()[ncells].z()) &&
	       ((mesh_.C()[i].z() - mesh_.C()[ncells].z()) > 1e-8)
	      ) break;
	   ncells = ncells + 1;
	}
	ncells = ncells + 1;

	//  Find Average Velocity
	//volScalarField* U_avg = new volScalarField(0);
	volScalarField U_avg = U_mag;
	double sum_rhoUV = 0;
	double sum_rhoV = 0;
	double sum_avg = 0;
	int ndiv = mesh_.V().size()/ncells;
	for (int i = 0; i < ndiv; i++)
	{
	    int index = i*ncells;
	    // for all cells in radial direction
	    for (int n = 0; n < ncells; n++)
	    {
		int celli = index + n;
		sum_rhoUV = sum_rhoUV + U_mag.internalField()[celli]*mesh_.V()[celli]
			*rho.internalField()[celli];
		sum_rhoV = sum_rhoV + mesh_.V()[celli]*rho.internalField()[celli];
	    }
	    sum_avg = sum_rhoUV/sum_rhoV;
	    for (int n = 0; n < ncells; n++)
	    {
		U_avg[index+n] = sum_avg;
	    }
	    sum_rhoUV = 0;
	    sum_rhoV = 0;
	}
	return store
        (
            resultName_,
            Foam::mag(U_avg),
	    mesh_.changing() && mesh_.cache(resultName_)
        );
    }

    return false;
}


// ************************************************************************* //
