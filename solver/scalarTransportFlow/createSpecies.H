 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Read transport properties dictionnary
    IOdictionary transportPropertiesDict
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info << "Reading solution species dictionary" << endl;
	const dictionary& solutionSpeciesDict = transportPropertiesDict.subDict("solutionSpecies");
	// const dictionary& liquidPhaseSpeciesDict = solutionSpeciesDict.subDict("liquid");
	// const dictionary& gasPhaseSpeciesDict = solutionSpeciesDict.subDict("gas");
	Info << "Reading solution species name\n" << endl;
	speciesTable solutionSpecies
	(
		solutionSpeciesDict.toc()
	);

	// speciesTable liquidSpecies
	// (
	// 	liquidPhaseSpeciesDict.toc()
	// );
	// speciesTable gasSpecies
	// (
	// 	gasPhaseSpeciesDict.toc()
	// );
	Info << "Create composition vector\n" << endl;
	// Liquid concentration averaged over the volume of liquid
	PtrList<volScalarField> CLVL(solutionSpecies.size());
	// Liquid concentration averaged over the volume of the cell (FVM)
	PtrList<volScalarField> CLV(solutionSpecies.size());
	// Gas concentration averaged over the volume of liquid
    PtrList<volScalarField> CGVG(solutionSpecies.size());
	// Gas concentration averaged over the volume of the cell (FVM)
    PtrList<volScalarField> CGV(solutionSpecies.size());

	// Create liquid phase volume and phase averaged concentration fields
	forAll(solutionSpecies, i)
    {
		IOobject headerLiqAveraged
		(
			"CLVL"+solutionSpecies[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		IOobject headerLiquidVolumeAveraged
		(
			"CLV"+solutionSpecies[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

    	// check if field exists and can be read
		if (headerLiqAveraged.typeHeaderOk<volScalarField>(true))
		{
			CLVL.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CLVL"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			volScalarField CLVLDefault
			(
				IOobject
				(
					"CLVLDefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);
			CLVL.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CLVL"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					CLVLDefault
				)
			);	
		}
    	// check if field exists and can be read
		if (headerLiquidVolumeAveraged.typeHeaderOk<volScalarField>(true))
		{
			CLV.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CLV"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			volScalarField CLVDefault
			(
				IOobject
				(
					"CLVDefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);
			CLV.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CLV"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					CLVDefault
				)
			);
		}
    	mesh.setFluxRequired(solutionSpecies[i]);
    }

	// Create liquid phase volume and phase averaged concentration fields

	forAll(solutionSpecies, i)
    {
		IOobject headerGasAveraged
		(
			"CGVG"+solutionSpecies[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		IOobject headerGasVolumeAveraged
		(
			"CGV"+solutionSpecies[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

    	// check if field exists and can be read
		if (headerGasAveraged.typeHeaderOk<volScalarField>(true))
		{
			CGVG.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CGVG"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			volScalarField CGVGDefault
			(
				IOobject
				(
					"CGVGDefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);
			CGVG.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CGVG"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					CGVGDefault
				)
			);
		}
    	// check if field exists and can be read
		if (headerGasVolumeAveraged.typeHeaderOk<volScalarField>(true))
		{
			CGV.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CGV"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			volScalarField CGVDefault
			(
				IOobject
				(
					"CGVDefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);
			CGV.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"CGV"+solutionSpecies[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					CGVDefault
				)
			);
		}
    	mesh.setFluxRequired(solutionSpecies[i]);
    }
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info << "Read species diffusion coefficients\n" << endl;
	PtrList<dimensionedScalar> DL(solutionSpecies.size());
	PtrList<dimensionedScalar> DG(solutionSpecies.size());
	forAll(solutionSpecies, i)
	{
		const dictionary& subDict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		DL.set
		(
			i,
			new dimensionedScalar("DL",dimensionSet(0,2,-1,0,0,0,0),subDict.lookupOrDefault("DL", 0.0))
		);
		DG.set
		(
			i,
			new dimensionedScalar("DG",dimensionSet(0,2,-1,0,0,0,0),subDict.lookupOrDefault("DG", 0.0))
		);
	}

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info << "Read species molecular weight\n" << endl;
	PtrList<dimensionedScalar> Mw(solutionSpecies.size());
	// PtrList<dimensionedScalar> MwG(gasSpecies.size());

	forAll(solutionSpecies, i)
	{
		const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		Mw.set
		(
			i,
			new dimensionedScalar("Mw",dimensionSet(1,0,0,0,-1,0,0),subdict.lookupOrDefault("Mw", 1.0))
		);
	}

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info << "Read species Henry's constant\n" << endl;
	PtrList<dimensionedScalar> H(solutionSpecies.size());

	forAll(solutionSpecies, i)
	{
		const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		H.set
		(
			i,
			new dimensionedScalar("H",dimensionSet(0,0,0,0,0,0,0),subdict.lookupOrDefault("H", 1.0))
		);
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info << "Read species liquid and gas saturations \n" << endl;
	PtrList<dimensionedScalar> CLsat(solutionSpecies.size());
	PtrList<dimensionedScalar> CGsat(solutionSpecies.size());

	forAll(solutionSpecies, i)
	{
		const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		CLsat.set
		(
			i,
			new dimensionedScalar("CLsat",dimensionSet(1,-3,0,0,0,0,0),subdict.lookupOrDefault("CLsat", 0.0))
		);
	}

	forAll(solutionSpecies, i)
	{
		const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		CGsat.set
		(
			i,
			new dimensionedScalar("CGsat",dimensionSet(1,-3,0,0,0,0,0),subdict.lookupOrDefault("CGsat", 0.0))
		);
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info << "Read stoechiometric coefficients\n" << endl;
	PtrList<dimensionedScalar> stoechioCoefficients(solutionSpecies.size());
	forAll(solutionSpecies, i)
	{
		const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		stoechioCoefficients.set
		(
			i,
			new dimensionedScalar("coeff",dimensionSet(0,0,0,0,0,0,0),subdict.lookupOrDefault("coeff", 0.0))
		);
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info << "Read reaction rate constants\n" << endl;
	PtrList<dimensionedScalar> kArrhenius(solutionSpecies.size());
	forAll(solutionSpecies, i)
	{
		const dictionary& subdict = solutionSpeciesDict.subDict(solutionSpecies[i]);
		kArrhenius.set
		(
			i,
			new dimensionedScalar("kArrhenius",dimensionSet(-1,3,-1,0,0,0,0),subdict.lookupOrDefault("kArrhenius", 0.0))
		);
	}

	//**************************************************************************
