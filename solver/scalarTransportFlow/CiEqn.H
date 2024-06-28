{
    // Numerical constants
    dimensionedScalar SMALL
    (
        "SMALL",
        dimVol/dimTime,
        1e-50
    );
    // isoAdvector phase flux
    alphaPhi = advector->alphaPhi();
    // volScalarField alphaGas = 1 - alpha1.oldTime();
    // area of the face immersed in the liquid / total area of the face
    alphaf = alphaPhi/stabilise(phi.oldTime(), SMALL);
    // area of the face immersed in the gas / total area of the face
    alpha2f = 1 - alphaf;

    // Face interpolate density
    surfaceScalarField rhof = alphaf*rho1 + (1-alphaf)*rho2;
    volScalarField alphaOld = alpha1.oldTime();
    volScalarField alpha2Old = 1 - alphaOld;
    surfaceScalarField alpha2Phi = phi.oldTime() - alphaPhi;
    volVectorField interfaceVelocity = (interfaceCenters - interfaceCenters_0) / mesh.time().deltaT();

    typedef surfaceInterpolationScheme<scalar> schemeType;
    IStringStream schemeDataUD("upwind");
    tmp<schemeType> schemeTmpUD = schemeType::New(
        mesh,
        phi,
        schemeDataUD // Constructor requires a ref to Istream.
    );
    const schemeType& schemeUD = schemeTmpUD();

    IStringStream schemeDataVanLeer("vanLeer01");
    tmp<schemeType> schemeTmpVL = schemeType::New(
        mesh,
        phi,
        schemeDataVanLeer
    );
    const schemeType& schemeVL = schemeTmpVL();

    IStringStream schemeDataGasUD("upwind");
    tmp<schemeType> schemeTmpGas = schemeType::New(
        mesh,
        phi,
        schemeDataGasUD
    );
    const schemeType& schemeGas = schemeTmpGas();

    IStringStream schemeDataGasVanLeer("vanLeer");
    tmp<schemeType> schemeTmpGasVL = schemeType::New(
        mesh,
        phi,
        schemeDataGasVanLeer
    );
    const schemeType& schemeGasVL = schemeTmpGasVL();
    forAll(solutionSpecies, i)
    {
        // Henry's constant to turn on or off the mass transfer
        dimensionedScalar Hi = H[i];
        // Liquid diffusivity of species i
        dimensionedScalar DLi = DL[i];
        // Gas diffusivity of species i
        dimensionedScalar DGi = DG[i];

        // Liquid saturation concentration of species i
        dimensionedScalar CLsati = CLsat[i];
        // Gas saturation concentration of species i
        dimensionedScalar CGsati = CGsat[i];
        // Liquid concentration of species i
        volScalarField& CLVi = CLV[i];
        // Gas concentration of species i
        volScalarField& CGVi = CGV[i];
        // Liquid concentration averaged over the volume of liquid
        volScalarField& CLVLi = CLVL[i];
        // Gas concentration averaged over the volume of gas
        volScalarField& CGVGi = CGVG[i];

        CGVGi = CGVi / (alpha2.oldTime() + 1e-012);
        // Stoechiometric coefficient of species i in chemical reaction
        dimensionedScalar stoechioCoefficient = stoechioCoefficients[i];
        // Rate constant for reaction j
        dimensionedScalar kArrheniusi = kArrhenius[i];

        // Face interpolate diffisusivities
        surfaceScalarField DLeffi = fvc::interpolate(alpha1)*DLi;
        surfaceScalarField DGeffi = fvc::interpolate(alpha2)*DGi;
        // Interface diffusivities for the mass transfer
        volScalarField DLInt = 0*alpha1*DLi.value();
        volScalarField DGInt = 0*alpha2*DGi.value();
        // Weights for the interpolation (bulk)
        surfaceScalarField interpolationWeightsVAi = alphaf*0;
        surfaceScalarField interpolationWeightsVAVLi = alphaf*0;
        surfaceScalarField interpolationWeightsLAi = alphaf*0;
        surfaceScalarField interpolationWeightsLAVLi = alphaf*0;

        surfaceScalarField interpolationWeightsGasi = alpha2f*0;
        surfaceScalarField interpolationWeightsVLGasi = alpha2f*0;

    /************************ Face interpolation  *************************/

        /************************* Liquid phase  *************************/
        // Interpolation weights obtained for the volume averaged concentration with an UD scheme
        interpolationWeightsVAi = schemeUD.weights(CLVi);
        // Interpolation weights obtained for the liquid phase averaged concentration with an UD scheme
        // interpolationWeightsLAi = schemeUD.weights(CLVLi);
        // Interpolation weights obtained for the volume averaged concentration with a VL scheme
        interpolationWeightsVAVLi = schemeVL.weights(CLVi);
        // Interpolation weights obtained for the liquid phase averaged concentration with a VL scheme
        // interpolationWeightsLAVLi = schemeVL.weights(CLVLi);
        // Volume averaged liquid concentration interpolated with an UD scheme onto cell faces
        surfaceScalarField CfLVUDi = schemeUD.interpolate(CLVi, interpolationWeightsVAi);
        // Volume averaged liquid concentration interpolated with a VL scheme onto cell faces
        surfaceScalarField CfLVVLi = schemeVL.interpolate(CLVi, interpolationWeightsVAVLi);
        // Liquid phase concentration interpolated with an UD scheme onto cell faces
        // surfaceScalarField CfLVLUDi = schemeUD.interpolate(CLVLi, interpolationWeightsLAi);
        // Liquid phase concentration interpolated with a VL scheme onto cell faces
        // surfaceScalarField CfLVLVLi = schemeVL.interpolate(CLVLi, interpolationWeightsLAVLi);

        /************************* Gas phase  *************************/
        // Interpolation weights obtained for an UD scheme
        interpolationWeightsGasi = schemeUD.weights(CGVi);
        // Interpolation weights obtained for a VL scheme
        interpolationWeightsVLGasi = schemeVL.weights(CGVi);
        // Volume averaged gas concentration interpolated with UD scheme onto cell faces
        surfaceScalarField CfGVi = schemeUD.interpolate(CGVi, interpolationWeightsGasi);
        // Volume averaged gas concentration interpolated with VL scheme onto cell faces
        surfaceScalarField CfGVVLi = schemeVL.interpolate(CGVi, interpolationWeightsVLGasi);

        surfaceScalarField CfLVi = CfLVUDi;
        // surfaceScalarField CfLVLi = CfLVLUDi;

    /**********************************************************************/
        // functions that interpolate in the bulk
        // The volume averages of gas and liquid concentrations
        bulkInterpolate(CLVi, interfaceCells_0, CfLVVLi, CfLVi);
        bulkInterpolate(CGVi, interfaceCells_0, CfGVVLi, CfGVi);
        // The volume averages of gas and liquid concentrations
        // bulkInterpolate(CLVLi, interfaceCells_0, CfLVLVLi, CfLVLi);

        // Compute the diffusion coefficient in the bulk
        bulkDiffusivity(alphaOld, interfaceCells_0, DLi, DLeffi);
        bulkDiffusivity(alpha2Old, interfaceCells_0, DGi, DGeffi);

        // Compute the diffusion coefficient in the right phase at the interface
        interfaceDiffusivity(DLi, DLInt, interfaceCells_0);
        interfaceDiffusivity(DGi, DGInt, interfaceCells_0);

        // Interpolates in the interface cells
        // The volume averages of gas and liquid concentrations
        interfaceInterpolate(CLVi, alphaOld, interfaceCells_0, alphaf, alphaPhi, CfLVi);
        interfaceInterpolate(CGVi, alpha2Old, interfaceCells_0, alpha2f, alpha2Phi, CfGVi);
        // The volume averages of gas and liquid concentrations
        // phaseAverageInterpolate(CLVLi, alphaOld, interfaceCells_0, alphaf, alphaPhi, CfLVLi);
        // interfaceInterpolate(CGVGi, alpha2Old, interfaceCells_0, alpha2f, alpha2Phi, CfGVi);

        // Parallel exchanges between processors
        syncProcPatches(CfLVi, interfaceCells_0);
        syncProcPatches(CfGVi, interfaceCells_0);
        // limitFaceValues(CfLVi);
        // limitFaceValues(CfGVi);
        /*************************** Advection ****************************/
        // CLVi.oldTime() = CLVi.oldTime() - fvc::surfaceIntegrate(phi*CfLVi) * mesh.time().deltaT();
        // CGVi.oldTime() = CGVi.oldTime() - fvc::surfaceIntegrate(phi*CfGVi) * mesh.time().deltaT();
        // Create a tmp field to do the advection in order to keep only the desired fluxes
        volScalarField CLViOldTmp = CLVi.oldTime();
        volScalarField CLVLiOldTmp = CLVLi.oldTime();
        volScalarField CGViOldTmp = CGVi.oldTime();
        volScalarField CGVGiOldTmp = CGVGi.oldTime();
        // advect these tmp fields
        CLViOldTmp = CLVi.oldTime() - fvc::surfaceIntegrate(phi*CfLVi) * mesh.time().deltaT();
        // CGVi.oldTime() = alpha2 * rho2;
        // CGViOldTmp = CGVi.oldTime();// / alpha2;
        // CGVi.oldTime() = CGVi.oldTime() - fvc::surfaceIntegrate(phi*CfGVi) * mesh.time().deltaT();
        CGViOldTmp = CGVi.oldTime() - fvc::surfaceIntegrate(phi*CfGVi) * mesh.time().deltaT();
        // CGViOldTmp = rho2 * alpha2;

        // Temporary liquid phase average computation
        forAll(CLVLiOldTmp, cellI)
        {
            if (alpha1[cellI] == 1)
            {
                CLVLiOldTmp[cellI] = CLViOldTmp[cellI];
            }
            else if (mag(alpha1[cellI]) < 1e-012)
            {
                CLVLiOldTmp[cellI] = 0;
            }
            else
            {
                CLVLiOldTmp[cellI] = CLViOldTmp[cellI] / alpha1[cellI];
            }
        }
        // Temporary gas phase average computation
        forAll(CGViOldTmp, cellI)
        {
            if (alpha2[cellI] == 1)
            {
                CGVGiOldTmp[cellI] = CGViOldTmp[cellI];
            }
            else if (mag(alpha2[cellI]) < 1e-012)
            {
                CGVGiOldTmp[cellI] = 0;
            }
            else
            {
                CGVGiOldTmp[cellI] = CGViOldTmp[cellI] / (alpha2[cellI]);
            }
        }
        CGVGi.oldTime() = CGVGiOldTmp;
        // Compute the available space for species i in the liquid phase
        // volScalarField saturationLiquid = degreeOfFilling(CLVLiOldTmp, CGVGiOldTmp, Hi);
        // saturationLiquid should be bounded between 0 and 1
        volScalarField saturationLiquid = CLVLiOldTmp / CLsati;
        // Compute the available space for species i in the gas phase
        // volScalarField saturationGas = degreeOfFilling(CGVGiOldTmp,nCLVLiOldTmp, Hi);
        surfaceScalarField advectionError
        (
            IOobject
            (
                "advectionError",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("advectionError", dimVolume / dimVolume, 1)
        );
        advectionError = alphaf * advectionError;
        // Compute the product Saturation x alphaf on faces for the liquid phase, it must be bounded between 0 and 1
        phaseAverageInterpolate(saturationLiquid, alphaOld, interfaceCells_0, alphaf, alphaPhi, advectionError);
        syncProcPatches(advectionError, interfaceCells_0);
        // redistribute the over / undershoots in neighbouring cells (to be implemented)
        surfaceScalarField advectionCorrection
        (
            IOobject
            (
                "advectionCorrection",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("advectionCorrection", dimVolume / dimVolume, 1)
        );
        forAll(advectionError, faceI)
        {
            if (advectionError[faceI] > mag(alphaf[faceI] - 1e-08))
            {
                advectionCorrection[faceI] = alphaf[faceI] / advectionError[faceI];
            }
            else if (advectionError[faceI] < 0)
            {
                advectionCorrection[faceI] = 0;
            }
            else
            {
                advectionCorrection[faceI] = 1;
            }
        
        //     if (advectionError[faceI] > (1 - 1e-08))
        //     {
        //         advectionCorrection[faceI] = 1 / advectionError[faceI];
        //     }
        //     else if (advectionError[faceI] < 0)
        //     {
        //         advectionCorrection[faceI] = 0;
        //     }
        //     else
        //     {
        //         advectionCorrection[faceI] = 1;
        //     }
        }
        surfaceScalarField CfLViTmp = CfLVi * advectionCorrection;
        // Do the actual advection
        CLVi.oldTime() = CLVi.oldTime() - fvc::surfaceIntegrate(phi * CfLViTmp) * mesh.time().deltaT();
        // volScalarField CLVOldi = CLVi.oldTime();
        // volScalarField CGVOldi = CGVi.oldTime();
        // Final computation of the liquid phase average 
        forAll(CLVLi.oldTime(), cellI)
        {
            if (alpha1[cellI] == 1)
            {
                CLVLi.oldTime()[cellI] = CLVi.oldTime()[cellI];
            }
            else if (mag(alpha1[cellI]) < 1e-08)
            {
                CLVLi.oldTime()[cellI] = 0;
            }
            else
            {
                CLVLi.oldTime()[cellI] = CLVi.oldTime()[cellI] / alpha1[cellI];
            }
        }
        forAll(saturationLiquid, cellI)
        {
            if (saturationLiquid[cellI] > mag(1 - 1e-08))
            {
                CLVLi.oldTime()[cellI] = /*saturationLiquid[cellI] */ CLsati.value();
            }
        }
        /******************************************************************/
        // Compute the diffusion coefficient in the bulk
        volScalarField CLVLOldi = CLVLi.oldTime();
        volScalarField CGVGOldi = CGVGi.oldTime();
        // Cut off the diffusivity regarding the emptiness of the cell
        correctDiffusivity(alphaOld, alpha1, alphaf, CLVLOldi, CGVGOldi, DLi, DLeffi, interfaceCells_0, Hi);
        correctDiffusivity(alpha2Old, alpha2, alpha2f, CGVGOldi, CLVLOldi, DGi, DGeffi, interfaceCells_0, Hi);
        /********************** Chemical reaction *************************/
        /***************** Diffusion + Chemical reaction ******************/

        CLVLi.correctBoundaryConditions();
        fvScalarMatrix CLVEqni
        (
            fvm::ddt(CLVLi)
          - fvm::laplacian(DLeffi, CLVLi)
        );

        CGVi.correctBoundaryConditions();
        fvScalarMatrix CGEqni
        (
            fvm::ddt(CGVGi)
          - fvm::laplacian(DGeffi,CGVGi)
        );

        CLVEqni.relax();
        CLVEqni.solve();
	    // CGEqni.relax();
        // CGEqni.solve();
        CLVi = CLVLi * alpha1;
        CGVi = CGVGi * alpha2;
 

        /*********************** Total mass flux **************************/

    }
}