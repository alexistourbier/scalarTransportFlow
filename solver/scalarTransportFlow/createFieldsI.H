#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

#include "createSpecies.H"

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

/**************************************************************/

surfaceScalarField alphaf
(
    IOobject
    (
        "alphaf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0*fvc::interpolate(alpha1)
);

surfaceScalarField alpha2f
(
    IOobject
    (
        "alpha2f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1 - alphaf
);

surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

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

/**************************************************************/

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        U.time().constant(),
        U.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);



surfaceForces surfForces(alpha1,phi,U,transportProperties);


//#include "readGravitationalAcceleration.H"
//#include "readhRef.H"
//#include "gh.H"
const volScalarField& gh = surfForces.acc();



volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"

autoPtr<advectionSchemes> advector
(
    advectionSchemes::New(alpha1,phi,U)
);

volVectorField interfaceNormals_0
(
    IOobject
    (
        "interfaceNormals_0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector->normal()
);

volVectorField interfaceNormals
(
    IOobject
    (
        "interfaceNormals",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector->normal()
);

volVectorField interfaceCenters_0
(
    IOobject
    (
        "interfaceCenters_0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector->centre()
);

volVectorField interfaceCenters
(
    IOobject
    (
        "interfaceCenters",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    advector->centre()
);

const dictionary& advectionDict = advector->surf().modelDict();
autoPtr<reconstructionSchemes> PLIC
(
    reconstructionSchemes::New(alpha1,phi,U,advectionDict)
);

