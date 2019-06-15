/* 
 * File:   driver.hh
 * 
 */

#ifndef DRIVER_HH
#define	DRIVER_HH

#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>

#include <dune/pdelab/localoperator/taylorhoodnavierstokes.hh>
#include <dune/pdelab/newton/newton.hh>

template<typename GV, typename IF, typename PARAMS>
void driver(const GV& gv, std::string filename, PARAMS & parameters, IF & bdry_solution)
{  
    using namespace Dune::PDELab;  // skrati imena

    static const unsigned int dim = GV::dimension;

    // Konstrukcija prostora konačnih elemenata
    const int k = 2;
    const int q = 2 * k;  // preciznost integracijske formule lokalnog operatora

    // Taylor-Hoodovi elementi -- P2 za brzinu P1 za tlak 
    typedef PkLocalFiniteElementMap<GV, double, double, k>      V_FEM;  // komponenta brzine
    typedef PkLocalFiniteElementMap<GV, double, double, k - 1 > P_FEM;  // tlak
    V_FEM vFem(gv);
    P_FEM pFem(gv);
    // za brzinu polinomi 2 stupnja (P2), a za tlak 1 stupnja (P1)


    using CDC = ConformingDirichletConstraints;
    using VB = ISTL::VectorBackend<>;
    // Ova klasa direktno konstruira vektorske elemente u R^dim. 
    // Prostor mrežnih funkcija za brzinu (vektorski):  V_h 
    using V_GFS = VectorGridFunctionSpace<GV, V_FEM, dim, VB, VB, CDC>;
    V_GFS VGfs(gv, vFem);
    VGfs.name("velocity");
    // Prostor mrežnih funkcija za tlak (skalarni): W_h
    using P_GFS = GridFunctionSpace<GV, P_FEM, CDC, VB>;
    P_GFS pGfs(gv, pFem);
    pGfs.name("pressure");
    // Prostor V_h x W_h 
    // LexicographicOrderingTag daje poredak varijabli: v1,v2,p
    // KARTEZIJEV PRODUKT
    using GFS = CompositeGridFunctionSpace<VB, LexicographicOrderingTag, V_GFS, P_GFS> ;
    GFS gfs(VGfs, pGfs);

    // Primjena Dirichletovih ograničenja
    using C = typename GFS::template ConstraintsContainer<double>::Type;
    C cg;
    cg.clear();

    // Određivanje Dirichletove granice. Ovdje se koriste pomoćne klase koje određuju
    // Dirichletovu granicu za svaku komponentu vektorske funkcije (v1,v2,p).
    using ScalarVelConstraints = StokesVelocityDirichletConstraints<PARAMS>;
    using VelocityConstraints = PowerConstraintsParameters<ScalarVelConstraints, dim>;
    using PressureConstraints = StokesPressureDirichletConstraints<PARAMS>;
    using Constraints = CompositeConstraintsParameters<VelocityConstraints, PressureConstraints>;

    ScalarVelConstraints scalarvelocity_constraints(parameters);
    VelocityConstraints  velocity_constraints(scalarvelocity_constraints);
    PressureConstraints  pressure_constraints(parameters);
    Constraints          bconst(velocity_constraints, pressure_constraints); // odredi koji su dirichletovi

    // Odredi Dirichletova ograničenja
    constraints(bconst, gfs, cg);

    // Prostorni lokalni operator - definiran u Dune::PDELab-u.
    using LOP = TaylorHoodNavierStokes<PARAMS>;
    LOP lop(parameters, q);

    // Mrežni operator
    using MBE = ISTL::BCRSMatrixBackend<>;
    MBE mbe(5); // maksimalan broj ne-nul elemenat u retku (samo pretpostavka)
    using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, double, double, double, C, C>;
    GO go(gfs, cg, gfs, cg, lop, mbe);

    // Vektor koeficijenata i interpolacija rubnog uvjeta
    using U = typename GO::Traits::Domain;
    U x0(gfs);
    x0 = 0.0;
    interpolate(bdry_solution, gfs, x0);

    // Postavi sve stupnjeve slobode koji nisu Dirichletovi na nulu.
    set_shifted_dofs(cg, 0.0, x0);

    // Linear solver
    using LS = ISTLBackend_SEQ_SuperLU;
    LS ls(false);

    // Riješi sustav.
    Newton<GO, LS, U> newton(go, x0, ls);
    newton.setReassembleThreshold(0.0);
    newton.setVerbosityLevel(2);
    newton.setMaxIterations(25);
    newton.setLineSearchMaxIterations(30);
    newton.apply();

    // Izračunaj rezidual i ispiši njegovu normu.
    U r(gfs);
    r = 0.;
    go.residual(x0, r);
    std::cout << "Final Residual: " << r.two_norm() << std::endl;

    // Ispis rješenja.
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, Dune::RefinementIntervals{k});
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x0); // damo samo koeficijente a ona doda komponente
    vtkwriter.write(filename, Dune::VTK::ascii);  
}

#endif	/* DRIVER_HH */

