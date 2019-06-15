#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <dune/common/parametertreeparser.hh>

#include "bc_extension.hh"
#include "driver.hh"

int main(int argc, char **argv) {

#if !HAVE_SUPERLU
  std::cerr << "Error: These examples work only if SuperLU is available."
            << std::endl;
  exit(1);
#endif
  Dune::MPIHelper::instance(argc, argv);
  // Initialize Navier Stokes parameter class from file
  Dune::ParameterTree params;

  std::string input_filename("stokes.input");
  if (argc > 2)
    input_filename = argv[2];

  std::cout << "Reading input file \"" << input_filename << "\"" << std::endl;

  try {
    Dune::ParameterTreeParser::readINITree(input_filename, params);
  } catch (...) {
    std::cerr << "The configuration file \"" << argv[2] << "\" "
                 "could not be read. Exiting..."
              << std::endl;
    exit(1);
  }

std::string domain_filename("1.msh");
std::string domain_filename_noext = "1";

if (argc > 1) {
	
	domain_filename_noext = argv[1];
	domain_filename = argv[1];
	domain_filename.append(".msh");
	
}

  // 1. Konstrukcija mreže. Koristimo  UG Grid  u 2D
  const int dim = 2;
  typedef Dune::UGGrid<dim> GridType;
  typedef double RF;

  // Pročitaj sve parametre iz ulazne datoteke
  int refinement   = params.get<int>("domain.level");
  std::string outf = params.get<std::string>("OutputFile"); // izlazni VTK file
  RF mu            = params.get<double>("physics.mu");
  RF rho           = params.get<double>("physics.rho");
  RF intensity      = params.get<double>("velmax");  // maksimalna ulazna brzina

  // Mreža se čita iz grids/lshape.msh datoteke
  std::unique_ptr<GridType> pgrid{Dune::GmshReader<GridType>::read(domain_filename, true, false)};
  pgrid->globalRefine(refinement);

  using GV =  GridType::LeafGridView;
  const GV &gv = pgrid->leafGridView();

  // 2. Konstrukcija klase paramatara
  // Dune::PDELab::NavierStokesDefaultParameters.
  // Ovdje moramo konstruirati Dirichletov rubni uvjet jer ga očekuje klasa
  // Dune::PDELab::NavierStokesDefaultParameters<>
  using BdryPressure = ZeroFunction<GV, RF, 1>; // dummy r.u. // TU KLASU MI PISEMO
  using BdryVelocity = Velocity<GV, RF, 2>; // TU KLASU MI PISEMO
  using BdrySolution =  Dune::PDELab::CompositeGridFunction<BdryVelocity, BdryPressure>; // tenzorski produkt

  BdryVelocity bdry_velocity(gv, intensity);
  BdryPressure bdry_pressure(gv);
  BdrySolution bdry_solution(bdry_velocity, bdry_pressure);

  BCTypeParam bdry_type;

  using SourceFunction = ZeroFunction<GV, RF, dim>;
  using NeumannFlux    = ZeroFunction<GV, RF, dim>; // artificijelna definicija. Nemamo Neumannovog rubnog uvjeta

  NeumannFlux    neumann_flux(gv); // Ne koristi se ali mora biti prisutan
  SourceFunction source_function(gv);

 // const bool navier = true; // Treba li asemblirati nelinearni konvektivni član ?
  const bool navier = true;
  const bool tensor = true; // Treba li raditi sa simetriziranim gradijentom ili ne
                            // (utječe na interpretaciju rubnih uvjeta)

  using Parameters = Dune::PDELab::NavierStokesDefaultParameters<
      GV, RF, SourceFunction, BCTypeParam, BdrySolution, NeumannFlux, navier, tensor>;

  Parameters parameters(mu, rho, source_function, bdry_type, bdry_solution,
                           neumann_flux);
	outf = domain_filename_noext; 
  //driver<GV, BdrySolution, Parameters>(gv, outf, parameters, bdry_solution);
  driver(gv, outf, parameters, bdry_solution); // DOVOLJNO OVAKO

  return 0;
}
