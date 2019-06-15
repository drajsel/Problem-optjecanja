#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <dune/common/parallel/mpihelper.hh> 
#include <dune/grid/uggrid.hh>             // Koristimo UGGrid
#include <dune/grid/io/file/gmshreader.hh> // GmshReader klasa
#include <dune/grid/common/gridinfo.hh>

#include "driver.hh"

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    
    using GridType = Dune::UGGrid<2	>;
    using GridView = GridType::LeafGridView;

    bool verbosity = true;
    bool insertBoundarySegments = false; //Bez toga GmshReader zna podbaciti (barem u 3D)
    std::vector<int> boundarySegmentToPhysicalEntity;  // To ne koristimo
    // Vektor u kojem je za svaki element označeno kojoj poddomeni pripada.
    std::vector<int> materijali;

std::string domain_filename("1.msh");
std::string domain_filename_noext = "1";

if (argc > 1) {
	
	domain_filename_noext = argv[1];
	domain_filename = argv[1];
	domain_filename.append(".msh");
	
}

    GridType* pgrid = Dune::GmshReader<GridType>::read(domain_filename, boundarySegmentToPhysicalEntity,
                                                       materijali, verbosity, insertBoundarySegments);
    
    // loadBalance() distribuira mrežu po procesorima
    pgrid->loadBalance();

    Dune::gridinfo(*pgrid);

    auto gv = pgrid->leafGridView();

	std::string outf = domain_filename_noext; 
	
    driver(gv, outf);
    
    return 0;
}
