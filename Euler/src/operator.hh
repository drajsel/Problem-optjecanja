#ifndef OPERATOR_HH
#define OPERATOR_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

/** Lokalni operator za zadaću :
 *
 *   - div( a(x) grad u) + b(x) u = f(x)   u \Omega
 *                   u = g(x)   na \Gamma_D\subseteq\partial\Omega
 *        -a(x) grad u . n = j(x)   na \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * sa konformnim konačnim elementima svih tipova u svim dimenzijama
 *
 * \tparam BCType klasa koja indicira rubni uvjet
 */



template<typename BCType, typename FEM>
class DiffusionLocalOperator : // derivacijska lista -- jakobijan i pattern računa PDELab
  public Dune::PDELab::NumericalJacobianApplyVolume  <DiffusionLocalOperator<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianVolume       <DiffusionLocalOperator<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<DiffusionLocalOperator<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianBoundary     <DiffusionLocalOperator<BCType, FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // Zastavice koje signaliziraju da na svakom elementu treba zvati:
  enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
  enum { doAlphaVolume = true };    // alpha_volume
  enum { doAlphaBoundary = true };  // alpha_boundary
  using  LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType ;

  DiffusionLocalOperator(const BCType& bctype_, // boundary cond.type
                         const FEM & fem_,
                         unsigned int intorder_=2) :
    bctype( bctype_ ), fem(fem_), intorder( intorder_ )
  {}

  // Računanje volumnog integrala
  // eg   = element (geometry)
  // lfsu = lokalni prostor funkcija za rješenje
  // lfsv = lokalni prostor funkcija za test funkciju
  // x    = vektor koeficijenata rješenja
  // r    = lokalni rezidual
  // int K a(x) grad u . grad phi_i  +  b(x) u phi_i  -   f(x) phi_i dx
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // dimenzije
    const int dim  = EG::Geometry::mydimension;
    const int dimw = EG::Geometry::coorddimension;

    // tipovi
    using  Traits = typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits;
    typedef typename Traits::DomainFieldType DF;
    typedef typename Traits::RangeFieldType  RF;
    typedef typename Traits::JacobianType    Jacobian;
    typedef typename Traits::RangeType       Range;
    typedef Dune::FieldVector<double,dimw> Gradient;
    typedef typename LFSU::Traits::SizeType size_type;

    // integracijska formula
    auto gt = eg.geometry().type();
    auto& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // petlja po svim integracijskim točkama
    for (auto qpoint : rule)
      {
        // računanje baznih funckcija na referentnom elementu
        auto& phi = cache.evaluateFunction(qpoint.position(), lfsu.finiteElement().localBasis());

        // rješenje u integracijskoj točki
        double u=0.0;
        for (size_type i=0; i<lfsu.size(); ++i) u += x(lfsu,i)*phi[i];

        // gradijent baznih funkcija
        auto & gradphihat = cache.evaluateJacobian(qpoint.position(), lfsu.finiteElement().localBasis());

        // transformacija gradijenata s referentnog na fizički element
        auto const & jac = eg.geometry().jacobianInverseTransposed(qpoint.position());
        std::vector<Gradient> gradphi(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
          jac.mv(gradphihat[i][0],gradphi[i]);

        //gradijent rješenja u integracijskoj točki
        Gradient gradu(0.0);
        for (size_type i=0; i<lfsu.size(); ++i)
          gradu.axpy(x(lfsu,i),gradphi[i]); // grad u = sum_i x_i grad phi_i


        auto qglobal = eg.geometry().global(qpoint.position());

        // integriramo :  grad u * grad phi_i + a*u*phi_i - f phi_i
        RF factor = qpoint.weight() * eg.geometry().integrationElement(qpoint.position());

        for (size_type i=0; i<lfsu.size(); ++i)
          r.accumulate(lfsu, i, (gradu*gradphi[i]) * factor);
      }
  }

  // integral po rubu
  // ig     = intersection (= stranica elementa)
  // lfsu_s = lokalni prostor funkcija na stranici za rješenje
  // lfsu_v = lokalni prostor funkcija na stranici za test funkciju
  // x_s    = vektor koeficijenata rješenja (na stranici)
  // r_s    = rezidual (na stranici)
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    // tipovi
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimenzije (u Dune::Geometry imamo mydimension i coorddimension)
    const int dim = IG::Geometry::coorddimension;

    // integracijska formula na stranici
    auto gtface = ig.geometryInInside().type();
    const auto& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // petlja po svim integracijskim točkama
    for (auto qpoint : rule)
      {
        // Ako smo na Dirichletovoj granici preskačemo petlju
        if ( bctype.isDirichlet( ig, qpoint.position() ) )
          continue;

        // pozicija int. točke u lokalnim koordinatam elementa
        auto local = ig.geometryInInside().global(qpoint.position());

        // izračunaj bazne funkcije u integracijskoj točki
//        std::vector<Range> phi(lfsu_s.size());
//        lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);
        auto& phi = cache.evaluateFunction(local, lfsu_s.finiteElement().localBasis());
        // rješenje u integracijskoj točki
        RF u=0.0;
        for (size_type i=0; i<lfsu_s.size(); ++i)
          u += x_s(lfsu_s,i)*phi[i];

        /*// računanje Neumannovog rubnog uvjeta
        Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(qpoint.position());
        RF j = 0.0;
        if (globalpos[1]<0.5)
          j = 1.0;
        else
          j = -1.0;

        // integracija
        RF factor = qpoint.weight()*ig.geometry().integrationElement(qpoint.position());

        for (size_type i=0; i<lfsu_s.size(); ++i)
          r_s.accumulate(lfsu_s,i, j*phi[i]*factor);*/


        //računanje Neumannovog rubnog uvjeta
        Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(qpoint.position());
        RF j = fun_j(globalpos);
        // integracija
        RF factor =  qpoint.weight()*ig.geometry().integrationElement(qpoint.position());

        for (size_type i=0; i<lfsu_s.size(); ++i)
        r_s.accumulate(lfsu_s,i, j*phi[i]*factor);
      }
  }

private:
  BCType const & bctype;
  FEM const & fem;
  unsigned int intorder;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};
#endif
