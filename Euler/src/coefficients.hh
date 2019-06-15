#ifndef COEFFICIENTS
#define COEFFICIENTS

#include <dune/pdelab/common/function.hh>

#include <math.h>

//mijenja se s obzirom na normalu, ovo je za Neumanna na desnoj strani
template <typename Point>
double fun_j(Point const & globalPoint){
    double x = globalPoint[0];
    double val = -100.0;

    if(x < -3 + 1E-6 )
        return val;
    else if(x > 3 - 1E-6 )
        return -val;
    else 
		return 0.0;
}


template <typename GV>
class ExactSolution : public Dune::PDELab::AnalyticGridFunctionBase<
        Dune::PDELab::AnalyticGridFunctionTraits<GV, double, 1>,
        ExactSolution<GV> > {

   public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, double,1> Traits;


    ExactSolution(GV const & gv) :
        Dune::PDELab::AnalyticGridFunctionBase<Traits,ExactSolution<GV> >(gv){}


    void evaluate (const typename Traits::ElementType &e,
                      const typename Traits::DomainType &x,
                      typename Traits::RangeType &y) const
            {
                auto xglobal = e.geometry().global(x);
                y = exact(xglobal);
            }


};


template <typename GV>
class ExactSolutionX : public Dune::PDELab::AnalyticGridFunctionBase<
        Dune::PDELab::AnalyticGridFunctionTraits<GV, double, 1>,
        ExactSolutionX<GV> > {

   public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, double,1> Traits;


    ExactSolutionX(GV const & gv) :
        Dune::PDELab::AnalyticGridFunctionBase<Traits,ExactSolutionX<GV> >(gv){}


    void evaluate (const typename Traits::ElementType &e,
                      const typename Traits::DomainType &x,
                      typename Traits::RangeType &y) const
            {
                auto xglobal = e.geometry().global(x);
                y = grad_exact(xglobal)[0];
            }

};

template <typename GV>
class ExactSolutionY : public Dune::PDELab::AnalyticGridFunctionBase<
        Dune::PDELab::AnalyticGridFunctionTraits<GV, double, 1>,
        ExactSolutionY<GV> > {

   public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, double,1> Traits;


    ExactSolutionY(GV const & gv) :
        Dune::PDELab::AnalyticGridFunctionBase<Traits,ExactSolutionY<GV> >(gv){}


    void evaluate (const typename Traits::ElementType &e,
                      const typename Traits::DomainType &x,
                      typename Traits::RangeType &y) const
            {
                auto xglobal = e.geometry().global(x);
                y = grad_exact(xglobal)[1];
            }

};


template <typename GV>
class ExactSolutionZ : public Dune::PDELab::AnalyticGridFunctionBase<
        Dune::PDELab::AnalyticGridFunctionTraits<GV, double, 1>,
        ExactSolutionZ<GV> > {

   public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, double,1> Traits;


    ExactSolutionZ(GV const & gv) :
        Dune::PDELab::AnalyticGridFunctionBase<Traits,ExactSolutionZ<GV> >(gv){}


    void evaluate (const typename Traits::ElementType &e,
                      const typename Traits::DomainType &x,
                      typename Traits::RangeType &y) const
            {
                auto xglobal = e.geometry().global(x);
                y = grad_exact(xglobal)[1];
            }

};


template <typename GV>
class ExactSolutionP : public Dune::PDELab::AnalyticGridFunctionBase<
        Dune::PDELab::AnalyticGridFunctionTraits<GV, double, 1>,
        ExactSolutionP<GV> > {

   public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, double,1> Traits;


    ExactSolutionP(GV const & gv) :
        Dune::PDELab::AnalyticGridFunctionBase<Traits,ExactSolutionP<GV> >(gv){}


    void evaluate (const typename Traits::ElementType &e,
                      const typename Traits::DomainType &x,
                      typename Traits::RangeType &y) const
            {
                auto xglobal = e.geometry().global(x);
                y = p_exact(xglobal);
            }

};

#endif
