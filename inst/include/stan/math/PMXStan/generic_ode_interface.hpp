#ifndef __STAN__MATH__FUNCTIONS__GENERIC_ODE_INTERFACE_HPP__
#define __STAN__MATH__FUNCTIONS__GENERIC_ODE_INTERFACE_HPP__

#include <vector>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/PMXStan/lsoda_tmpl_class.hpp>
#include <ctime>


namespace stan {
namespace math {

typedef stan::math::var AVAR;
typedef double ADBL;

static const size_t neq = 2;
static const size_t npar= 2;
static const size_t ncov=0;
static std::vector<AVAR> pars_var(npar);
static std::vector<ADBL> pars_dbl(npar, 0.0);
static std::vector<ADBL> ptr_cov(ncov, 0.0);
static std::vector<double> InfusionRate(neq, 0.0);

void dydt(double t, AVAR *x, AVAR *dxdt, void *data) {
  AVAR
  depot,
  ka,
  centr,
  ke;

  ka = pars_var[0];
  ke = pars_var[1];
  depot = x[0];
  centr = x[1];

  dxdt[0] = InfusionRate[0] + - ka * depot;
  dxdt[1] = InfusionRate[1] + ka * depot - ke * centr;
}

void dydt(double t, ADBL *x, ADBL *dxdt, void *data) {
  ADBL
  depot,
  ka,
  centr,
  ke;

  ka = pars_dbl[0];
  ke = pars_dbl[1];
  depot = x[0];
  centr = x[1];

  dxdt[0] = InfusionRate[0] + - ka * depot;
  dxdt[1] = InfusionRate[1] + ka * depot - ke * centr;
}

void jex(const double tn,
         ADBL  *x,
         Eigen::Matrix<ADBL, Eigen::Dynamic, Eigen::Dynamic>& wm) {
}
void jex(const double tn,
         AVAR  *x,
         Eigen::Matrix<AVAR, Eigen::Dynamic, Eigen::Dynamic>& wm) {
}

template <class T, class T0>
Eigen::Matrix<T,Eigen::Dynamic,1>
generic_ode_interface(Eigen::Matrix<T,Eigen::Dynamic,1>& params,
                      const Eigen::VectorXd& inits,
                      const Eigen::VectorXd& time,
                      const Eigen::VectorXd& evid,
                      const Eigen::VectorXd& amt,
                      const T0& absolute_tolerance,
                      const T0& relative_tolerance,
                      const int nobs,
                      const int which) {
  const char *err_msg[] = {
    "excess work done on this call (perhaps wrong jt).",
    "excess accuracy requested (tolerances too small).",
    "illegal input detected (see printed message).",
    "repeated error test failures (check all inputs).",
    "repeated convergence failures (perhaps bad jacobian supplied or wrong choice of jt or tolerances).",
    "error weight became zero during problem. (solution component i vanished, and atol or atol(i) = 0.)",
    "work space insufficient to finish (see messages)."
  };

#ifdef __TRACK_ELAPSE__
  std::clock_t start;
  double duration;

  start = std::clock();
#endif


  for (int i = 0; i < params.size(); i++) {
	if (std::isinf(params[i]) || params[i] > 1.0E200 || params[i] < -1.0E200) {
	  std::cerr
	    << "WARNING!!!\n"
	    << "extreme parameter value:  "
	    << "parameter " << i+1
	    << ", value = " << params[i]
	    << std::endl;
	}
  }


  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXd;
  VectorXd state(neq+1);
  LSODA<T> ode_solver;

  for (size_t i = 0; i < neq; i++)
  {
    state[i+1] = inits[i];
  }
  for (size_t i = 0; i < npar; i++)
  {
    pars_var[i] = params[i];
    pars_dbl[i] = pars_var[i].val();
  }

  VectorXd C(nobs);
  double t0, t1;
  double rwork1, rwork5, rwork6, rwork7;
  double atol=absolute_tolerance, rtol=relative_tolerance;
  int    iwork1, iwork2, iwork5, iwork6, iwork7, iwork8, iwork9;
  int    itol=1, itask=1, istate=1, iopt=0, jt=2;
  iwork1 = iwork2 = iwork5 = iwork6 = iwork7 = iwork8 = iwork9 = 0;
  rwork1 = rwork5 = rwork6 = rwork7 = 0.0;
  int wh, cmt;

  t0 = time[0];
  for (int i = 0, n = 0; i < time.size(); i++) {
    t1 = time[i];

    if(t1>t0) {
      ode_solver.lsoda(dydt, neq, state, &t0, t1, itol, &rtol-1, &atol-1,
             itask, &istate, iopt, jt, jex, iwork1, iwork2, iwork5, iwork6,
             iwork7, iwork8, iwork9, rwork1, rwork5, rwork6, rwork7, 0);
      if (istate<0) {
        std::cerr << "LSODA exception: " << err_msg[-istate-1] << std::endl;
        std::cerr << params << std::endl;
        std::cerr << time << std::endl;
        std::cerr << amt << std::endl;
        std::cerr << evid << std::endl;
	  }
    }

    wh = evid[i];

    if (wh) {    //dosing events
      cmt = (wh % 10000) / 100 - 1;
      if (wh > 10000)
        InfusionRate[cmt] += amt[i];
      else
        state[cmt+1] += amt[i];    //dosing before obs
      istate = 1;
    }
    else {      //observation events
      C[n] = state[which];
      n++;
    }

    t0 = t1;
  }


#ifdef __TRACK_ELAPSE__
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  if (duration>__TRACK_ELAPSE__) {
    std::cout<<"lsoda duration: " << duration <<'\n';
    std::cout<<"params:\n" << params <<'\n';
  }
#endif

  //ode_solver.n_lsoda_terminate();
  return C;
}


struct lin_cmt_ode {
	const Eigen::VectorXd inits_, time_, evid_, amt_;
	const double absolute_tolerance_, relative_tolerance_;
	const int nobs_, wh_;

	lin_cmt_ode(const Eigen::VectorXd& inits,
			const Eigen::VectorXd& time,
			const Eigen::VectorXd& evid,
			const Eigen::VectorXd& amt,
			const double absolute_tolerance,
			const double relative_tolerance,
			const int nobs,
			const int wh) :
			inits_(inits),
			time_(time),
			evid_(evid),
			amt_(amt),
			absolute_tolerance_(absolute_tolerance),
			relative_tolerance_(relative_tolerance),
			nobs_(nobs),
			wh_(wh)
			{ }



	template <typename T>
	Eigen::Matrix<T, -1, 1> operator()(Eigen::Matrix<T, -1, 1>& params) const {
		return generic_ode_interface<T, double>(
		    params,
		    inits_,
		    time_,
		    evid_,
		    amt_,
		    absolute_tolerance_,
		    relative_tolerance_,
		    nobs_,
		    wh_);
	}
};


} // ns math
}// ns stan
#endif
