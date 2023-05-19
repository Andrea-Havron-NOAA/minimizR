

#include <Rcpp.h>
#include <valarray>
#include <chrono>

const double mgc(const std::valarray<double>& gr){
  double maxgc = std::numeric_limits<double>::min();
  for(int i = 0; i < gr.size(); i++){
    if(std::fabs(gr[i]) >= maxgc){
      maxgc = std::fabs(gr[i]);
    }
  }
  return maxgc;
}

const double norm(const std::valarray<double> &v) {

  double ret = 0.0;
  unsigned int i;
  for (i = 0; i < v.size(); i++) {

    ret += v[i] * v[i];

  }
  return std::sqrt(ret);
}

const double norm(const Rcpp::NumericVector &v) {

  double ret = 0.0;
  unsigned int i;
  for (i = 0; i < v.size(); i++) {

    ret += v[i] * v[i];

  }
  return std::sqrt(ret);
}

Rcpp::NumericVector to_nvector(const std::valarray<double> &v){
  Rcpp::NumericVector ret(v.size());
  for(int i =0; i < v.size(); i++){
    ret[i] = v[i];
  }
  return ret;
}

std::valarray<double> from_nvector(const Rcpp::NumericVector &v){
  std::valarray<double> ret(v.size());
  for(int i =0; i < v.size(); i++){
    ret[i] = v[i];
  }
  return ret;
}

/**
 * Compute the dot product of two vectors.
 * @param a
 * @param b
 * @return
 */
const double Dot(const std::valarray<double> &a, const std::valarray<double> &b) {
  double ret = 0;
  for (size_t i = 0; i < a.size(); i++) {

    ret += a[i] * b[i];
  }
  return ret;
}


/**
 * returns the a column of a matrix as a std::valarray.
 * @param matrix
 * @param column
 * @return
 */
const std::valarray<double> Column(std::valarray<std::valarray<double> > &matrix, size_t column, size_t length) {

  std::valarray<double> ret(length);

  for (int i = 0; i < ret.size(); i++) {

    ret[i] = matrix[i][column];
  }
  return ret;
}


bool line_search(     Rcpp::Function fn,
                      Rcpp::Function  gr,
                      double& fx,
                      double& function_value,
                      std::valarray<double>& x,
                      std::valarray<double>& best,
                      std::valarray<double>& z,
                      std::valarray<double>& gradient,
                      std::valarray<double>& wg,
                      double& maxgc, int& i,
                      int& max_iterations,
                      bool inner = true) {

  //    int max_iterations = 1000;
  double tolerance = 1e-4;
  int max_line_searches = 1000;
  double descent = 0;

  int nops = x.size();

  std::valarray<double> ng(nops);

  for (size_t j = 0; j < nops; j++) {
    descent += z[j] * wg[j];
  }//end for

  double norm_g = norm(gradient);
  double relative_tolerance = tolerance * std::max<double > ((1.0), norm_g);

  descent *= -1.0; // * Dot(z, g);
  if ((descent > (-0.00000001) * relative_tolerance /* tolerance relative_tolerance*/)) {
    z = wg + .001;
    if (!inner) {
      max_iterations -= i;
      i = 0;
    }
    descent = -1.0 * Dot(z, wg);
  }//end if

  double step = i ? 1.0 : (1.0 / norm_g);

  if (step != step || step == 0.0) {
    step = 1.0;
  }

  bool down = false;

  int ls;




  for (int j = 0; j < nops; j++) {
    best[j] = x[j];
  }


  for (ls = 0; ls < max_line_searches; ++ls) {

    // Tentative solution, gradient and loss
    std::valarray<double> nx = x - step * z;



    //line_search:
    fx = Rcpp::as<double>(fn(to_nvector(nx)));

    if (fx <= function_value + tolerance * (10e-4) * step * descent) { // First Wolfe condition

      for (size_t j = 0; j < nops; j++) {
        best[j] = nx[j];
      }





      ng = from_nvector(gr(to_nvector(nx)));
      maxgc = mgc(ng);

      if (down || (-1.0 * Dot(z, ng) >= 0.9 * descent)) { // Second Wolfe condition
        x = nx;
        gradient = ng;
        function_value = fx;
        return true;
      } else {

        step *= 2.0; //*= 10.0; //2.0; //10.0;
      }
    } else {
      step *= .5; /// /= 10.0; //*= .5; ///
      down = true;
    }
  }

  for (size_t j = 0; j < nops; j++) {
    x[j] = best[j];
  }
  fx = Rcpp::as<double>(fn(to_nvector(x)));

  return false;

}


//l-bfgs minimizer

// [[Rcpp::export]]

Rcpp::List minimize(
    Rcpp::NumericVector par,  Rcpp::Function fn,
    Rcpp::Function gr,
    Rcpp::Nullable<Rcpp::List> control = R_NilValue) {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  int max_history = 200;
  int max_iterations = 1000;
  double tolerance = 1e-4;
  double function_value;
  double maxgc;
  bool verbose = true;
  int iprint = 10;
  bool converged = false;
  Rcpp::List results;


  Rcpp::List ctrl(control);
  if (control.isNotNull()) {
    if (ctrl.containsElementNamed("max_iterations")) {
      double maxi = ctrl["max_iterations"];
      if (maxi != 0) {
        max_iterations = maxi;
      }
    }

    if (ctrl.containsElementNamed("tolerance")) {
      double tol = ctrl["tolerance"];
      tolerance = tol;
    }

    if (ctrl.containsElementNamed("iprint")) {
      int print_interval = ctrl["iprint"];
      if (print_interval != 0) {
        iprint = print_interval;
      }
    }
    if (ctrl.containsElementNamed("verbose")) {
      bool is_verbose = ctrl["verbose"];
      if (!is_verbose) {
        verbose = is_verbose;
      }
    }
  }


  int nops = par.size();

  std::valarray<double> x, best, gradient;

  x.resize(nops);
  best.resize(nops);
  gradient.resize(nops);

  for (int i = 0; i < nops; i++) {
    x[i] = par[i];
    gradient[i] = 0;
  }


  std::valarray<double> wg(nops);
  std::valarray<double> nwg(nops);
  std::valarray<double> ng(nops);
  Rcpp::NumericVector rgrad(nops);
  Rcpp::NumericVector rx(nops);

  //initial evaluation
  double fx(0.0);
  fx = Rcpp::as<double>(fn(to_nvector(x)));//internal_evaluate(x);
  function_value = fx;

  //Historical evaluations
  std::valarray<double> px(nops);
  std::valarray<double> pg(nops);
  std::valarray<std::valarray<double> > dxs(std::valarray<double > (max_history), nops);
  std::valarray<std::valarray<double> > dgs(std::valarray<double > (max_history), nops);
  //search direction
  std::valarray<double> z(nops);



  gradient = from_nvector(gr(to_nvector(x)));
  maxgc = mgc(gradient);



  std::valarray<double> p(max_history);
  std::valarray<double> a(max_history);
  int no_progress_count = 0;
  int i;
  for (int iteration = 0; iteration < max_iterations; iteration++) {
    i = iteration;

    for (int j = 0; j < nops; j++) {
      wg[j] = gradient[j];
    }

    if (((i % iprint) == 0) && verbose) {
      std::cout << "Iteration " << i << "\n";
      std::cout << "f = " << fx << ", maxgc = " << maxgc << "\n";
    }


    if (maxgc < tolerance) {
      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      for (size_t j = 0; j < nops; j++) {
        rgrad[j] = gradient[j];
        rx[j] = x[j];
      }
      results["method"] = "l-bfgs";
      results["converged"] = true;
      results["message"] = "NA";
      results["iterations"] = iteration;
      results["runtime (seconds)"] = elapsed_seconds.count();
      results["function value"] = fx;
      results["norm gradient"] = norm(rgrad);
      results["max gradient component"] = maxgc;
      results["gradient"] = rgrad;
      results["parameter values"] = rx;

      return results;
    }

    z = wg;

    if (i > 0 && max_history > 0) {

      size_t h = std::min<size_t > (i, max_history);
      size_t end = (i - 1) % h;

      //update histories
      for (size_t r = 0; r < nops; r++) {
        dxs[r][end] = x[r] - px[r];
        dgs[r][end] = wg[r] - pg[r];
      }



      for (size_t j = 0; j < h; ++j) {
        const size_t k = (end - j + h) % h;
        p[k] = 1.0 / Dot(Column(dxs, k, nops), Column(dgs, k, nops));

        a[k] = p[k] * Dot(Column(dxs, k, nops), z);
        z -= a[k] * Column(dgs, k, nops);
      }
      // Scaling of initial Hessian (identity matrix)
      z *= Dot(Column(dxs, end, nops), Column(dgs, end, nops)) / Dot(Column(dgs, end, nops), Column(dgs, end, nops));

      for (size_t j = 0; j < h; ++j) {
        const size_t k = (end + j + 1) % h;
        const double b = p[k] * Dot(Column(dgs, k, nops), z);
        z += Column(dxs, k, nops) * (a[k] - b);
      }

    }//end if(i>0)

    for (size_t j = 0; j < nops; j++) {
      px[j] = x[j];
      //            x[j] = px[j];
      pg[j] = wg[j];


    }//end for





    double fv = function_value;
    if (!line_search( fn,
                      gr,
                      fx,
                      function_value,
                      x,
                      best,
                      z,
                      gradient,
                      wg,
                      maxgc,
                      iteration,
                      max_iterations,
                      false)) {
      std::cout << "Max line searches\n\n";
      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      for (size_t j = 0; j < nops; j++) {
        rgrad[j] = gradient[j];
        rx[j] = x[j];
      }
      maxgc = mgc(gradient);
      if(maxgc<= tolerance){
        converged = true;
      }

      results["mehtod"] ="l-bfgs";
      results["converged"] = converged;
      results["message"] = "max line searches";
      results["iterations"] = iteration;
      results["runtime (seconds)"] = elapsed_seconds.count();
      results["function value"] = fx;
      results["norm gradient"] = norm(rgrad);
      results["max gradient component"] = maxgc;
      results["gradient"] = rgrad;
      results["parameter values"] = rx;

      return results;

    }

    if ((fv - function_value) == 0.0 && no_progress_count == 15) {
      std::cout << "Not progressing...bailing out!\n";
      end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      for (size_t j = 0; j < nops; j++) {
        rgrad[j] = gradient[j];
        rx[j] = x[j];
      }


      if(maxgc<= tolerance){
        converged = true;
      }
      results["method"] = "l-bfgs";
      results["converged"] = converged;
      results["message"] = "no progress";
      results["iterations"] = iteration;
      results["runtime (seconds)"] = elapsed_seconds.count();
      results["function value"] = fx;
      results["norm gradient"] = norm(rgrad);
      results["max gradient component"] = maxgc;
      results["gradient"] = rgrad;
      results["parameter values"] = rx;
      return results;
    } else {
      no_progress_count++;
    }

  }

  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  for (size_t j = 0; j < nops; j++) {
    rgrad[j] = gradient[j];
    rx[j] = x[j];
  }

  if(maxgc<= tolerance){
    converged = true;
  }

  results["method"] ="l-bfgs";
  results["converged"] = converged;
  results["message"] = "max iterations";
  results["iterations"] = i;
  results["runtime (seconds)"] = elapsed_seconds.count();
  results["function value"] = fx;
  results["norm gradient"] = norm(rgrad);
  results["max gradient component"] = maxgc;
  results["gradient"] = rgrad;
  results["parameter values"] = rx;

  std::cout << "Max iterations!\n\n";

  return results;
}

