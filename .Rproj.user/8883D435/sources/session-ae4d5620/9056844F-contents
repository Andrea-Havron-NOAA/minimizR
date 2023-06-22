

#include <Rcpp.h>
#include <valarray>
#include <chrono>
#include <sstream>

const double mgc(const std::valarray<double>& gr) {
    double maxgc = std::numeric_limits<double>::min();
    for (int i = 0; i < gr.size(); i++) {
        if (std::fabs(gr[i]) >= maxgc) {
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

Rcpp::NumericVector to_nvector(const std::valarray<double> &v) {
    Rcpp::NumericVector ret(v.size());
    for (int i = 0; i < v.size(); i++) {
        ret[i] = v[i];
    }
    return ret;
}

std::valarray<double> from_nvector(const Rcpp::NumericVector &v) {
    std::valarray<double> ret(v.size());
    for (int i = 0; i < v.size(); i++) {
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
const double Dot(const std::valarray<double> &a,
        const std::valarray<double> &b) {
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
const std::valarray<double> Column(
        std::valarray<std::valarray<double> > &matrix,
        size_t column,
        size_t length) {

    std::valarray<double> ret(length);

    for (int i = 0; i < ret.size(); i++) {

        ret[i] = matrix[i][column];
    }
    return ret;
}



//transformations for value bounded parameters

/**
 * Sine transformation for internal to external.
 */
inline double transform_internal_2_external(
        const double& value,
        const double&minb,
        const double&maxb) {

    return minb + (std::sin(value) + 1.0)*((maxb - minb) / 2.0);

    //    return minb + .5 * (maxb - minb)*(1.0 + std::tanh(value));
}

/**
 * Sine transformation for external to internal.
 */
inline double transform_external_2_internal(
        const double& value,
        const double&minb,
        const double&maxb) {
    return std::asin((2.0 * (value - minb) / (maxb - minb)) - 1.0);
    //    return std::atanh(2.0 * (value - minb) / (maxb - minb) - 1.0);
}

/**
 * Sine transformation for internal derivative to external.
 */
inline double transform_derivative_internal_2_external(
        const double& value,
        const double&minb,
        const double&maxb) {

    return 0.5 * ((maxb - minb) * std::cos(value));
    //    return 2.0 / ((maxb - minb) * std::pow((1.0 - ((2.0 * (value - minb)) / maxb - minb - 1.0)), 2.0));
}

//void find_r_xmin_xmax(double& xmin,
//        double& xmax) {
//
//    Rcpp::Environment env = Rcpp::Environment::global_env();
//    Rcpp::List l = Rcpp::as<Rcpp::List>(env.ls(true));
//    SEXP e, E, EE;
//    std::stringstream ss;
//    ss << "capture.output(show(.Machine$double.xmax))"; //, file = NULL, append = FALSE, type = c(\"output\", \"message\"), split = FALSE)";
//
//
//    SEXP expression, result;
//    ParseStatus status;
//
//    PROTECT(expression = R_ParseVector(Rf_mkString(ss.str().c_str()), 1, &status, R_NilValue));
//
//    if (status != PARSE_OK) {
//        std::cout << "Error parsing expression" << std::endl;
//        UNPROTECT(1);
//    }
//
//    PROTECT(result = Rf_eval(VECTOR_ELT(expression, 0), R_GlobalEnv));
//
//    double xmax_r  = Rcpp::as<double>(result);
//
//
//}

void update_value(
        double& parameter,
        double& value,
        const double& minb,
        const double& maxb,
        const bool& bounded) {

    if (bounded) {
        parameter = transform_internal_2_external(value, minb, maxb);
    } else {
        parameter = value;
    }
}

void get_working_derivative(double& value,
        const double& minb,
        const double& maxb,
        const bool& bounded) {

}

Rcpp::NumericMatrix calculateHessian(Rcpp::Function gradFunc, Rcpp::NumericVector x, double h = 1e-4) {
    int n = x.size();
    Rcpp::NumericMatrix hessian(n, n);

    for (int i = 0; i < n; i++) {
        Rcpp::NumericVector x_plus_h = clone(x);
        Rcpp::NumericVector x_minus_h = clone(x);
        double xi = x[i];

        x_plus_h[i] = xi + h;
        x_minus_h[i] = xi - h;

        Rcpp::NumericVector grad_plus_h = gradFunc(x_plus_h);
        Rcpp::NumericVector grad_minus_h = gradFunc(x_minus_h);

        for (int j = 0; j < n; j++) {
            double hessian_val = (grad_plus_h[j] - grad_minus_h[j]) / (2 * h);
            hessian(i, j) = hessian_val;
        }
    }

    return hessian;
}

bool line_search(Rcpp::Function fn,
        Rcpp::Function gr,
        double& fx,
        double& function_value,
        std::valarray<double>& x,
        std::valarray<double>& wx,
        Rcpp::NumericVector& minb,
        Rcpp::NumericVector& maxb,
        std::vector<bool> bounded_parameter,
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
        std::valarray<double> nx = x - step * z; //(nops); // = x - step * z;

//                for (size_t j = 0; j < nops; j++) {
//                    nx[j] = wx[j] - step * z[j];
//                    double temp = nx[j];
//                    if (bounded_parameter[j]) {
//                        nx[j] = transform_internal_2_external(temp, minb[j], maxb[j]);
//                    } else {
//                        nx[j] = temp;
//                    }
//                }


        //line_search:
        fx = Rcpp::as<double>(fn(to_nvector(nx)));

        if (fx <= function_value + tolerance * (10e-4) * step * descent) { // First Wolfe condition

            for (size_t j = 0; j < nops; j++) {
                if (bounded_parameter[j]) {
                    best[j] = transform_internal_2_external(nx[j], minb[j], maxb[j]);
                } else {
                    best[j] = nx[j];
                }
            }





            ng = from_nvector(gr(to_nvector(nx)));
            maxgc = mgc(ng);

            if (down || (-1.0 * Dot(z, ng) >= 0.9 * descent)) { // Second Wolfe condition
                x = nx;
                for (size_t j = 0; j < nops; j++) {

                    if (bounded_parameter[j]) {
                        wx[j] = transform_external_2_internal(nx[j], minb[j], maxb[j]);
                    } else {
                        wx[j] = nx[j];
                    }
                }

                gradient = ng;
                function_value = fx;
                return true;
            } else {

                step *= 10.0; //2.0; //10.0;
            }
        } else {
            step /= 10.0; //*= .5; ///
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

Rcpp::List minimizR(
        Rcpp::NumericVector par,
        Rcpp::Function fn,
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
    bool bounded = false;
    bool error = false;
    double xmin = -1.0 * std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::max();
    std::vector<bool> bounded_parameter(par.size(), false);
    Rcpp::List results;
    Rcpp::NumericVector minb; //(par.size(), xmin);
    Rcpp::NumericVector maxb; //(par.size(), xmax);
    bool do_hessian = false;


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

        if (ctrl.containsElementNamed("lb")) {
            bounded = true;
            minb = ctrl["lb"];
        }

        if (ctrl.containsElementNamed("ub")) {
            bounded = true;
            maxb = ctrl["ub"];
        }
        
        if (ctrl.containsElementNamed("hessian")) {
            do_hessian = Rcpp::as<bool>(ctrl["hessian"]);
        }

    }




    //error checking
    std::stringstream error_stream;
    if (bounded) {
        if (minb.size() != par.size()) {
            Rcpp::Rcout << "Error: minb not equal to size of parameter vector.\n";
            error_stream << "Error: minb not equal to size of parameter vector.\n";
            error = true;
        }

        if (maxb.size() != par.size()) {
            Rcpp::Rcout << "Error: maxb not equal to size of parameter vector.\n";
            error_stream << "Error: maxb not equal to size of parameter vector.\n";
            error = true;
        }
    }


    //begin initialization
    int nops = par.size();

    std::valarray<double> x, wx, best, gradient;

    x.resize(nops);
    wx.resize(nops);
    best.resize(nops);
    gradient.resize(nops);


    //if bounded, specify which parameters have numeric bounds
    if (bounded) {
        //        find_r_xmin_xmax(xmin, xmax);

        bounded_parameter.resize(par.size());
        //not all parameters are necessarily bounded.
        std::fill(bounded_parameter.begin(), bounded_parameter.end(), false);

        //find the bounded parameters.
        for (size_t i = 0; i < par.size(); i++) {
            if (minb[i] != xmin || maxb[i] != xmax) {
                bounded_parameter[i] = true;
                std::cout <<"parameter "<<i<<" is bounded!\n";
            }
        }

    }


    for (int i = 0; i < nops; i++) {

        if (bounded_parameter[i]) {

            if (par[i] < minb[i]) {
                error_stream << "Error: parameter[" << i + 1 << "] less than minb[" << i + 1 << "]. "<<std::endl;
                error = true;

            }

            if (par[i] > maxb[i]) {
                error_stream << "Error: parameter[" << i + 1 << "] greater than maxb[" << i + 1 << "]. "<<std::endl;
                error = true;
            }

            wx[i] = transform_external_2_internal(par[i], minb[i], maxb[i]);

        } else {
            wx[i] = par[i];
        }
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
    fx = Rcpp::as<double>(fn(par)); //internal_evaluate(x);
    function_value = fx;

    //Historical evaluations
    std::valarray<double> px(nops);
    std::valarray<double> pg(nops);
    std::valarray<std::valarray<double> > dxs(std::valarray<double > (max_history), nops);
    std::valarray<std::valarray<double> > dgs(std::valarray<double > (max_history), nops);
    //search direction
    std::valarray<double> z(nops);



    //initial gradient
    gradient = from_nvector(gr(par)); //to_nvector(par)));
    maxgc = mgc(gradient);

    if (error) {
        results["method"] = "l-bfgs";
        results["converged"] = false;
        results["bounded problem"] = bounded;
        results["message"] = error_stream.str();
        results["iterations"] = 0;
        results["runtime (seconds)"] = 0;
        results["function value"] = fx;
        results["norm gradient"] = norm(rgrad);
        results["max gradient component"] = maxgc;
        results["gradient"] = rgrad;
        if(do_hessian){
            results["hessian"] = calculateHessian(gr, rx);
        }
        results["parameter values"] = par;

        return results;
    }




    std::valarray<double> p(max_history);
    std::valarray<double> a(max_history);
    int no_progress_count = 0;
    int i;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        i = iteration;

        for (int j = 0; j < nops; j++) {
            if (bounded_parameter[i]) {
                wx[j] = transform_external_2_internal(x[j], minb[j], maxb[j]);
                wg[j] = transform_derivative_internal_2_external(wx[j],
                        minb[j], maxb[j]) *
                        gradient[j];
            } else {
                wg[j] = gradient[j];
                wx[j] = x[j];
            }
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
            results["bounded problem"] = bounded;
            results["message"] = "converged";
            results["iterations"] = iteration;
            results["runtime (seconds)"] = elapsed_seconds.count();
            results["function value"] = fx;
            results["norm gradient"] = norm(rgrad);
            results["max gradient component"] = maxgc;
            results["gradient"] = rgrad;
            if(do_hessian){
                results["hessian"] = calculateHessian(gr, rx);
            }
            results["parameter values"] = rx;

            return results;
        }

        z = wg;

        if (i > 0 && max_history > 0) {

            size_t h = std::min<size_t > (i, max_history);
            size_t end = (i - 1) % h;

            //update histories
            for (size_t r = 0; r < nops; r++) {
                dxs[r][end] = wx[r] - px[r];
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
            px[j] = wx[j];
            //            x[j] = px[j];
            pg[j] = wg[j];


        }//end for




        double fv = function_value;

        if (!line_search(fn,
                gr,
                fx,
                function_value,
                x,
                wx,
                minb,
                maxb,
                bounded_parameter,
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
            if (maxgc <= tolerance) {
                converged = true;
            }

            results["method"] = "l-bfgs";
            results["converged"] = converged;
            results["bounded problem"] = bounded;
            results["message"] = "max line searches";
            results["iterations"] = iteration;
            results["runtime (seconds)"] = elapsed_seconds.count();
            results["function value"] = fx;
            results["norm gradient"] = norm(rgrad);
            results["max gradient component"] = maxgc;
            results["gradient"] = rgrad;
            if(do_hessian){
                results["hessian"] = calculateHessian(gr, rx);
            }
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


            if (maxgc <= tolerance) {
                converged = true;
            }
            results["method"] = "l-bfgs";
            results["converged"] = converged;
            results["bounded problem"] = bounded;
            results["message"] = "no progress";
            results["iterations"] = iteration;
            results["runtime (seconds)"] = elapsed_seconds.count();
            results["function value"] = fx;
            results["norm gradient"] = norm(rgrad);
            results["max gradient component"] = maxgc;
            results["gradient"] = rgrad;
            if(do_hessian){
                results["hessian"] = calculateHessian(gr, rx);
            }
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

    if (maxgc <= tolerance) {
        converged = true;
    }

    results["method"] = "l-bfgs";
    results["converged"] = converged;
    results["bounded problem"] = bounded;
    results["message"] = "max iterations";
    results["iterations"] = i;
    results["runtime (seconds)"] = elapsed_seconds.count();
    results["function value"] = fx;
    results["norm gradient"] = norm(rgrad);
    results["max gradient component"] = maxgc;
    results["gradient"] = rgrad;
    if(do_hessian){
        results["hessian"] = calculateHessian(gr, rx);
    }
    results["parameter values"] = rx;

    std::cout << "Max iterations!\n\n";

    return results;
}
