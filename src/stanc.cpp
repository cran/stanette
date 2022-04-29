// This file is part of RStan
// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
//
// RStan is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// RStan is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

// #include <Rcpp.h>
// #include <string>
// #include <iostream>
#include <stan/version.hpp>
#include <stan/lang/compiler.hpp>
#include <stan/lang/compile_functions.hpp>

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Rcpp.h>
#include <rstan/io/r_ostream.hpp>

RcppExport SEXP CPP_stanc280(SEXP model_stancode, SEXP model_name, 
                             SEXP allow_undefined, SEXP include_paths);
RcppExport SEXP CPP_stan_version();

SEXP CPP_stan_version() {
  BEGIN_RCPP
  std::string stan_version
    = stan::MAJOR_VERSION + "." +
      stan::MINOR_VERSION + "." +
      stan::PATCH_VERSION;
  SEXP __sexp_result;
  PROTECT(__sexp_result = Rcpp::wrap(stan_version));
  UNPROTECT(1);
  return __sexp_result;
  END_RCPP
}

/*
 * Split a string by the new line and put each line into a vector.
 */
void split_str_by_newline(const std::string& str,  std::vector<std::string>& v) {
  std::string delim = "\n";
  std::string::size_type start = 0;
  std::string::size_type end = str.find_first_of(delim, start);
  while (start < str.size()) {
      if (end != std::string::npos) {
         v.push_back(str.substr(start, end - start));
      } else {
         v.push_back(str.substr(start));
         break;
      }
     start = end + delim.size();
     end = str.find_first_of(delim, start);
  }
}

SEXP CPP_stanc280(SEXP model_stancode, 
                  SEXP model_name, 
                  SEXP allow_undefined,
                  SEXP include_paths) {
  BEGIN_RCPP;
  static const int SUCCESS_RC = 0;
  static const int EXCEPTION_RC = -1;
  static const int PARSE_FAIL_RC = -2;

  /*
  std::string stan_version
    = stan::MAJOR_VERSION + "." +
      stan::MINOR_VERSION + "." +
      stan::PATCH_VERSION;
  */

  std::string mcode_ = Rcpp::as<std::string>(model_stancode);
  std::string mname_ = Rcpp::as<std::string>(model_name);
  std::vector<std::string> paths_ = Rcpp::as<std::vector<std::string> >(include_paths);

  std::stringstream out;
  std::istringstream in(mcode_);
  try {
    bool valid_model
      = stan::lang::compile(&rstan::io::rcerr,in,out,mname_,
                            Rcpp::as<bool>(allow_undefined),
                            mname_, paths_);
    if (!valid_model) {
      return Rcpp::List::create(Rcpp::Named("status") = PARSE_FAIL_RC);

    }
  } catch(const std::exception& e) {
    // REprintf("\nERROR PARSING\n %s\n", e.what());
    std::string msgstr = e.what();
    std::vector<std::string> msgv; 
    split_str_by_newline(msgstr, msgv);
    return Rcpp::List::create(Rcpp::Named("status") = EXCEPTION_RC,
                              Rcpp::Named("msg") = Rcpp::List(msgv.begin(), msgv.end()));
  }

  Rcpp::List lst =
    Rcpp::List::create(Rcpp::Named("status") = SUCCESS_RC,
                       Rcpp::Named("model_cppname") = mname_,
                       Rcpp::Named("cppcode") = out.str());
  SEXP __sexp_result;
  PROTECT(__sexp_result = lst);
  UNPROTECT(1);
  return __sexp_result;
  END_RCPP
}

RcppExport SEXP stanfuncs(SEXP model_stancode, SEXP model_name, 
                          SEXP allow_undefined);

SEXP stanfuncs(SEXP model_stancode, 
               SEXP model_name, 
               SEXP allow_undefined) {
  BEGIN_RCPP;
  static const int SUCCESS_RC = 0;
  static const int EXCEPTION_RC = -1;
  static const int PARSE_FAIL_RC = -2;
  
  std::string mcode_ = Rcpp::as<std::string>(model_stancode);
  std::string mname_ = Rcpp::as<std::string>(model_name);
  std::vector<std::string> ns_(0);
  ns_.push_back(mname_);
  
  std::stringstream out;
  std::istringstream in(mcode_);
  try {
    bool valid_model
    = stan::lang::compile_functions(&rstan::io::rcerr,in,out,ns_,
                                    Rcpp::as<bool>(allow_undefined));
    if (!valid_model) {
      return Rcpp::List::create(Rcpp::Named("status") = PARSE_FAIL_RC);
      
    }
  } catch(const std::exception& e) {
    // REprintf("\nERROR PARSING\n %s\n", e.what());
    std::string msgstr = e.what();
    std::vector<std::string> msgv; 
    split_str_by_newline(msgstr, msgv);
    return Rcpp::List::create(Rcpp::Named("status") = EXCEPTION_RC,
                              Rcpp::Named("msg") = Rcpp::List(msgv.begin(), msgv.end()));
  }
  
  Rcpp::List lst =
    Rcpp::List::create(Rcpp::Named("status") = SUCCESS_RC,
                       Rcpp::Named("model_cppname") = mname_,
                       Rcpp::Named("cppcode") = out.str());
  SEXP __sexp_result;
  PROTECT(__sexp_result = lst);
  UNPROTECT(1);
  return __sexp_result;
  END_RCPP
}
