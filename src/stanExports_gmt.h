// Generated by rstantools.  Do not edit by hand.

/*
    titertools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    titertools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with titertools.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.26.1-4-gd72b68b7-dirty
#include <stan/model/model_header.hpp>
namespace model_gmt_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'gmt', line 53, column 2 to column 10)",
                                                      " (in 'gmt', line 54, column 2 to column 25)",
                                                      " (in 'gmt', line 58, column 2 to column 43)",
                                                      " (in 'gmt', line 59, column 2 to column 57)",
                                                      " (in 'gmt', line 62, column 4 to line 66, column 6)",
                                                      " (in 'gmt', line 61, column 17 to line 67, column 3)",
                                                      " (in 'gmt', line 61, column 2 to line 67, column 3)",
                                                      " (in 'gmt', line 44, column 2 to column 17)",
                                                      " (in 'gmt', line 45, column 9 to column 10)",
                                                      " (in 'gmt', line 45, column 2 to column 23)",
                                                      " (in 'gmt', line 46, column 9 to column 10)",
                                                      " (in 'gmt', line 46, column 2 to column 23)",
                                                      " (in 'gmt', line 47, column 2 to column 19)",
                                                      " (in 'gmt', line 48, column 2 to column 22)",
                                                      " (in 'gmt', line 49, column 2 to column 25)",
                                                      " (in 'gmt', line 50, column 2 to column 24)",
                                                      " (in 'gmt', line 6, column 2 to column 14)",
                                                      " (in 'gmt', line 33, column 4 to line 35, column 6)",
                                                      " (in 'gmt', line 31, column 54 to line 37, column 3)",
                                                      " (in 'gmt', line 31, column 9 to line 37, column 3)",
                                                      " (in 'gmt', line 27, column 4 to line 29, column 6)",
                                                      " (in 'gmt', line 25, column 54 to line 31, column 3)",
                                                      " (in 'gmt', line 25, column 9 to line 37, column 3)",
                                                      " (in 'gmt', line 20, column 4 to line 23, column 6)",
                                                      " (in 'gmt', line 18, column 55 to line 25, column 3)",
                                                      " (in 'gmt', line 18, column 9 to line 37, column 3)",
                                                      " (in 'gmt', line 14, column 4 to line 16, column 6)",
                                                      " (in 'gmt', line 12, column 37 to line 18, column 3)",
                                                      " (in 'gmt', line 12, column 9 to line 37, column 3)",
                                                      " (in 'gmt', line 10, column 4 to column 15)",
                                                      " (in 'gmt', line 8, column 46 to line 12, column 3)",
                                                      " (in 'gmt', line 8, column 2 to line 37, column 3)",
                                                      " (in 'gmt', line 39, column 2 to column 16)",
                                                      " (in 'gmt', line 4, column 89 to line 41, column 1)"};
template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
normal_int_censored_likelihood(const T0__& lower_lim, const T1__& upper_lim,
                               const T2__& mu, const T3__& sigma,
                               std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__>;
  const static bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  
  try {
    local_scalar_t__ result;
    result = DUMMY_VAR__;
    
    current_statement__ = 32;
    if ((primitive_value(is_inf(lower_lim)) && primitive_value(
        is_inf(upper_lim)))) {
      current_statement__ = 30;
      result = 0;
    } else {
      current_statement__ = 29;
      if (logical_eq(lower_lim, upper_lim)) {
        current_statement__ = 27;
        result = normal_lpdf<false>(lower_lim, mu, sigma);
      } else {
        current_statement__ = 26;
        if ((primitive_value(logical_negation(is_inf(lower_lim))) &&
            primitive_value(logical_negation(is_inf(upper_lim))))) {
          current_statement__ = 24;
          result = log_diff_exp(normal_lcdf(upper_lim, mu, sigma),
                     normal_lcdf(lower_lim, mu, sigma));
        } else {
          current_statement__ = 23;
          if ((primitive_value(logical_negation(is_inf(lower_lim))) &&
              primitive_value(is_inf(upper_lim)))) {
            current_statement__ = 21;
            result = normal_lccdf(lower_lim, mu, sigma);
          } else {
            current_statement__ = 20;
            if ((primitive_value(is_inf(lower_lim)) && primitive_value(
                logical_negation(is_inf(upper_lim))))) {
              current_statement__ = 18;
              result = normal_lcdf(upper_lim, mu, sigma);
            } 
          }
        }
      }
    }
    current_statement__ = 33;
    return result;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}
struct normal_int_censored_likelihood_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
operator()(const T0__& lower_lim, const T1__& upper_lim, const T2__& mu,
           const T3__& sigma, std::ostream* pstream__)  const 
{
return normal_int_censored_likelihood(lower_lim, upper_lim, mu, sigma,
         pstream__);
}
};
#include <stan_meta_header.hpp>
class model_gmt final : public model_base_crtp<model_gmt> {
private:
  int N;
  Eigen::Matrix<double, -1, 1> upper_lims;
  Eigen::Matrix<double, -1, 1> lower_lims;
  double mu_prior_mu;
  double mu_prior_sigma;
  double sigma_prior_alpha;
  double sigma_prior_beta;
 
public:
  ~model_gmt() { }
  
  inline std::string model_name() const final { return "model_gmt"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-4-gd72b68b7-dirty", "stancflags = "};
  }
  
  
  model_gmt(stan::io::var_context& context__, unsigned int random_seed__ = 0,
            std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_gmt_namespace::model_gmt";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 8;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      N = std::numeric_limits<int>::min();
      
      current_statement__ = 8;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 8;
      current_statement__ = 8;
      check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 9;
      validate_non_negative_index("upper_lims", "N", N);
      current_statement__ = 10;
      context__.validate_dims("data initialization","upper_lims","double",
          context__.to_vec(N));
      upper_lims = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(upper_lims, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> upper_lims_flat__;
        current_statement__ = 10;
        assign(upper_lims_flat__, nil_index_list(),
          context__.vals_r("upper_lims"),
          "assigning variable upper_lims_flat__");
        current_statement__ = 10;
        pos__ = 1;
        current_statement__ = 10;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 10;
          assign(upper_lims, cons_list(index_uni(sym1__), nil_index_list()),
            upper_lims_flat__[(pos__ - 1)], "assigning variable upper_lims");
          current_statement__ = 10;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 11;
      validate_non_negative_index("lower_lims", "N", N);
      current_statement__ = 12;
      context__.validate_dims("data initialization","lower_lims","double",
          context__.to_vec(N));
      lower_lims = Eigen::Matrix<double, -1, 1>(N);
      stan::math::fill(lower_lims, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> lower_lims_flat__;
        current_statement__ = 12;
        assign(lower_lims_flat__, nil_index_list(),
          context__.vals_r("lower_lims"),
          "assigning variable lower_lims_flat__");
        current_statement__ = 12;
        pos__ = 1;
        current_statement__ = 12;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 12;
          assign(lower_lims, cons_list(index_uni(sym1__), nil_index_list()),
            lower_lims_flat__[(pos__ - 1)], "assigning variable lower_lims");
          current_statement__ = 12;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 13;
      context__.validate_dims("data initialization","mu_prior_mu","double",
          context__.to_vec());
      mu_prior_mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 13;
      mu_prior_mu = context__.vals_r("mu_prior_mu")[(1 - 1)];
      current_statement__ = 14;
      context__.validate_dims("data initialization","mu_prior_sigma",
          "double",context__.to_vec());
      mu_prior_sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 14;
      mu_prior_sigma = context__.vals_r("mu_prior_sigma")[(1 - 1)];
      current_statement__ = 15;
      context__.validate_dims("data initialization","sigma_prior_alpha",
          "double",context__.to_vec());
      sigma_prior_alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 15;
      sigma_prior_alpha = context__.vals_r("sigma_prior_alpha")[(1 - 1)];
      current_statement__ = 16;
      context__.validate_dims("data initialization","sigma_prior_beta",
          "double",context__.to_vec());
      sigma_prior_beta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 16;
      sigma_prior_beta = context__.vals_r("sigma_prior_beta")[(1 - 1)];
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_gmt_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      local_scalar_t__ mu;
      mu = DUMMY_VAR__;
      
      current_statement__ = 1;
      mu = in__.scalar();
      local_scalar_t__ sigma;
      sigma = DUMMY_VAR__;
      
      current_statement__ = 2;
      sigma = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        sigma = stan::math::lb_constrain(sigma, 0.01, lp__);
      } else {
        current_statement__ = 2;
        sigma = stan::math::lb_constrain(sigma, 0.01);
      }
      {
        current_statement__ = 3;
        lp_accum__.add(normal_lpdf<propto__>(mu, mu_prior_mu, mu_prior_sigma));
        current_statement__ = 4;
        lp_accum__.add(
          inv_gamma_lpdf<propto__>(sigma, sigma_prior_alpha,
            sigma_prior_beta));
        current_statement__ = 7;
        for (int i = 1; i <= N; ++i) {
          current_statement__ = 5;
          lp_accum__.add(
            normal_int_censored_likelihood(lower_lims[(i - 1)],
              upper_lims[(i - 1)], mu, sigma, pstream__));}
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_gmt_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      double mu;
      mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      mu = in__.scalar();
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      sigma = in__.scalar();
      current_statement__ = 2;
      sigma = stan::math::lb_constrain(sigma, 0.01);
      vars__.emplace_back(mu);
      vars__.emplace_back(sigma);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      double mu;
      mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      mu = context__.vals_r("mu")[(1 - 1)];
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      sigma = context__.vals_r("sigma")[(1 - 1)];
      double sigma_free__;
      sigma_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      sigma_free__ = stan::math::lb_free(sigma, 0.01);
      vars__.emplace_back(mu);
      vars__.emplace_back(sigma_free__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("mu");
    names__.emplace_back("sigma");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "mu");
    param_names__.emplace_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "mu");
    param_names__.emplace_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_gmt_namespace::model_gmt;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_gmt_namespace::profiles__;
}
#endif
#endif
