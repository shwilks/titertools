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
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_gmt_me_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_gmt_me");
    reader.add_event(2, 2, "include", "/functions/normal_int_censored_likelihood.stan");
    reader.add_event(2, 0, "start", "/functions/normal_int_censored_likelihood.stan");
    reader.add_event(41, 39, "end", "/functions/normal_int_censored_likelihood.stan");
    reader.add_event(41, 3, "restart", "model_gmt_me");
    reader.add_event(92, 52, "end", "model_gmt_me");
    return reader;
}
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
normal_int_censored_likelihood(const T0__& lower_lim,
                                   const T1__& upper_lim,
                                   const T2__& mu,
                                   const T3__& sigma, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 6;
        local_scalar_t__ result(DUMMY_VAR__);
        (void) result;  // dummy to suppress unused var warning
        stan::math::initialize(result, DUMMY_VAR__);
        stan::math::fill(result, DUMMY_VAR__);
        current_statement_begin__ = 8;
        if (as_bool((primitive_value(is_inf(lower_lim)) && primitive_value(is_inf(upper_lim))))) {
            current_statement_begin__ = 10;
            stan::math::assign(result, 0);
        } else if (as_bool(logical_eq(lower_lim, upper_lim))) {
            current_statement_begin__ = 14;
            stan::math::assign(result, normal_log(lower_lim, mu, sigma));
        } else if (as_bool((primitive_value(logical_negation(is_inf(lower_lim))) && primitive_value(logical_negation(is_inf(upper_lim)))))) {
            current_statement_begin__ = 20;
            stan::math::assign(result, log_diff_exp(normal_cdf_log(upper_lim, mu, sigma), normal_cdf_log(lower_lim, mu, sigma)));
        } else if (as_bool((primitive_value(logical_negation(is_inf(lower_lim))) && primitive_value(is_inf(upper_lim))))) {
            current_statement_begin__ = 27;
            stan::math::assign(result, normal_ccdf_log(lower_lim, mu, sigma));
        } else if (as_bool((primitive_value(is_inf(lower_lim)) && primitive_value(logical_negation(is_inf(upper_lim)))))) {
            current_statement_begin__ = 33;
            stan::math::assign(result, normal_cdf_log(upper_lim, mu, sigma));
        }
        current_statement_begin__ = 39;
        return stan::math::promote_scalar<fun_return_scalar_t__>(result);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct normal_int_censored_likelihood_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const T0__& lower_lim,
                                   const T1__& upper_lim,
                                   const T2__& mu,
                                   const T3__& sigma, std::ostream* pstream__) const {
        return normal_int_censored_likelihood(lower_lim, upper_lim, mu, sigma, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_gmt_me
  : public stan::model::model_base_crtp<model_gmt_me> {
private:
        int N;
        int N_ags;
        int N_srs;
        vector_d upper_lims;
        vector_d lower_lims;
        std::vector<int> ag;
        std::vector<int> sr;
        double mu_prior_mu;
        double mu_prior_sigma;
        double sigma_prior_alpha;
        double sigma_prior_beta;
        double sigma_prior_ag_effect;
        double sigma_prior_sr_effect;
public:
    model_gmt_me(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_gmt_me(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_gmt_me_namespace::model_gmt_me";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 45;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 1);
            current_statement_begin__ = 46;
            context__.validate_dims("data initialization", "N_ags", "int", context__.to_vec());
            N_ags = int(0);
            vals_i__ = context__.vals_i("N_ags");
            pos__ = 0;
            N_ags = vals_i__[pos__++];
            check_greater_or_equal(function__, "N_ags", N_ags, 1);
            current_statement_begin__ = 47;
            context__.validate_dims("data initialization", "N_srs", "int", context__.to_vec());
            N_srs = int(0);
            vals_i__ = context__.vals_i("N_srs");
            pos__ = 0;
            N_srs = vals_i__[pos__++];
            check_greater_or_equal(function__, "N_srs", N_srs, 1);
            current_statement_begin__ = 48;
            validate_non_negative_index("upper_lims", "N", N);
            context__.validate_dims("data initialization", "upper_lims", "vector_d", context__.to_vec(N));
            upper_lims = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("upper_lims");
            pos__ = 0;
            size_t upper_lims_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < upper_lims_j_1_max__; ++j_1__) {
                upper_lims(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 49;
            validate_non_negative_index("lower_lims", "N", N);
            context__.validate_dims("data initialization", "lower_lims", "vector_d", context__.to_vec(N));
            lower_lims = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("lower_lims");
            pos__ = 0;
            size_t lower_lims_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < lower_lims_j_1_max__; ++j_1__) {
                lower_lims(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 50;
            validate_non_negative_index("ag", "N", N);
            context__.validate_dims("data initialization", "ag", "int", context__.to_vec(N));
            ag = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("ag");
            pos__ = 0;
            size_t ag_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < ag_k_0_max__; ++k_0__) {
                ag[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 51;
            validate_non_negative_index("sr", "N", N);
            context__.validate_dims("data initialization", "sr", "int", context__.to_vec(N));
            sr = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("sr");
            pos__ = 0;
            size_t sr_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < sr_k_0_max__; ++k_0__) {
                sr[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 52;
            context__.validate_dims("data initialization", "mu_prior_mu", "double", context__.to_vec());
            mu_prior_mu = double(0);
            vals_r__ = context__.vals_r("mu_prior_mu");
            pos__ = 0;
            mu_prior_mu = vals_r__[pos__++];
            current_statement_begin__ = 53;
            context__.validate_dims("data initialization", "mu_prior_sigma", "double", context__.to_vec());
            mu_prior_sigma = double(0);
            vals_r__ = context__.vals_r("mu_prior_sigma");
            pos__ = 0;
            mu_prior_sigma = vals_r__[pos__++];
            current_statement_begin__ = 54;
            context__.validate_dims("data initialization", "sigma_prior_alpha", "double", context__.to_vec());
            sigma_prior_alpha = double(0);
            vals_r__ = context__.vals_r("sigma_prior_alpha");
            pos__ = 0;
            sigma_prior_alpha = vals_r__[pos__++];
            current_statement_begin__ = 55;
            context__.validate_dims("data initialization", "sigma_prior_beta", "double", context__.to_vec());
            sigma_prior_beta = double(0);
            vals_r__ = context__.vals_r("sigma_prior_beta");
            pos__ = 0;
            sigma_prior_beta = vals_r__[pos__++];
            current_statement_begin__ = 56;
            context__.validate_dims("data initialization", "sigma_prior_ag_effect", "double", context__.to_vec());
            sigma_prior_ag_effect = double(0);
            vals_r__ = context__.vals_r("sigma_prior_ag_effect");
            pos__ = 0;
            sigma_prior_ag_effect = vals_r__[pos__++];
            current_statement_begin__ = 57;
            context__.validate_dims("data initialization", "sigma_prior_sr_effect", "double", context__.to_vec());
            sigma_prior_sr_effect = double(0);
            vals_r__ = context__.vals_r("sigma_prior_sr_effect");
            pos__ = 0;
            sigma_prior_sr_effect = vals_r__[pos__++];
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 61;
            num_params_r__ += 1;
            current_statement_begin__ = 62;
            num_params_r__ += 1;
            current_statement_begin__ = 63;
            validate_non_negative_index("ag_effects", "N_ags", N_ags);
            num_params_r__ += N_ags;
            current_statement_begin__ = 64;
            validate_non_negative_index("sr_effects", "N_srs", N_srs);
            num_params_r__ += N_srs;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_gmt_me() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 61;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu", "double", context__.to_vec());
        double mu(0);
        mu = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 62;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec());
        double sigma(0);
        sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0.01, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 63;
        if (!(context__.contains_r("ag_effects")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable ag_effects missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("ag_effects");
        pos__ = 0U;
        validate_non_negative_index("ag_effects", "N_ags", N_ags);
        context__.validate_dims("parameter initialization", "ag_effects", "vector_d", context__.to_vec(N_ags));
        Eigen::Matrix<double, Eigen::Dynamic, 1> ag_effects(N_ags);
        size_t ag_effects_j_1_max__ = N_ags;
        for (size_t j_1__ = 0; j_1__ < ag_effects_j_1_max__; ++j_1__) {
            ag_effects(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(ag_effects);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable ag_effects: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 64;
        if (!(context__.contains_r("sr_effects")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sr_effects missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sr_effects");
        pos__ = 0U;
        validate_non_negative_index("sr_effects", "N_srs", N_srs);
        context__.validate_dims("parameter initialization", "sr_effects", "vector_d", context__.to_vec(N_srs));
        Eigen::Matrix<double, Eigen::Dynamic, 1> sr_effects(N_srs);
        size_t sr_effects_j_1_max__ = N_srs;
        for (size_t j_1__ = 0; j_1__ < sr_effects_j_1_max__; ++j_1__) {
            sr_effects(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(sr_effects);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sr_effects: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 61;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.scalar_constrain(lp__);
            else
                mu = in__.scalar_constrain();
            current_statement_begin__ = 62;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.scalar_lb_constrain(0.01, lp__);
            else
                sigma = in__.scalar_lb_constrain(0.01);
            current_statement_begin__ = 63;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> ag_effects;
            (void) ag_effects;  // dummy to suppress unused var warning
            if (jacobian__)
                ag_effects = in__.vector_constrain(N_ags, lp__);
            else
                ag_effects = in__.vector_constrain(N_ags);
            current_statement_begin__ = 64;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> sr_effects;
            (void) sr_effects;  // dummy to suppress unused var warning
            if (jacobian__)
                sr_effects = in__.vector_constrain(N_srs, lp__);
            else
                sr_effects = in__.vector_constrain(N_srs);
            // model body
            current_statement_begin__ = 70;
            lp_accum__.add(normal_log<propto__>(mu, mu_prior_mu, mu_prior_sigma));
            current_statement_begin__ = 71;
            lp_accum__.add(inv_gamma_log<propto__>(sigma, sigma_prior_alpha, sigma_prior_beta));
            current_statement_begin__ = 72;
            lp_accum__.add(normal_log<propto__>(ag_effects, 0, sigma_prior_ag_effect));
            current_statement_begin__ = 73;
            lp_accum__.add(normal_log<propto__>(sr_effects, 0, sigma_prior_sr_effect));
            current_statement_begin__ = 76;
            for (int i = 1; i <= N; ++i) {
                {
                current_statement_begin__ = 79;
                local_scalar_t__ logtiter(DUMMY_VAR__);
                (void) logtiter;  // dummy to suppress unused var warning
                stan::math::initialize(logtiter, DUMMY_VAR__);
                stan::math::fill(logtiter, DUMMY_VAR__);
                stan::math::assign(logtiter,((mu + get_base1(ag_effects, get_base1(ag, i, "ag", 1), "ag_effects", 1)) + get_base1(sr_effects, get_base1(sr, i, "sr", 1), "sr_effects", 1)));
                current_statement_begin__ = 82;
                lp_accum__.add(normal_int_censored_likelihood(get_base1(lower_lims, i, "lower_lims", 1), get_base1(upper_lims, i, "upper_lims", 1), logtiter, sigma, pstream__));
                }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("mu");
        names__.push_back("sigma");
        names__.push_back("ag_effects");
        names__.push_back("sr_effects");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_ags);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_srs);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_gmt_me_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double mu = in__.scalar_constrain();
        vars__.push_back(mu);
        double sigma = in__.scalar_lb_constrain(0.01);
        vars__.push_back(sigma);
        Eigen::Matrix<double, Eigen::Dynamic, 1> ag_effects = in__.vector_constrain(N_ags);
        size_t ag_effects_j_1_max__ = N_ags;
        for (size_t j_1__ = 0; j_1__ < ag_effects_j_1_max__; ++j_1__) {
            vars__.push_back(ag_effects(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> sr_effects = in__.vector_constrain(N_srs);
        size_t sr_effects_j_1_max__ = N_srs;
        for (size_t j_1__ = 0; j_1__ < sr_effects_j_1_max__; ++j_1__) {
            vars__.push_back(sr_effects(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_gmt_me";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        size_t ag_effects_j_1_max__ = N_ags;
        for (size_t j_1__ = 0; j_1__ < ag_effects_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "ag_effects" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t sr_effects_j_1_max__ = N_srs;
        for (size_t j_1__ = 0; j_1__ < sr_effects_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sr_effects" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        size_t ag_effects_j_1_max__ = N_ags;
        for (size_t j_1__ = 0; j_1__ < ag_effects_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "ag_effects" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t sr_effects_j_1_max__ = N_srs;
        for (size_t j_1__ = 0; j_1__ < sr_effects_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sr_effects" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_gmt_me_namespace::model_gmt_me stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
