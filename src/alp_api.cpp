#include "alp_api.h"

#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

#include "sls_alignment_evaluer.hpp"

namespace {

void set_error(char *errbuf, int errbuf_len, const std::string &msg) {
    if (errbuf == nullptr || errbuf_len <= 0) return;
    std::snprintf(errbuf, static_cast<size_t>(errbuf_len), "%s", msg.c_str());
}

void copy_out(const Sls::ALP_set_of_parameters &p, gumbel_params *out) {
    out->lambda  = p.lambda;  out->lambda_error  = p.lambda_error;
    out->K       = p.K;       out->K_error       = p.K_error;
    out->a_I     = p.a_I;     out->a_I_error     = p.a_I_error;
    out->b_I     = p.b_I;     out->b_I_error     = p.b_I_error;
    out->a_J     = p.a_J;     out->a_J_error     = p.a_J_error;
    out->b_J     = p.b_J;     out->b_J_error     = p.b_J_error;
    out->alpha_I = p.alpha_I; out->alpha_I_error = p.alpha_I_error;
    out->beta_I  = p.beta_I;  out->beta_I_error  = p.beta_I_error;
    out->alpha_J = p.alpha_J; out->alpha_J_error = p.alpha_J_error;
    out->beta_J  = p.beta_J;  out->beta_J_error  = p.beta_J_error;
    out->sigma   = p.sigma;   out->sigma_error   = p.sigma_error;
    out->tau     = p.tau;     out->tau_error     = p.tau_error;
}

}  // namespace

extern "C" {

int gumbel_compute(const long   *mat,
                   int           n,
                   const double *freqs1,
                   const double *freqs2,
                   int           gapped,
                   long          gap_open1,
                   long          gap_extend1,
                   long          gap_open2,
                   long          gap_extend2,
                   int           insertions_after_deletions,
                   double        eps_lambda,
                   double        eps_K,
                   double        max_seconds,
                   double        max_mem_mb,
                   long          seed,
                   long          n_samples,
                   long          n_prelim_samples,
                   long         *samples_used,
                   gumbel_params *out,
                   char         *errbuf,
                   int           errbuf_len) {
    if (mat == nullptr || freqs1 == nullptr || freqs2 == nullptr || out == nullptr || n <= 0) {
        set_error(errbuf, errbuf_len, "invalid arguments");
        return GUMBEL_ERR_UNKNOWN;
    }
    std::memset(out, 0, sizeof(*out));
    if (samples_used) samples_used[0] = samples_used[1] = 0;

    /* ALP wants long** rows. */
    std::vector<const long *> rows(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) rows[i] = mat + static_cast<size_t>(i) * n;

    try {
        Sls::AlignmentEvaluer evaluer;
        if (gapped) {
            /* Adaptive mode needs a positive budget (ALP's own default is 60 s). */
            double alp_max_time = max_seconds > 0 ? max_seconds : 60.0;
            if (n_samples > 0) {
                /* Fixed plan: max_time <= 0 makes ALP replay the given sample
                 * counts; the time limit only acts as a hard abort. */
                evaluer.set_gapped_computation_parameters_simplified(
                    max_seconds > 0 ? max_seconds : -1.0,  /* <= 0: no time limit */
                    n_samples, n_prelim_samples > 0 ? n_prelim_samples : 200);
                alp_max_time = -1.0;
            }
            evaluer.initGapped(n, rows.data(), freqs1, freqs2,
                               gap_open1, gap_extend1, gap_open2, gap_extend2,
                               insertions_after_deletions != 0,
                               eps_lambda, eps_K, alp_max_time, max_mem_mb, seed);
            if (samples_used) {
                const auto &cp = evaluer.gapped_computation_parameters();
                samples_used[0] = cp.d_total_realizations_number_with_ALP;
                samples_used[1] = cp.d_total_realizations_number_with_killing;
            }
        } else {
            /* ALP replaces <= 0 by 60 s; a huge value means no limit. */
            evaluer.initGapless(n, rows.data(), freqs1, freqs2,
                                max_seconds > 0 ? max_seconds : 1e300);
        }
        if (!evaluer.isGood()) {
            set_error(errbuf, errbuf_len, "ALP finished without producing a complete parameter set");
            return GUMBEL_ERR_NOT_GOOD;
        }
        copy_out(evaluer.parameters(), out);
        return GUMBEL_OK;
    } catch (const Sls::error &e) {
        set_error(errbuf, errbuf_len, e.st);
        return GUMBEL_ERR_ALP;
    } catch (const std::exception &e) {
        set_error(errbuf, errbuf_len, e.what());
        return GUMBEL_ERR_UNKNOWN;
    } catch (...) {
        set_error(errbuf, errbuf_len, "unknown error in ALP");
        return GUMBEL_ERR_UNKNOWN;
    }
}

int gumbel_evalue(const gumbel_params *params,
                  double  score,
                  double  len1,
                  double  len2,
                  double *evalue,
                  double *pvalue,
                  double *area,
                  char   *errbuf,
                  int     errbuf_len) {
    if (params == nullptr || evalue == nullptr) {
        set_error(errbuf, errbuf_len, "invalid arguments");
        return GUMBEL_ERR_UNKNOWN;
    }
    /* ALP's AlignmentEvaluerParameters uses suffix 1 for J and 2 for I. */
    Sls::AlignmentEvaluerParameters p;
    p.d_lambda = params->lambda;
    p.d_k      = params->K;
    p.d_a1     = params->a_J;     p.d_b1    = params->b_J;
    p.d_a2     = params->a_I;     p.d_b2    = params->b_I;
    p.d_alpha1 = params->alpha_J; p.d_beta1 = params->beta_J;
    p.d_alpha2 = params->alpha_I; p.d_beta2 = params->beta_I;
    p.d_sigma  = params->sigma;
    p.d_tau    = params->tau;

    try {
        Sls::AlignmentEvaluer evaluer;
        evaluer.initParameters(p);
        const double a = evaluer.area(score, len1, len2);
        const double e = a * evaluer.evaluePerArea(score);
        *evalue = e;
        if (pvalue) *pvalue = Sls::AlignmentEvaluer::pvalue(e);
        if (area) *area = a;
        return GUMBEL_OK;
    } catch (const Sls::error &e) {
        set_error(errbuf, errbuf_len, e.st);
        return GUMBEL_ERR_ALP;
    } catch (const std::exception &e) {
        set_error(errbuf, errbuf_len, e.what());
        return GUMBEL_ERR_UNKNOWN;
    } catch (...) {
        set_error(errbuf, errbuf_len, "unknown error in ALP");
        return GUMBEL_ERR_UNKNOWN;
    }
}

}  /* extern "C" */
