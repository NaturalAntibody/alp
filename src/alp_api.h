/*
 * alp_api.h - minimal C API over the NCBI ALP library
 * (Ascending Ladder Program, Sheetlin & Frith) for estimating the Gumbel
 * parameters of local alignment score distributions.
 *
 * All functions are thread-safe (no shared state); errors are reported via
 * a return code plus an optional caller-provided message buffer.
 */
#ifndef ALP_API_H
#define ALP_API_H

#ifdef __cplusplus
extern "C" {
#endif

/* The 12 ALP parameters plus their estimated errors.
 * Index I refers to sequence #1 (query), J to sequence #2 (subject). */
typedef struct {
    double lambda,  lambda_error;
    double K,       K_error;
    double a_I,     a_I_error;
    double b_I,     b_I_error;
    double a_J,     a_J_error;
    double b_J,     b_J_error;
    double alpha_I, alpha_I_error;
    double beta_I,  beta_I_error;
    double alpha_J, alpha_J_error;
    double beta_J,  beta_J_error;
    double sigma,   sigma_error;
    double tau,     tau_error;
} gumbel_params;

#define GUMBEL_OK            0
#define GUMBEL_ERR_ALP       1   /* ALP raised an error (message in errbuf)   */
#define GUMBEL_ERR_NOT_GOOD  2   /* ALP finished but parameters are undefined */
#define GUMBEL_ERR_UNKNOWN   3   /* any other C++ exception                   */

/**
 * Estimate Gumbel parameters.
 *
 * Gap costs follow the NCBI convention used by ALP: a gap of length k costs
 * gap_open + k * gap_extend.  Penalties are positive numbers.
 *
 * @param mat          Row-major n x n integer score matrix.
 * @param n            Alphabet size.
 * @param freqs1       Background letter frequencies of sequence #1 (length n).
 * @param freqs2       Background letter frequencies of sequence #2 (length n).
 * @param gapped       0 = gapless (analytic, gap args ignored), 1 = gapped (simulation).
 * @param gap_open1/gap_extend1   Gap penalties for sequence #1.
 * @param gap_open2/gap_extend2   Gap penalties for sequence #2.
 * @param insertions_after_deletions  Allow an insertion directly after a deletion.
 * @param eps_lambda   Target relative error of lambda (e.g. 0.01).
 * @param eps_K        Target relative error of K (e.g. 0.05).
 * @param max_seconds  Time limit.  Gapless and fixed-plan mode: a hard limit
 *                     (exceeding it is an error) and <= 0 means no limit.
 *                     Adaptive mode: the time budget ALP plans its sampling
 *                     around (ALP uses 60 s if <= 0).
 * @param max_mem_mb   Memory limit in megabytes.
 * @param seed         Random seed for the importance-sampling simulation.
 * @param n_samples    > 0: fixed-plan mode.  Use exactly this many realizations
 *                     for the main ALP and killing stages, and
 *                     n_prelim_samples for each preliminary stage.  No
 *                     decision depends on elapsed time, so the result is a pure
 *                     function of the inputs and the seed; max_seconds is only
 *                     a hard timeout (error, never a different result) and
 *                     eps_lambda / eps_K are ignored.
 *                     <= 0: adaptive mode (ALP default).  Sample sizes are
 *                     chosen from eps_* and from measured run time, so results
 *                     can vary with machine speed / load when max_seconds binds.
 * @param n_prelim_samples  Preliminary-stage realizations in fixed-plan mode.
 * @param samples_used Optional long[2]: receives the realizations actually used
 *                     for the ALP and killing stages (0 for gapless).
 * @param out          Receives the parameters on success.
 * @param errbuf       Optional buffer for an error message (may be NULL).
 * @param errbuf_len   Size of errbuf.
 * @return GUMBEL_OK on success, otherwise a GUMBEL_ERR_* code.
 */
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
                   int           errbuf_len);

/**
 * Compute an E-value / P-value with ALP's finite-size ("area") correction.
 *
 * @param params  Parameters (error fields are ignored).
 * @param score   Raw alignment score.
 * @param len1    Length of sequence #1 (query).
 * @param len2    Length of sequence #2 (subject / database size).
 * @param evalue  Receives the E-value.
 * @param pvalue  Receives the P-value (may be NULL).
 * @param area    Receives the corrected search-space area (may be NULL).
 */
int gumbel_evalue(const gumbel_params *params,
                  double  score,
                  double  len1,
                  double  len2,
                  double *evalue,
                  double *pvalue,
                  double *area,
                  char   *errbuf,
                  int     errbuf_len);

#ifdef __cplusplus
}
#endif

#endif /* ALP_API_H */
