"""Frozen copy of the five ``utils.py`` functions the equivalence check needs.

This is a **verbatim** extraction, not a reimplementation. It exists so that
``check_utils_vs_lntest.py`` keeps running after ``utils.py`` is deleted in
stage 3d, which is what makes the RECOMB-era estimator reproducible by a
referee who only ever sees the post-refactor tree.

Do not edit these functions. They are the estimator that produced every
published RECOMB number, defects included -- in particular
``trigamma(x) = 1/x``, which is an integral approximation to psi_1 and not
psi_1 itself. See decision-retire-utils-py and handoff-math-questions.

Source: utils.py at sha256
    18040fec581a69a93c830ce28679badb4d2fe0d2941dd316058c7fb5ea280a6a
Extracted by AST, from these line ranges:
    trigamma           utils.py:9-10
    get_DELN_lfcs      utils.py:134-224
    get_LN_lfcs        utils.py:227-326
    compute_p_vals     utils.py:376-392
    get_t_statistic    utils.py:394-413

``check_utils_vs_lntest.py`` re-verifies this file against ``utils.py`` on
every run for as long as ``utils.py`` exists, and says so loudly when it is
gone and the check can no longer be made.
"""

import numpy as np
import scipy.stats as stats
from scipy.stats.distributions import t

def trigamma(x):
    return 1 / x  # + 0.5 / (x ** 2)


def get_DELN_lfcs(Y_, X_, normalize=True, test='t', normalization='CP10K', return_standard_error=False, return_log_abs_statistic=False):
    # Y is (n_cells, n_genes)
    eps = 1e-9

    Y = Y_.astype(float).copy()
    n = Y.shape[0]
    Y[Y <= 0] = np.nan  # Replace all non-positive with NaN
    n_plus = n - np.sum(np.isnan(Y), 0)

    X = X_.astype(float).copy()
    n_prime = X.shape[0]
    X[X <= 0] = np.nan
    n_plus_prime = n_prime - np.sum(np.isnan(X), 0)

    if normalize and (normalization == 'CP10K'):
        X = 1e4 * X / np.nansum(X, 1, keepdims=True)
        Y = 1e4 * Y / np.nansum(Y, 1, keepdims=True)

    elif normalize and (normalization == 'median-of-ratios'):
        # the normalization scheme proposed in DESeq2
        denom_Y = np.exp(np.nanmean(np.log(Y), 0))
        denom_Y[np.isnan(denom_Y)] = 1  # Avoid division by NaN for unexpressed genes
        c_Y = np.nanmedian(Y / denom_Y, 1, keepdims=True)
        Y /= c_Y

        denom_X = np.exp(np.nanmean(np.log(X), 0))
        denom_X[np.isnan(denom_X)] = 1  # Avoid division by NaN for unexpressed genes
        c_X = np.nanmedian(X / denom_X, 1, keepdims=True)
        X /= c_X


    pos_mean_Y = np.nanmean(Y, axis=0)
    pos_mean_X = np.nanmean(X, axis=0)

    # \hat{a}
    a_hat_Y = n_plus + (eps ** (1 + n_plus))
    a_hat_X = n_plus_prime + (eps ** (1 + n_plus_prime))

    # compute \log2\hat{theta} for each gene
    log2_theta_hat_Y = np.log2(a_hat_Y / n)
    log2_theta_hat_X = np.log2(a_hat_X / n_prime)

    # compute sample mean of positive counts
    log2_m_Y = np.log2(pos_mean_Y)
    log2_m_X = np.log2(pos_mean_X)

    lfc = (log2_theta_hat_Y + log2_m_Y) - (log2_theta_hat_X + log2_m_X)

    # compute standard errors
    se_Y_1 = trigamma(a_hat_Y) - trigamma(n)
    se_Y_2 = np.log(1 + np.nanvar(Y, axis=0) / (n_plus * (2 ** log2_m_Y) ** 2))
    se_Y = np.sqrt(se_Y_1 + se_Y_2) / np.log(2)

    se_X_1 = trigamma(a_hat_X) - trigamma(n_prime)
    se_X_2 = np.log(1 + np.nanvar(X, axis=0) / (n_plus_prime * (2 ** log2_m_X) ** 2))
    se_X = np.sqrt(se_X_1 + se_X_2) / np.log(2)

    if return_log_abs_statistic:
        # Return test statistics, log(|statistic|), AND p-values
        if test == 't':
            statistic, log_abs_statistic = get_t_statistic(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X,
                                                            se_Y, se_X, n, n_prime, return_log_abs_statistic=True)
            # Also compute p-values from the test statistic
            nu1, nu2 = n - 1, n_prime - 1
            se_combined_sq = se_Y ** 2 + se_X ** 2
            df = se_combined_sq ** 2 / (se_Y ** 4 / nu1 + se_X ** 4 / nu2)
            t_dist = t(df)
            p_vals = 2 * t_dist.sf(np.abs(statistic))
        else:
            # z-test
            statistic, log_abs_statistic = compute_p_vals(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X, 
                                                          se_Y, se_X, return_log_abs_statistic=True)
            # Also compute p-values from the test statistic
            p_vals = 2 * (1 - stats.norm.cdf(np.abs(statistic)))
        
        if return_standard_error:
            return lfc, statistic, log_abs_statistic, p_vals, np.sqrt(se_X ** 2 + se_Y ** 2)
        return lfc, statistic, log_abs_statistic, p_vals
    else:
        # Return p-values (original behavior)
        if test == 't':
            statistic, p_vals = get_t_statistic(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X,
                                                se_Y, se_X, n, n_prime, return_log_abs_statistic=False)
        else:
            # z-test
            statistic, p_vals = compute_p_vals(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X, 
                                               se_Y, se_X, return_log_abs_statistic=False)
        
        if return_standard_error:
            return lfc, p_vals, np.sqrt(se_X ** 2 + se_Y ** 2)
        return lfc, p_vals


def get_LN_lfcs(Y_, X_, normalize=True, test='t', normalization='CP10K', return_standard_error=False, return_log_abs_statistic=False):
    # Y is (n_cells, n_genes)

    G = Y_.shape[1]
    Y = Y_.astype(float).copy()
    n = Y.shape[0]
    Y[Y <= 0] = np.nan  # Replace all non-positive with NaN
    n_plus = n - np.sum(np.isnan(Y), 0)
    all_zeros_Y = (n_plus == 0)

    X = X_.astype(float).copy()
    n_prime = X.shape[0]
    X[X <= 0] = np.nan
    n_plus_prime = n_prime - np.sum(np.isnan(X), 0)
    all_zeros_X = (n_plus_prime == 0)

    if normalize and (normalization == 'CP10K'):
        X = 1e4 * X / np.nansum(X, 1, keepdims=True)
        Y = 1e4 * Y / np.nansum(Y, 1, keepdims=True)

    elif normalize and (normalization == 'median-of-ratios'):
        # the normalization scheme proposed in DESeq2
        denom_Y = np.exp(np.nanmean(np.log(Y), 0))
        denom_Y[np.isnan(denom_Y)] = 1  # Avoid division by NaN for unexpressed genes
        c_Y = np.nanmedian(Y / denom_Y, 1, keepdims=True)
        Y /= c_Y

        denom_X = np.exp(np.nanmean(np.log(X), 0))
        denom_X[np.isnan(denom_X)] = 1  # Avoid division by NaN for unexpressed genes
        c_X = np.nanmedian(X / denom_X, 1, keepdims=True)
        X /= c_X


    pos_mean_Y = np.ones(G, dtype=np.float32)  # to avoid NaNs in LFC when all counts are zero
    pos_mean_Y[~all_zeros_Y] = np.nanmean(Y[:, ~all_zeros_Y], axis=0)
    pos_mean_X = np.ones(G, dtype=np.float32)
    pos_mean_X[~all_zeros_X] = np.nanmean(X[:, ~all_zeros_X], axis=0)

    # \hat{a}
    a_hat_Y = np.ones(G, dtype=np.float32)  # to avoid NaNs in LFC when all counts are zero
    a_hat_Y[~all_zeros_Y] = n_plus[~all_zeros_Y]

    a_hat_X = np.ones(G, dtype=np.float32)
    a_hat_X[~all_zeros_X] = n_plus_prime[~all_zeros_X]

    # compute \log2\hat{theta} for each gene
    log2_theta_hat_Y = np.log2(a_hat_Y / n)
    log2_theta_hat_X = np.log2(a_hat_X / n_prime)

    # compute sample mean of positive counts
    log2_m_Y = np.log2(pos_mean_Y)
    log2_m_X = np.log2(pos_mean_X)

    lfc = (log2_theta_hat_Y + log2_m_Y) - (log2_theta_hat_X + log2_m_X)

    # compute standard errors
    se_Y_1 = trigamma(a_hat_Y) - trigamma(n)
    se_Y_2 = np.ones(G, dtype=np.float32)  # to avoid NaNs in SE when all counts are zero
    se_Y_2[~all_zeros_Y] = np.log(1 + np.nanvar(Y[:, ~all_zeros_Y], axis=0) / (n_plus[~all_zeros_Y] * (2 ** log2_m_Y[~all_zeros_Y]) ** 2))
    se_Y = np.sqrt(se_Y_1 + se_Y_2) / np.log(2)

    se_X_1 = trigamma(a_hat_X) - trigamma(n_prime)
    se_X_2 = np.ones(G, dtype=np.float32)
    se_X_2[~all_zeros_X] = np.log(1 + np.nanvar(X[:, ~all_zeros_X], axis=0) / (n_plus_prime[~all_zeros_X] * (2 ** log2_m_X[~all_zeros_X]) ** 2))
    se_X = np.sqrt(se_X_1 + se_X_2) / np.log(2)

    if return_log_abs_statistic:
        # Return test statistics, log(|statistic|), AND p-values
        if test == 't':
            statistic, log_abs_statistic = get_t_statistic(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X,
                                                            se_Y, se_X, n, n_prime, return_log_abs_statistic=True)
            # Also compute p-values from the test statistic
            nu1, nu2 = n - 1, n_prime - 1
            se_combined_sq = se_Y ** 2 + se_X ** 2
            df = se_combined_sq ** 2 / (se_Y ** 4 / nu1 + se_X ** 4 / nu2)
            t_dist = t(df)
            p_vals = 2 * t_dist.sf(np.abs(statistic))
        else:
            # z-test
            statistic, log_abs_statistic = compute_p_vals(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X, 
                                                          se_Y, se_X, return_log_abs_statistic=True)
            # Also compute p-values from the test statistic
            p_vals = 2 * (1 - stats.norm.cdf(np.abs(statistic)))
        
        if return_standard_error:
            return lfc, statistic, log_abs_statistic, p_vals, np.sqrt(se_X ** 2 + se_Y ** 2)
        return lfc, statistic, log_abs_statistic, p_vals
    else:
        # Return p-values (original behavior)
        if test == 't':
            statistic, p_vals = get_t_statistic(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X,
                                                se_Y, se_X, n, n_prime, return_log_abs_statistic=False)
        else:
            # z-test
            statistic, p_vals = compute_p_vals(log2_theta_hat_Y + log2_m_Y, log2_theta_hat_X + log2_m_X, 
                                               se_Y, se_X, return_log_abs_statistic=False)
        
        if return_standard_error:
            return lfc, p_vals, np.sqrt(se_X ** 2 + se_Y ** 2)
        return lfc, p_vals


def compute_p_vals(mean1, mean2, se1, se2, return_log_abs_statistic=False):
    # Compute the test statistic
    d = mean1 - mean2
    se_combined_sq = se1 ** 2 + se2 ** 2
    z_stat = d / (se_combined_sq ** 0.5)

    if return_log_abs_statistic:
        # Compute log(|z_stat|) using log-space arithmetic for numerical stability
        # log(|z_stat|) = log(|d|) - 0.5 * log(se_combined_sq)
        log_abs_d = np.log(np.abs(d) + np.finfo(float).tiny)  # Add tiny to avoid log(0)
        log_se_combined_sq = np.log(se_combined_sq + np.finfo(float).tiny)
        log_abs_z_stat = log_abs_d - 0.5 * log_se_combined_sq
        return z_stat, log_abs_z_stat
    else:
        # Compute the p-value for the two-tailed test
        p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))
        return z_stat, p_value


def get_t_statistic(mean1, mean2, se1, se2, n1, n2, return_log_abs_statistic=False):
    # implements two-sided t-test
    nu1, nu2 = n1 - 1, n2 - 1
    se_combined_sq = se1 ** 2 + se2 ** 2
    df = se_combined_sq ** 2 / (se1 ** 4 / nu1 + se2 ** 4 / nu2)
    d = mean1 - mean2
    denom = (se_combined_sq ** 0.5)
    t_statistic = d / denom
    t_dist = t(df)
    
    if return_log_abs_statistic:
        # Compute log(|t_statistic|) using log-space arithmetic for numerical stability
        # log(|t_statistic|) = log(|d|) - 0.5 * log(se_combined_sq)
        log_abs_d = np.log(np.abs(d) + np.finfo(float).tiny)  # Add tiny to avoid log(0)
        log_se_combined_sq = np.log(se_combined_sq + np.finfo(float).tiny)
        log_abs_t_stat = log_abs_d - 0.5 * log_se_combined_sq
        return t_statistic, log_abs_t_stat
    else:
        p_value = 2 * t_dist.sf(np.abs(t_statistic))
        return t_statistic, p_value
