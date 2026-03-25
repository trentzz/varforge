//! Chi-squared goodness-of-fit test for stochastic VAF sampling.
//!
//! This validates the core paper claim: per-read Bernoulli sampling produces
//! observed alt-read counts that follow a binomial distribution.
//!
//! For each target VAF in {0.05, 0.1, 0.25, 0.5}, we draw N = 500 independent
//! samples from Binomial(100, vaf) using `sample_alt_count`, bin the counts,
//! and test the binned distribution against the theoretical binomial PMF using
//! a chi-squared goodness-of-fit test.  We assert that the p-value exceeds 0.01.
//!
//! The test uses a fixed seed so it is deterministic and never flaky.

use rand::rngs::StdRng;
use rand::SeedableRng;
use varforge::variants::vaf::sample_alt_count;

// ---------------------------------------------------------------------------
// Incomplete gamma function for chi-squared p-values
// ---------------------------------------------------------------------------

/// Regularised lower incomplete gamma function P(a, x).
///
/// Uses the series expansion for x < a + 1, and the continued fraction
/// expansion otherwise.  Accuracy is sufficient for hypothesis testing
/// (absolute error < 1e-7 for the parameter ranges used here).
fn reg_lower_inc_gamma(a: f64, x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        gamma_series(a, x)
    } else {
        1.0 - gamma_continued_fraction(a, x)
    }
}

/// Series expansion for the regularised lower incomplete gamma.
fn gamma_series(a: f64, x: f64) -> f64 {
    let max_iter = 300;
    let epsilon = 1e-12;
    let log_gamma_a = ln_gamma(a);
    let mut term = 1.0 / a;
    let mut sum = term;
    for n in 1..=max_iter {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < epsilon * sum.abs() {
            break;
        }
    }
    (-x + a * x.ln() - log_gamma_a).exp() * sum
}

/// Continued fraction expansion for the regularised upper incomplete gamma.
///
/// Returns Q(a, x) = 1 - P(a, x).
fn gamma_continued_fraction(a: f64, x: f64) -> f64 {
    let max_iter = 300;
    let epsilon = 1e-12;
    let fpmin = f64::MIN_POSITIVE / epsilon;
    let log_gamma_a = ln_gamma(a);

    // Modified Lentz method.
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / fpmin;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..=max_iter {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = b + an / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < epsilon {
            break;
        }
    }
    (-x + a * x.ln() - log_gamma_a).exp() * h
}

/// Natural logarithm of the gamma function, using Lanczos approximation.
fn ln_gamma(x: f64) -> f64 {
    // Lanczos coefficients (g=7, n=9).
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_571e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let mut x = x;
    if x < 0.5 {
        // Reflection formula: ln(Gamma(x)) = ln(pi) - ln(sin(pi*x)) - ln(Gamma(1-x))
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x);
    }
    x -= 1.0;
    let t = x + G + 0.5;
    let mut ser = C[0];
    for (i, &c) in C[1..].iter().enumerate() {
        ser += c / (x + (i + 1) as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + ser.ln()
}

/// Chi-squared survival function: Pr(X > stat) where X ~ chi^2(df).
///
/// Returns the p-value for the chi-squared goodness-of-fit test.
fn chi2_pvalue(stat: f64, df: f64) -> f64 {
    1.0 - reg_lower_inc_gamma(df / 2.0, stat / 2.0)
}

// ---------------------------------------------------------------------------
// Binomial PMF (exact, using log-factorial)
// ---------------------------------------------------------------------------

/// Binomial PMF: Pr(X = k | n, p).
fn binom_pmf(n: u32, k: u32, p: f64) -> f64 {
    if k > n {
        return 0.0;
    }
    if p == 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    if p >= 1.0 {
        return if k == n { 1.0 } else { 0.0 };
    }
    let log_coef = log_binom_coef(n, k);
    (log_coef + (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln()).exp()
}

/// Log binomial coefficient ln(C(n, k)) using Stirling-like log-factorial.
fn log_binom_coef(n: u32, k: u32) -> f64 {
    ln_gamma((n + 1) as f64) - ln_gamma((k + 1) as f64) - ln_gamma((n - k + 1) as f64)
}

// ---------------------------------------------------------------------------
// Chi-squared goodness-of-fit test
// ---------------------------------------------------------------------------

/// Run a chi-squared goodness-of-fit test of observed alt-read counts against
/// Binomial(depth, vaf).
///
/// Bins are constructed by merging adjacent PMF bins until the expected count
/// in each bin is at least `min_expected`.  This ensures the chi-squared
/// approximation is valid.
///
/// Returns (chi-squared statistic, p-value, degrees of freedom).
fn chi2_goodness_of_fit(
    counts: &[u32],
    depth: u32,
    vaf: f64,
    min_expected: f64,
) -> (f64, f64, u32) {
    let n = counts.len() as f64;

    // Compute the raw PMF for every possible outcome in [0, depth].
    let pmf: Vec<f64> = (0..=depth).map(|k| binom_pmf(depth, k, vaf)).collect();

    // Merge bins to ensure each has expected count >= min_expected.
    // We scan from 0 upward, accumulating until the bin is large enough,
    // then scan from depth downward to merge the right tail.
    // Finally combine: bins are (expected_count, observed_count).
    let mut bins: Vec<(f64, u32)> = Vec::new();

    // Build a frequency table for observed counts.
    let mut freq = vec![0u32; (depth + 1) as usize];
    for &c in counts {
        if (c as usize) < freq.len() {
            freq[c as usize] += 1;
        }
    }

    // Merge bins left-to-right until expected >= min_expected.
    let mut exp_acc = 0.0;
    let mut obs_acc = 0u32;
    for k in 0..=depth {
        exp_acc += pmf[k as usize] * n;
        obs_acc += freq[k as usize];
        if exp_acc >= min_expected {
            bins.push((exp_acc, obs_acc));
            exp_acc = 0.0;
            obs_acc = 0;
        }
    }
    // Any remaining mass goes into the last bin.
    if exp_acc > 0.0 {
        if let Some(last) = bins.last_mut() {
            last.0 += exp_acc;
            last.1 += obs_acc;
        }
    }

    // Compute chi-squared statistic.
    let stat: f64 = bins
        .iter()
        .map(|&(exp, obs)| {
            let diff = obs as f64 - exp;
            diff * diff / exp
        })
        .sum();

    let df = (bins.len().saturating_sub(1)) as u32;
    let pval = if df == 0 {
        1.0
    } else {
        chi2_pvalue(stat, df as f64)
    };

    (stat, pval, df)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Validate that `sample_alt_count` produces alt-read counts consistent with
/// the Binomial(depth, vaf) distribution at each VAF level.
///
/// For each VAF in {0.05, 0.1, 0.25, 0.5}:
/// - Draw N = 500 samples at depth = 100.
/// - Compute chi-squared goodness-of-fit against Binomial(100, vaf).
/// - Assert p-value > 0.01.
///
/// A fixed seed makes the test deterministic.
#[test]
fn test_vaf_follows_binomial_distribution() {
    const DEPTH: u32 = 100;
    const N_SAMPLES: usize = 500;
    const SEED: u64 = 12345;
    const ALPHA: f64 = 0.01;
    const MIN_EXPECTED: f64 = 5.0;

    let vafs = [0.05_f64, 0.1, 0.25, 0.5];

    for &vaf in &vafs {
        let mut rng = StdRng::seed_from_u64(SEED);
        let counts: Vec<u32> = (0..N_SAMPLES)
            .map(|_| sample_alt_count(DEPTH, vaf, &mut rng))
            .collect();

        let (stat, pval, df) = chi2_goodness_of_fit(&counts, DEPTH, vaf, MIN_EXPECTED);

        assert!(
            pval > ALPHA,
            "VAF={vaf}: chi-squared goodness-of-fit failed (stat={stat:.4}, df={df}, p={pval:.4}). \
             Alt-read counts deviate significantly from Binomial({DEPTH}, {vaf}). \
             This would invalidate the stochastic VAF paper claim."
        );
    }
}

/// Sanity check: a deliberately biased sampler (always returning floor(vaf * depth))
/// should fail the chi-squared test, confirming the test has power.
///
/// This guards against the test being trivially non-discriminating.
#[test]
fn test_chi2_rejects_deterministic_spike_in() {
    const DEPTH: u32 = 100;
    const N_SAMPLES: usize = 500;
    const ALPHA: f64 = 0.01;
    const MIN_EXPECTED: f64 = 5.0;

    // Simulate a deterministic spike-in: every sample returns exactly floor(vaf * depth).
    let vaf = 0.1_f64;
    let fixed_count = (vaf * DEPTH as f64).floor() as u32; // always 10
    let counts: Vec<u32> = vec![fixed_count; N_SAMPLES];

    let (_stat, pval, _df) = chi2_goodness_of_fit(&counts, DEPTH, vaf, MIN_EXPECTED);

    assert!(
        pval <= ALPHA,
        "The chi-squared test should reject a deterministic spike-in (p={pval:.6}), \
         but it did not. The test lacks statistical power."
    );
}
