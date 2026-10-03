//! Stochastic methylation kinetics. Biological states and observation likelihoods stay separate.
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::f64::consts::PI;

pub fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}
pub fn variance(x: &[f64]) -> f64 {
    let m = mean(x);
    x.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (x.len().saturating_sub(1).max(1)) as f64
}
pub fn normal_logpdf(x: f64, m: f64, v: f64) -> f64 {
    let v = v.max(1e-12);
    -0.5 * ((2.0 * PI * v).ln() + (x - m).powi(2) / v)
}
pub fn logsumexp(x: &[f64]) -> f64 {
    let m = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !m.is_finite() {
        return m;
    }
    m + x.iter().map(|a| (a - m).exp()).sum::<f64>().ln()
}

/// Deterministic bounded coordinate search with several starts supplied by the caller.
/// Returns the objective and parameters; this is maximum likelihood, not posterior inference.
pub fn optimize(
    f: impl Fn(&[f64]) -> f64,
    starts: &[Vec<f64>],
    bounds: &[(f64, f64)],
) -> (f64, Vec<f64>) {
    let mut best = (f64::INFINITY, vec![]);
    for start in starts {
        let mut x = start.clone();
        let mut score = f(&x);
        let mut steps: Vec<f64> = bounds.iter().map(|(l, h)| (h - l) * 0.15).collect();
        for _ in 0..100 {
            let mut improved = false;
            for j in 0..x.len() {
                for sign in [-1.0, 1.0] {
                    let mut y = x.clone();
                    y[j] = (y[j] + sign * steps[j]).clamp(bounds[j].0, bounds[j].1);
                    let s = f(&y);
                    if s < score {
                        score = s;
                        x = y;
                        improved = true;
                    }
                }
            }
            if !improved {
                for s in &mut steps {
                    *s *= 0.5;
                }
                if steps
                    .iter()
                    .zip(bounds)
                    .all(|(s, (l, h))| *s < (h - l) * 1e-5)
                {
                    break;
                }
            }
        }
        if score < best.0 {
            best = (score, x);
        }
    }
    best
}

/// Exact marginal probability for a two-state continuous-time chain.
pub fn probability(initial: f64, gain: f64, loss: f64, t: f64) -> f64 {
    assert!((0.0..=1.0).contains(&initial) && gain >= 0.0 && loss >= 0.0 && t >= 0.0);
    let rate = gain + loss;
    if rate == 0.0 {
        return initial;
    }
    let equilibrium = gain / rate;
    equilibrium + (initial - equilibrium) * (-rate * t).exp()
}

/// Microscopically derived drift and instantaneous aggregate variance for N independent alleles.
pub fn aggregate_moments(x: f64, gain: f64, loss: f64, n: usize) -> (f64, f64) {
    assert!((0.0..=1.0).contains(&x) && n > 0 && gain >= 0.0 && loss >= 0.0);
    (
        gain * (1.0 - x) - loss * x,
        (gain * (1.0 - x) + loss * x) / n as f64,
    )
}

/// Rounded methylated count is never inferred from array beta values.
#[derive(Clone, Debug)]
pub enum Observation {
    Binary(bool),
    Counts {
        methylated: u32,
        total: u32,
        concentration: Option<f64>,
    },
    Array {
        beta: f64,
        sd: f64,
    },
    Missing,
}
fn log_gamma(x: f64) -> f64 {
    // Lanczos approximation; positive arguments only.
    let c = [
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        return PI.ln() - (PI * x).sin().ln() - log_gamma(1.0 - x);
    }
    let z = x - 1.0;
    let mut a = 0.9999999999998099;
    for (i, c) in c.iter().enumerate() {
        a += c / (z + i as f64 + 1.0);
    }
    let t = z + 7.5;
    0.5 * (2.0 * PI).ln() + (z + 0.5) * t.ln() - t + a.ln()
}
fn log_beta(a: f64, b: f64) -> f64 {
    log_gamma(a) + log_gamma(b) - log_gamma(a + b)
}
pub fn log_likelihood(obs: &Observation, p: f64) -> Result<f64, String> {
    if !p.is_finite() || !(0.0..=1.0).contains(&p) {
        return Err("invalid probability".into());
    }
    let ll = match *obs {
        Observation::Missing => 0.0,
        Observation::Binary(b) => {
            if b {
                p.ln()
            } else {
                (-p).ln_1p()
            }
        }
        Observation::Array { beta, sd } => {
            if !(0.0..=1.0).contains(&beta) || !sd.is_finite() || sd <= 0.0 {
                return Err("invalid array observation".into());
            }
            normal_logpdf(beta, p, sd * sd)
        }
        Observation::Counts {
            methylated: k,
            total: n,
            concentration: c,
        } => {
            if k > n {
                return Err("methylated count exceeds coverage".into());
            }
            if n == 0 {
                return Ok(0.0);
            }
            if p == 0.0 {
                return Ok(if k == 0 { 0.0 } else { f64::NEG_INFINITY });
            }
            if p == 1.0 {
                return Ok(if k == n { 0.0 } else { f64::NEG_INFINITY });
            }
            let choose = log_gamma(n as f64 + 1.0)
                - log_gamma(k as f64 + 1.0)
                - log_gamma((n - k) as f64 + 1.0);
            if let Some(c) = c {
                if !c.is_finite() || c <= 0.0 {
                    return Err("invalid beta-binomial concentration".into());
                }
                let a = p * c;
                let b = (1.0 - p) * c;
                choose + log_beta(k as f64 + a, (n - k) as f64 + b) - log_beta(a, b)
            } else {
                choose + k as f64 * p.ln() + (n - k) as f64 * (-p).ln_1p()
            }
        }
    };
    Ok(ll)
}

/// Published slow/fast Polycomb baseline, with the switch time integrated out.
/// Continuous switch-time quadrature and Gaussian approximation to weekly Bernoulli gains.
/// Values are beta fractions; increment events add one percentage point.
pub fn polycomb_logpdf(y: f64, t: f64, initial: f64, initial_var: f64, p: &[f64]) -> f64 {
    let slow = p[0];
    let fast = if p.len() > 2 { p[1] } else { slow };
    let rho = if p.len() > 2 { p[2] } else { 0.0 };
    let extra = *p.last().unwrap();
    let mut terms = vec![
        -rho * t
            + normal_logpdf(
                y,
                initial + 0.01 * slow * t,
                initial_var + extra * extra + 0.0001 * slow * (1.0 - slow) * t,
            ),
    ];
    if rho > 0.0 && t > 0.0 {
        let n = 48;
        let dt = t / n as f64;
        for j in 0..n {
            let s = (j as f64 + 0.5) * dt;
            let mass = (-rho * (j as f64 * dt)).exp() - (-rho * ((j + 1) as f64 * dt)).exp();
            let m = initial + 0.01 * (slow * s + fast * (t - s));
            let v = initial_var
                + extra * extra
                + 0.0001 * (slow * (1.0 - slow) * s + fast * (1.0 - fast) * (t - s));
            terms.push(mass.ln() + normal_logpdf(y, m, v));
        }
    }
    logsumexp(&terms)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn chain_semigroup() {
        let p = probability(0.2, 0.1, 0.3, 7.0);
        assert!((probability(p, 0.1, 0.3, 4.0) - probability(0.2, 0.1, 0.3, 11.0)).abs() < 1e-12);
        assert_eq!(probability(0.3, 0.0, 0.0, 100.0), 0.3);
    }
    #[test]
    fn measurement_normalizes() {
        for c in [None, Some(2.0), Some(100.0)] {
            let s = (0..=10)
                .map(|k| {
                    log_likelihood(
                        &Observation::Counts {
                            methylated: k,
                            total: 10,
                            concentration: c,
                        },
                        0.3,
                    )
                    .unwrap()
                    .exp()
                })
                .sum::<f64>();
            assert!((s - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn missing_and_boundaries() {
        assert_eq!(log_likelihood(&Observation::Missing, 0.5).unwrap(), 0.0);
        assert_eq!(
            log_likelihood(&Observation::Binary(false), 0.0).unwrap(),
            0.0
        );
        assert!(log_likelihood(
            &Observation::Counts {
                methylated: 2,
                total: 1,
                concentration: None
            },
            0.5
        )
        .is_err());
    }
    #[test]
    fn optimizer_recovers() {
        let (s, x) = optimize(|p| (p[0] - 0.37).powi(2), &[vec![0.1]], &[(0.0, 1.0)]);
        assert!(s < 1e-8 && (x[0] - 0.37).abs() < 1e-4);
    }
    #[test]
    fn aggregate_boundary_noise() {
        let (d, v) = aggregate_moments(0.0, 0.1, 0.2, 100);
        assert_eq!(d, 0.1);
        assert_eq!(v, 0.001);
    }
}

#[derive(Clone, Debug)]
pub struct Site {
    pub target: bool,
    pub initial: f64,
    /// Independent youthful reference for identity distance, separate from the binary kinetic target.
    pub reference: f64,
    pub q: f64,
    /// Frozen mechanistic/context covariates; never derived from held-out aging outcomes.
    pub context: Vec<f64>,
    pub away: f64,
    pub restoration: f64,
    /// Multiplicative rate increase after a discrete maintenance-state switch.
    pub old_multiplier: f64,
    /// Low-rank log-rate factor loadings (M5).
    pub loadings: Vec<f64>,
}
#[derive(Clone, Debug)]
pub enum Maintenance {
    Fixed,
    Switching {
        rho: f64,
    },
    /// Explicit exploratory M4b: OU process on a log-rate modifier, initialized at zero.
    Continuous {
        reversion: f64,
        mean: f64,
        noise: f64,
    },
}
#[derive(Clone, Debug)]
pub struct Mechanism {
    pub protection: f64,
    pub recovery: f64,
    pub selection_baseline: f64,
    pub selection_strength: f64,
    pub maintenance: Maintenance,
    /// M5 factors are independent OU processes; shared across loci within a cell.
    pub factor_reversion: f64,
    pub factor_noise: f64,
}
impl Default for Mechanism {
    fn default() -> Self {
        Self {
            protection: 0.0,
            recovery: 0.0,
            selection_baseline: 0.0,
            selection_strength: 0.0,
            maintenance: Maintenance::Fixed,
            factor_reversion: 1.0,
            factor_noise: 0.0,
        }
    }
}
#[derive(Clone, Debug)]
pub struct Path {
    pub states: Vec<bool>,
    pub death: Option<f64>,
    pub first_passage: Option<f64>,
    pub switched: Option<f64>,
    pub endpoint_time: f64,
    /// None means an exact jump-process path; Some means grid-resolved first passage.
    pub time_resolution: Option<f64>,
}
fn validate(sites: &[Site], m: &Mechanism, horizon: f64, threshold: f64) -> Result<(), String> {
    if sites.is_empty()
        || !horizon.is_finite()
        || horizon < 0.0
        || !threshold.is_finite()
        || threshold < 0.0
    {
        return Err("invalid simulation horizon/threshold/sites".into());
    }
    for s in sites {
        if !(0.0..=1.0).contains(&s.initial)
            || !(0.0..=1.0).contains(&s.reference)
            || !(0.0..=1.0).contains(&s.q)
            || !s.away.is_finite()
            || s.away < 0.0
            || !s.restoration.is_finite()
            || s.restoration < 0.0
            || !s.old_multiplier.is_finite()
            || s.old_multiplier <= 0.0
            || s.context.iter().chain(&s.loadings).any(|v| !v.is_finite())
        {
            return Err("invalid site configuration".into());
        }
    }
    for x in [
        m.protection,
        m.recovery,
        m.selection_baseline,
        m.selection_strength,
        m.factor_reversion,
        m.factor_noise,
    ] {
        if !x.is_finite() || x < 0.0 {
            return Err("invalid mechanism parameter".into());
        }
    }
    match m.maintenance {
        Maintenance::Switching { rho } if !rho.is_finite() || rho < 0.0 => {
            return Err("invalid switching rate".into())
        }
        Maintenance::Continuous {
            reversion,
            mean,
            noise,
        } if !reversion.is_finite()
            || reversion < 0.0
            || !mean.is_finite()
            || !noise.is_finite()
            || noise < 0.0 =>
        {
            return Err("invalid OU maintenance parameters".into())
        }
        _ => {}
    }
    let rank = sites[0].loadings.len();
    if sites.iter().any(|s| s.loadings.len() != rank) {
        return Err("inconsistent factor rank".into());
    }
    Ok(())
}
pub fn distance(states: &[bool], sites: &[Site]) -> f64 {
    assert_eq!(states.len(), sites.len());
    states
        .iter()
        .zip(sites)
        .map(|(x, s)| s.q * (f64::from(*x) - s.reference).powi(2))
        .sum()
}
/// Selection is defined against the binary kinetic target; identity distance can
/// instead use a fractional, independently estimated youthful reference.
fn selection_deviation(states: &[bool], sites: &[Site]) -> f64 {
    states
        .iter()
        .zip(sites)
        .map(|(x, s)| s.q * f64::from(*x != s.target))
        .sum()
}
fn rates(s: &Site, m: &Mechanism, old: bool, latent: f64, factors: &[f64]) -> (f64, f64) {
    let modifier = (latent
        + s.loadings
            .iter()
            .zip(factors)
            .map(|(b, f)| b * f)
            .sum::<f64>())
    .exp()
        * if old { s.old_multiplier } else { 1.0 };
    (
        s.away * (-m.protection * s.q).exp() * modifier,
        s.restoration * (m.recovery * s.q).exp() * modifier,
    )
}
fn exponential(rng: &mut StdRng, rate: f64) -> f64 {
    if rate == 0.0 {
        f64::INFINITY
    } else {
        -rng.gen::<f64>().max(f64::MIN_POSITIVE).ln() / rate
    }
}
fn gaussian(rng: &mut StdRng) -> f64 {
    (-2.0 * rng.gen::<f64>().max(f64::MIN_POSITIVE).ln()).sqrt()
        * (2.0 * PI * rng.gen::<f64>()).cos()
}
fn ou_step(rng: &mut StdRng, x: f64, k: f64, mu: f64, eta: f64, dt: f64) -> f64 {
    let a = (-k * dt).exp();
    let v = if k == 0.0 {
        eta * eta * dt
    } else {
        eta * eta * (-(-2.0 * k * dt).exp_m1()) / (2.0 * k)
    };
    mu + (x - mu) * a + v.sqrt() * gaussian(rng)
}

/// Exact Gillespie simulation for M1/M2/M3/M4a. Includes state-dependent killing as
/// a competing event and records first passage immediately, including at time zero.
/// No clone proliferation is implied by effective removal.
pub fn simulate(
    sites: &[Site],
    m: &Mechanism,
    horizon: f64,
    threshold: f64,
    seed: u64,
) -> Result<Path, String> {
    validate(sites, m, horizon, threshold)?;
    if matches!(m.maintenance, Maintenance::Continuous { .. }) || m.factor_noise > 0.0 {
        return Err(
            "continuous maintenance/factors require simulate_grid and an explicit time step".into(),
        );
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let states = sites.iter().map(|s| rng.gen::<f64>() < s.initial).collect();
    let mut path = Path {
        states,
        death: None,
        first_passage: None,
        switched: None,
        endpoint_time: horizon,
        time_resolution: None,
    };
    let mut t = 0.0;
    if distance(&path.states, sites) > threshold {
        path.first_passage = Some(0.0);
    }
    loop {
        let old = path.switched.is_some();
        let mut events: Vec<f64> = sites
            .iter()
            .zip(&path.states)
            .map(|(s, x)| {
                let (a, b) = rates(s, m, old, 0.0, &[]);
                if *x == s.target {
                    a
                } else {
                    b
                }
            })
            .collect();
        let switch = match m.maintenance {
            Maintenance::Switching { rho } if !old => rho,
            _ => 0.0,
        };
        let kill = m.selection_baseline
            * (m.selection_strength * selection_deviation(&path.states, sites)).exp();
        events.push(switch);
        events.push(kill);
        let total = events.iter().sum::<f64>();
        if !total.is_finite() {
            return Err("rate overflow; reduce parameters or normalize fixed site weights".into());
        }
        t += exponential(&mut rng, total);
        if t > horizon {
            break;
        }
        let u = rng.gen::<f64>() * total;
        let mut cum = 0.0;
        let event = events
            .iter()
            .position(|r| {
                cum += r;
                u < cum
            })
            .ok_or("failed event sampling")?;
        if event < sites.len() {
            path.states[event] = !path.states[event];
            if path.first_passage.is_none() && distance(&path.states, sites) > threshold {
                path.first_passage = Some(t);
            }
        } else if event == sites.len() {
            path.switched = Some(t);
        } else {
            path.death = Some(t);
            path.endpoint_time = t;
            break;
        }
    }
    Ok(path)
}

/// M4b/M5 approximation: exact frozen-rate endpoint CTMC updates and exact OU
/// updates per time step. Killing and first passage are grid resolved. Run dt
/// sensitivity checks; excursions that return within a step may be missed.
pub fn simulate_grid(
    sites: &[Site],
    m: &Mechanism,
    horizon: f64,
    threshold: f64,
    dt: f64,
    seed: u64,
) -> Result<Path, String> {
    validate(sites, m, horizon, threshold)?;
    if !dt.is_finite() || dt <= 0.0 {
        return Err("invalid time step".into());
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let states = sites.iter().map(|s| rng.gen::<f64>() < s.initial).collect();
    let mut p = Path {
        states,
        death: None,
        first_passage: None,
        switched: None,
        endpoint_time: horizon,
        time_resolution: Some(dt),
    };
    let mut latent = 0.0;
    let mut factors = vec![0.0; sites[0].loadings.len()];
    let mut t = 0.0;
    if distance(&p.states, sites) > threshold {
        p.first_passage = Some(0.0);
    }
    while t < horizon {
        let step = dt.min(horizon - t);
        let hazard = m.selection_baseline
            * (m.selection_strength * selection_deviation(&p.states, sites)).exp();
        if !hazard.is_finite() {
            return Err("selection hazard overflow".into());
        }
        if rng.gen::<f64>() < -(-hazard * step).exp_m1() {
            p.death = Some(t + step);
            p.endpoint_time = t + step;
            break;
        }
        if let Maintenance::Switching { rho } = m.maintenance {
            if p.switched.is_none() && rng.gen::<f64>() < -(-rho * step).exp_m1() {
                p.switched = Some(t + step);
            }
        }
        for (s, x) in sites.iter().zip(&mut p.states) {
            let (a, b) = rates(s, m, p.switched.is_some(), latent, &factors);
            if !a.is_finite() || !b.is_finite() {
                return Err("latent rate overflow".into());
            }
            let (gain, loss) = if s.target { (b, a) } else { (a, b) };
            *x = rng.gen::<f64>() < probability(f64::from(*x), gain, loss, step);
        }
        t += step;
        if p.first_passage.is_none() && distance(&p.states, sites) > threshold {
            p.first_passage = Some(t);
        }
        if let Maintenance::Continuous {
            reversion,
            mean,
            noise,
        } = m.maintenance
        {
            latent = ou_step(&mut rng, latent, reversion, mean, noise, step);
        }
        for f in &mut factors {
            *f = ou_step(&mut rng, *f, m.factor_reversion, 0.0, m.factor_noise, step);
        }
    }
    Ok(p)
}

/// Exact one-locus killed-CTMC endpoint probabilities, conditional on survival.
/// Unlike an unnormalized selection weight, this accounts for exposure throughout the path.
pub fn survivor_probability(
    initial: f64,
    gain: f64,
    loss: f64,
    h0: f64,
    h1: f64,
    t: f64,
) -> Result<(f64, f64), String> {
    if !(0.0..=1.0).contains(&initial)
        || [gain, loss, h0, h1, t]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err("invalid killed CTMC".into());
    }
    let a = -gain - h0;
    let d = -loss - h1;
    let tr = (a + d) / 2.0;
    let delta = (((a - d) / 2.0).powi(2) + gain * loss).sqrt();
    let ep = ((tr + delta) * t).exp();
    let em = ((tr - delta) * t).exp();
    let c = (ep + em) / 2.0;
    let ss = if delta < 1e-12 {
        t * (tr * t).exp()
    } else {
        (ep - em) / (2.0 * delta)
    };
    let u0 = c * (1.0 - initial) + ss * ((a - tr) * (1.0 - initial) + loss * initial);
    let u1 = c * initial + ss * (gain * (1.0 - initial) + (d - tr) * initial);
    let survival = u0 + u1;
    if survival <= 0.0 || !survival.is_finite() {
        return Err("survival underflow".into());
    }
    Ok(((u1 / survival).clamp(0.0, 1.0), survival.clamp(0.0, 1.0)))
}

/// Hierarchical penalized composite-likelihood kinetic fitting. Site effects are shrunk to fixed context predictions.
/// Targets, initial methylation and covariates must come from independent/training references.
#[derive(Clone)]
pub struct KineticRow {
    pub site: usize,
    pub time: f64,
    pub replicate: String,
    pub observation: Observation,
}
#[derive(Clone, Debug)]
pub struct KineticFit {
    pub parameters: Vec<f64>,
    pub objective: f64,
    pub protection: bool,
    pub context_count: usize,
    pub sites: usize,
    pub shrinkage: f64,
}
pub fn kinetic_prediction(
    fit: &KineticFit,
    sites: &[Site],
    row: &KineticRow,
) -> Result<f64, String> {
    let s = sites.get(row.site).ok_or("unknown site")?;
    let d = fit.context_count;
    let p = &fit.parameters;
    let qoffset = 2 + 2 * d;
    let aq = if fit.protection { p[qoffset] } else { 0.0 };
    let bq = if fit.protection { p[qoffset + 1] } else { 0.0 };
    let offset = qoffset + if fit.protection { 2 } else { 0 };
    let a = p[0]
        + s.context
            .iter()
            .enumerate()
            .map(|(j, z)| p[2 + j] * z)
            .sum::<f64>()
        - aq * s.q
        + p[offset + 2 * row.site];
    let b = p[1]
        + s.context
            .iter()
            .enumerate()
            .map(|(j, z)| p[2 + d + j] * z)
            .sum::<f64>()
        + bq * s.q
        + p[offset + 2 * row.site + 1];
    let (gain, loss) = if s.target {
        (b.exp(), a.exp())
    } else {
        (a.exp(), b.exp())
    };
    Ok(probability(s.initial, gain, loss, row.time))
}
pub fn fit_kinetics(
    sites: &[Site],
    rows: &[KineticRow],
    protection: bool,
    shrinkage: f64,
) -> Result<KineticFit, String> {
    if sites.is_empty() || rows.is_empty() || !shrinkage.is_finite() || shrinkage <= 0.0 {
        return Err("invalid kinetic fit input".into());
    }
    validate(sites, &Mechanism::default(), 1.0, 1.0)?;
    let d = sites[0].context.len();
    if sites.iter().any(|s| s.context.len() != d)
        || rows
            .iter()
            .any(|r| r.site >= sites.len() || !r.time.is_finite() || r.time < 0.0)
    {
        return Err("invalid kinetic features/rows".into());
    }
    let globals = 2 + 2 * d + if protection { 2 } else { 0 };
    let mut bounds = vec![(-12.0, 2.0), (-12.0, 2.0)];
    bounds.extend(vec![(-3.0, 3.0); 2 * d]);
    if protection {
        bounds.extend([(0.0, 5.0), (0.0, 5.0)]);
    }
    bounds.extend(vec![(-4.0, 4.0); 2 * sites.len()]);
    let template = KineticFit {
        parameters: vec![],
        objective: 0.0,
        protection,
        context_count: d,
        sites: sites.len(),
        shrinkage,
    };
    let counts: BTreeMapForFit = rows
        .iter()
        .fold(std::collections::BTreeMap::new(), |mut m, r| {
            *m.entry(r.replicate.clone()).or_default() += 1;
            m
        });
    let starts: Vec<_> = [-5.0, -3.0]
        .iter()
        .map(|v| {
            let mut x = vec![0.0; bounds.len()];
            x[0] = *v;
            x[1] = *v;
            x
        })
        .collect();
    let (objective, parameters) = optimize(
        |p| {
            let mut fit = template.clone();
            fit.parameters = p.to_vec();
            let ll = rows
                .iter()
                .map(|r| {
                    kinetic_prediction(&fit, sites, r)
                        .and_then(|p| log_likelihood(&r.observation, p))
                        .unwrap_or(f64::NEG_INFINITY)
                        / counts[&r.replicate] as f64
                })
                .sum::<f64>();
            let penalty = p[globals..]
                .iter()
                .map(|v| v * v / (2.0 * shrinkage * shrinkage))
                .sum::<f64>();
            -ll + penalty
        },
        &starts,
        &bounds,
    );
    if !objective.is_finite() {
        return Err("nonfinite kinetic objective".into());
    }
    Ok(KineticFit {
        parameters,
        objective,
        ..template
    })
}
type BTreeMapForFit = std::collections::BTreeMap<String, usize>;

/// Empirical matched-mask null: permute weights only within prespecified strata.
/// Strata must jointly encode baseline bins, CpG density, sequence score, compartment,
/// chromosome etc. Matching quality and unsupported strata are the caller's responsibility.
pub fn matched_masks(
    q: &[f64],
    strata: &[String],
    seed: u64,
    count: usize,
) -> Result<Vec<Vec<f64>>, String> {
    use rand::seq::SliceRandom;
    if q.len() != strata.len() || q.iter().any(|v| !(0.0..=1.0).contains(v)) {
        return Err("invalid mask/strata".into());
    }
    let mut groups: std::collections::BTreeMap<&str, Vec<usize>> =
        std::collections::BTreeMap::new();
    for (i, s) in strata.iter().enumerate() {
        groups.entry(s).or_default().push(i);
    }
    let mut rng = StdRng::seed_from_u64(seed);
    Ok((0..count)
        .map(|_| {
            let mut mask = q.to_vec();
            for js in groups.values() {
                let mut vals: Vec<_> = js.iter().map(|j| q[*j]).collect();
                vals.shuffle(&mut rng);
                for (j, v) in js.iter().zip(vals) {
                    mask[*j] = v;
                }
            }
            mask
        })
        .collect())
}

#[cfg(test)]
mod stochastic_tests {
    use super::*;
    fn site() -> Site {
        Site {
            target: false,
            initial: 0.0,
            reference: 0.0,
            q: 1.0,
            context: vec![],
            away: 0.2,
            restoration: 0.1,
            old_multiplier: 4.0,
            loadings: vec![],
        }
    }
    #[test]
    fn gillespie_matches_analytic() {
        let s = vec![site()];
        let n = 12000;
        let endpoints = (0..n)
            .map(|seed| {
                f64::from(
                    simulate(&s, &Mechanism::default(), 3.0, 0.5, seed)
                        .unwrap()
                        .states[0],
                )
            })
            .sum::<f64>()
            / n as f64;
        assert!((endpoints - probability(0.0, 0.2, 0.1, 3.0)).abs() < 0.015);
    }
    #[test]
    fn killed_chain_matches_simulation() {
        let s = vec![site()];
        let m = Mechanism {
            selection_baseline: 0.05,
            selection_strength: 2.0,
            ..Mechanism::default()
        };
        let paths: Vec<_> = (0..16000)
            .map(|seed| simulate(&s, &m, 5.0, 0.5, seed).unwrap())
            .collect();
        let alive: Vec<_> = paths.iter().filter(|p| p.death.is_none()).collect();
        let expected =
            survivor_probability(0.0, 0.2, 0.1, 0.05, 0.05 * 2.0_f64.exp(), 5.0).unwrap();
        let actual = alive.iter().filter(|p| p.states[0]).count() as f64 / alive.len() as f64;
        assert!((actual - expected.0).abs() < 0.015);
        assert!((alive.len() as f64 / paths.len() as f64 - expected.1).abs() < 0.015);
    }
    #[test]
    fn first_passage_exponential() {
        let s = vec![site()];
        let n = 12000;
        let crossed = (0..n)
            .filter(|seed| {
                simulate(&s, &Mechanism::default(), 2.0, 0.5, *seed)
                    .unwrap()
                    .first_passage
                    .is_some()
            })
            .count() as f64
            / n as f64;
        assert!((crossed - (1.0 - (-0.4_f64).exp())).abs() < 0.015);
    }
    #[test]
    fn immediate_passage_and_reproducibility() {
        let s = vec![Site {
            initial: 1.0,
            ..site()
        }];
        assert_eq!(
            simulate(&s, &Mechanism::default(), 0.0, 0.5, 1)
                .unwrap()
                .first_passage,
            Some(0.0)
        );
        assert_eq!(
            simulate(&s, &Mechanism::default(), 3.0, 0.5, 2)
                .unwrap()
                .states,
            simulate(&s, &Mechanism::default(), 3.0, 0.5, 2)
                .unwrap()
                .states
        );
    }
    #[test]
    fn grid_and_latent_paths() {
        let s = vec![Site {
            loadings: vec![0.3],
            ..site()
        }];
        let m = Mechanism {
            maintenance: Maintenance::Continuous {
                reversion: 0.2,
                mean: 0.0,
                noise: 0.3,
            },
            factor_noise: 0.1,
            ..Mechanism::default()
        };
        assert!(simulate(&s, &m, 3.0, 0.5, 1).is_err());
        let p = simulate_grid(&s, &m, 3.0, 0.5, 0.01, 1).unwrap();
        assert_eq!(p.time_resolution, Some(0.01));
    }
    #[test]
    fn masks_keep_strata_counts() {
        let q = vec![1.0, 0.0, 1.0, 0.0];
        let strata = vec!["a".into(), "a".into(), "b".into(), "b".into()];
        for mask in matched_masks(&q, &strata, 3, 100).unwrap() {
            assert_eq!(mask[0] + mask[1], 1.0);
            assert_eq!(mask[2] + mask[3], 1.0);
        }
    }
    #[test]
    fn survivor_reduces_to_chain() {
        let p = survivor_probability(0.3, 0.2, 0.1, 0.05, 0.05, 3.0).unwrap();
        assert!((p.0 - probability(0.3, 0.2, 0.1, 3.0)).abs() < 1e-12);
        assert!((p.1 - (-0.15_f64).exp()).abs() < 1e-12);
    }
}

/// Generic simulation-based likelihood for the hierarchy. This is an explicitly
/// finite-particle approximation; common random numbers stabilize optimization.
/// Genomic features/site rates must be frozen before evaluating held-out subjects.
#[derive(Clone, Debug)]
pub enum ModelFamily {
    Independent,
    Protected,
    Selection,
    Switching,
    Continuous,
    Factors,
}
#[derive(Clone, Debug)]
pub enum SampleObservation {
    /// Single-cell monoallelic calls/counts with a fixed conversion/call error.
    Cell {
        sites: Vec<Observation>,
        call_error: f64,
    },
    /// Reads or arrays aggregated over surviving cells.
    Bulk(Vec<Observation>),
}
#[derive(Clone, Debug)]
pub struct Snapshot {
    pub time: f64,
    pub replicate: String,
    pub observation: SampleObservation,
}
#[derive(Clone, Debug)]
pub struct ParticlePrediction {
    pub log_likelihood: f64,
    pub survivors: usize,
    pub particles: usize,
    pub survivor_probabilities: Vec<f64>,
}
fn model_parameters(
    family: &ModelFamily,
    p: &[f64],
    sites: &[Site],
) -> Result<(Vec<Site>, Mechanism), String> {
    let expected = match family {
        ModelFamily::Independent => 2,
        _ => 4,
    };
    if p.len() != expected || p.iter().any(|x| !x.is_finite()) {
        return Err("incorrect model parameters".into());
    }
    let mut sites = sites.to_vec();
    for s in &mut sites {
        s.away *= p[0].exp();
        s.restoration *= p[1].exp();
    }
    let mut m = Mechanism::default();
    match family {
        ModelFamily::Independent => {}
        ModelFamily::Protected => {
            m.protection = p[2];
            m.recovery = p[3];
        }
        ModelFamily::Selection => {
            m.selection_baseline = p[2].exp();
            m.selection_strength = p[3];
        }
        ModelFamily::Switching => {
            m.maintenance = Maintenance::Switching { rho: p[2].exp() };
            for s in &mut sites {
                s.old_multiplier = p[3].exp();
            }
        }
        ModelFamily::Continuous => {
            m.maintenance = Maintenance::Continuous {
                reversion: 0.2,
                mean: p[2],
                noise: p[3],
            };
        }
        ModelFamily::Factors => {
            m.maintenance = Maintenance::Continuous {
                reversion: 0.2,
                mean: p[2],
                noise: 0.0,
            };
            m.factor_noise = p[3];
            if sites.iter().all(|s| s.loadings.is_empty()) {
                return Err("factor model requires fixed low-rank loadings".into());
            }
        }
    }
    Ok((sites, m))
}
pub fn particle_predict(
    family: &ModelFamily,
    parameters: &[f64],
    sites: &[Site],
    sample: &Snapshot,
    particles: usize,
    dt: f64,
    seed: u64,
) -> Result<ParticlePrediction, String> {
    if particles < 2 {
        return Err("at least two particles required".into());
    }
    let (sites, m) = model_parameters(family, parameters, sites)?;
    let obs = match &sample.observation {
        SampleObservation::Cell { sites, call_error } => {
            if !call_error.is_finite() || !(0.0..0.5).contains(call_error) {
                return Err("invalid call error".into());
            }
            sites
        }
        SampleObservation::Bulk(v) => v,
    };
    if obs.len() != sites.len() {
        return Err("observation and site dimensions differ".into());
    }
    let continuous = matches!(family, ModelFamily::Continuous | ModelFamily::Factors);
    let mut sums = vec![0.0; sites.len()];
    let mut cell_ll = vec![];
    let mut survivors = 0;
    for j in 0..particles {
        let particle_seed = seed.wrapping_add(j as u64 * 104729);
        let path = if continuous {
            simulate_grid(&sites, &m, sample.time, f64::MAX, dt, particle_seed)?
        } else {
            simulate(&sites, &m, sample.time, f64::MAX, particle_seed)?
        };
        if path.death.is_some() {
            continue;
        }
        survivors += 1;
        for (sum, x) in sums.iter_mut().zip(&path.states) {
            *sum += f64::from(*x);
        }
        if let SampleObservation::Cell { call_error, .. } = sample.observation {
            let ll = obs
                .iter()
                .zip(&path.states)
                .map(|(o, x)| log_likelihood(o, if *x { 1.0 - call_error } else { call_error }))
                .collect::<Result<Vec<_>, _>>()?
                .iter()
                .sum();
            cell_ll.push(ll);
        }
    }
    if survivors < 2 {
        return Err("too few surviving particles; increase particle count or reduce hazard".into());
    }
    let probabilities: Vec<_> = sums.iter().map(|x| x / survivors as f64).collect();
    let ll = match sample.observation {
        SampleObservation::Cell { .. } => logsumexp(&cell_ll) - (survivors as f64).ln(),
        SampleObservation::Bulk(_) => obs
            .iter()
            .zip(&probabilities)
            .map(|(o, p)| log_likelihood(o, *p))
            .collect::<Result<Vec<_>, _>>()?
            .iter()
            .sum(),
    };
    Ok(ParticlePrediction {
        log_likelihood: ll,
        survivors,
        particles,
        survivor_probabilities: probabilities,
    })
}
#[derive(Clone, Debug)]
pub struct ParticleFit {
    pub family: ModelFamily,
    pub parameters: Vec<f64>,
    pub objective: f64,
    pub particles: usize,
    pub dt: f64,
    pub seed: u64,
}
/// Penalized likelihood of approximate marginal observations. Neither posterior
/// inference nor proof of mechanism identifiability. Evaluate with fresh seeds and
/// larger particle counts and independent held-out biological replicates.
pub fn fit_particles(
    family: ModelFamily,
    sites: &[Site],
    samples: &[Snapshot],
    particles: usize,
    dt: f64,
    seed: u64,
) -> Result<ParticleFit, String> {
    if samples.is_empty() {
        return Err("empty fitting dataset".into());
    }
    let mut bounds = vec![(-2.0, 2.0); 2];
    let mut start = vec![0.0; 2];
    match family {
        ModelFamily::Independent => {}
        ModelFamily::Protected => {
            bounds.extend([(0.0, 3.0), (0.0, 3.0)]);
            start.extend([0.5, 0.5]);
        }
        ModelFamily::Selection => {
            bounds.extend([(-7.0, -1.0), (0.0, 3.0)]);
            start.extend([-4.0, 1.0]);
        }
        ModelFamily::Switching => {
            bounds.extend([(-7.0, -1.0), (0.0, 2.0)]);
            start.extend([-4.0, 1.0]);
        }
        ModelFamily::Continuous | ModelFamily::Factors => {
            bounds.extend([(-1.0, 1.0), (0.0, 1.0)]);
            start.extend([0.0, 0.1]);
        }
    }
    // Preflight validation prevents optimizing a permanently invalid configuration.
    particle_predict(&family, &start, sites, &samples[0], particles, dt, seed)?;
    let counts: BTreeMapForFit =
        samples
            .iter()
            .fold(std::collections::BTreeMap::new(), |mut m, s| {
                *m.entry(s.replicate.clone()).or_default() += 1;
                m
            });
    let (objective, parameters) = optimize(
        |p| {
            let mut ll = 0.0;
            for (j, s) in samples.iter().enumerate() {
                match particle_predict(
                    &family,
                    p,
                    sites,
                    s,
                    particles,
                    dt,
                    seed.wrapping_add(j as u64 * 10000019),
                ) {
                    Ok(pred) => ll += pred.log_likelihood / counts[&s.replicate] as f64,
                    Err(_) => return f64::INFINITY,
                }
            }
            -ll + 0.01 * p.iter().map(|x| x * x).sum::<f64>()
        },
        &[start],
        &bounds,
    );
    if !objective.is_finite() {
        return Err("nonfinite particle fit".into());
    }
    Ok(ParticleFit {
        family,
        parameters,
        objective,
        particles,
        dt,
        seed,
    })
}

/// Low-rank residual covariance via power iteration with deflation. Inputs must
/// already have confounders removed within training folds. Rank selection is left
/// to held-out likelihood; eigenvalues themselves are not proof of common fidelity.
pub fn residual_factors(
    residuals: &[Vec<f64>],
    rank: usize,
) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    if residuals.len() < 2 || residuals[0].is_empty() {
        return Err("empty residual matrix".into());
    }
    let d = residuals[0].len();
    if rank > d
        || residuals
            .iter()
            .any(|r| r.len() != d || r.iter().any(|v| !v.is_finite()))
    {
        return Err("invalid residual matrix/rank".into());
    }
    let means: Vec<_> = (0..d)
        .map(|j| residuals.iter().map(|r| r[j]).sum::<f64>() / residuals.len() as f64)
        .collect();
    let centered: Vec<Vec<f64>> = residuals
        .iter()
        .map(|r| r.iter().zip(&means).map(|(x, m)| x - m).collect())
        .collect();
    let mut loadings: Vec<Vec<f64>> = vec![];
    let mut eigenvalues = vec![];
    for k in 0..rank {
        let mut v: Vec<_> = (0..d).map(|j| (((j + 1) * (k + 1)) as f64).sin()).collect();
        for _ in 0..200 {
            let mut w = vec![0.0; d];
            for r in &centered {
                let dot = r.iter().zip(&v).map(|(a, b)| a * b).sum::<f64>();
                for (w, r) in w.iter_mut().zip(r) {
                    *w += r * dot / (centered.len() - 1) as f64;
                }
            }
            for old in &loadings {
                let dot = w.iter().zip(old).map(|(a, b)| a * b).sum::<f64>();
                for (w, o) in w.iter_mut().zip(old) {
                    *w -= dot * o;
                }
            }
            let norm = w.iter().map(|v| v * v).sum::<f64>().sqrt();
            if norm < 1e-12 {
                v.fill(0.0);
                break;
            }
            for w in &mut w {
                *w /= norm;
            }
            let change = w.iter().zip(&v).map(|(a, b)| (a - b).abs()).sum::<f64>();
            v = w;
            if change < 1e-9 {
                break;
            }
        }
        let eigen = centered
            .iter()
            .map(|r| r.iter().zip(&v).map(|(a, b)| a * b).sum::<f64>().powi(2))
            .sum::<f64>()
            / (centered.len() - 1) as f64;
        eigenvalues.push(eigen);
        loadings.push(v);
    }
    Ok((eigenvalues, loadings))
}

#[cfg(test)]
mod inference_tests {
    use super::*;
    fn site(q: f64) -> Site {
        Site {
            target: false,
            initial: 0.0,
            reference: 0.0,
            q,
            context: vec![],
            away: 0.1,
            restoration: 0.2,
            old_multiplier: 4.0,
            loadings: vec![1.0],
        }
    }
    #[test]
    fn protection_recovery_synthetic() {
        let sites = vec![site(0.0), site(1.0)];
        let mut rows = vec![];
        for (i, s) in sites.iter().enumerate() {
            for time in [1.0, 5.0, 15.0] {
                let p = probability(
                    0.0,
                    s.away * (-1.2 * s.q).exp(),
                    s.restoration * (0.7 * s.q).exp(),
                    time,
                );
                rows.push(KineticRow {
                    site: i,
                    time,
                    replicate: "training".into(),
                    observation: Observation::Array { beta: p, sd: 0.005 },
                });
            }
        }
        let fit = fit_kinetics(&sites, &rows, true, 0.1).unwrap();
        for r in &rows {
            let Observation::Array { beta, .. } = r.observation else {
                unreachable!()
            };
            assert!((kinetic_prediction(&fit, &sites, r).unwrap() - beta).abs() < 0.01);
        }
    }
    #[test]
    fn all_particle_families_run() {
        for family in [
            ModelFamily::Independent,
            ModelFamily::Protected,
            ModelFamily::Selection,
            ModelFamily::Switching,
            ModelFamily::Continuous,
            ModelFamily::Factors,
        ] {
            let p = match family {
                ModelFamily::Independent => vec![0.0, 0.0],
                ModelFamily::Selection | ModelFamily::Switching => vec![0.0, 0.0, -4.0, 1.0],
                _ => vec![0.0, 0.0, 0.1, 0.1],
            };
            let sample = Snapshot {
                time: 2.0,
                replicate: "mouse".into(),
                observation: SampleObservation::Bulk(vec![Observation::Array {
                    beta: 0.2,
                    sd: 0.1,
                }]),
            };
            assert!(
                particle_predict(&family, &p, &[site(0.5)], &sample, 1000, 0.05, 7)
                    .unwrap()
                    .log_likelihood
                    .is_finite()
            );
        }
    }
    #[test]
    fn residual_rank_one() {
        let x = vec![vec![1.0, 2.0], vec![2.0, 4.0], vec![3.0, 6.0]];
        let (e, _) = residual_factors(&x, 2).unwrap();
        assert!((e[0] - 5.0).abs() < 1e-8);
        assert!(e[1] < 1e-8);
    }
}

/// Distance for complete fixed-universe observations. Missingness requires a
/// separate model; silently dropping sites changes the first-passage boundary.
pub fn identity_distance(
    observed: &[Option<f64>],
    reference: &[f64],
    weights: &[f64],
) -> Result<f64, String> {
    if observed.len() != reference.len() || observed.len() != weights.len() || observed.is_empty() {
        return Err("invalid identity metric dimensions".into());
    }
    observed
        .iter()
        .zip(reference)
        .zip(weights)
        .try_fold(0.0, |sum, ((x, r), q)| {
            let x = x.ok_or("missing identity site; imputation/coverage model required")?;
            if !(0.0..=1.0).contains(&x) || !(0.0..=1.0).contains(r) || !(0.0..=1.0).contains(q) {
                return Err("invalid identity metric input".into());
            }
            Ok(sum + q * (x - r).powi(2))
        })
}
#[cfg(test)]
mod validation_tests {
    use super::*;
    #[test]
    fn identity_fractional_reference() {
        assert!(
            (identity_distance(&[Some(0.5), Some(1.0)], &[0.5, 0.0], &[1.0, 0.5]).unwrap() - 0.5)
                .abs()
                < 1e-12
        );
        assert!(identity_distance(&[None], &[0.5], &[1.0]).is_err());
    }
    #[test]
    fn particle_cell_likelihood_matches_analytic() {
        let s = Site {
            target: false,
            initial: 0.0,
            reference: 0.0,
            q: 0.0,
            context: vec![],
            away: 0.2,
            restoration: 0.1,
            old_multiplier: 1.0,
            loadings: vec![],
        };
        let sample = Snapshot {
            time: 3.0,
            replicate: "animal".into(),
            observation: SampleObservation::Cell {
                sites: vec![Observation::Binary(true)],
                call_error: 0.01,
            },
        };
        let pred = particle_predict(
            &ModelFamily::Independent,
            &[0.0, 0.0],
            &[s],
            &sample,
            20000,
            0.1,
            42,
        )
        .unwrap();
        let expected = 0.01 + 0.98 * probability(0.0, 0.2, 0.1, 3.0);
        assert!((pred.log_likelihood.exp() - expected).abs() < 0.015);
    }
    #[test]
    fn particle_fit_predicts_independent_snapshot() {
        let s = Site {
            target: false,
            initial: 0.0,
            reference: 0.0,
            q: 0.0,
            context: vec![],
            away: 0.2,
            restoration: 0.1,
            old_multiplier: 1.0,
            loadings: vec![],
        };
        let samples: Vec<_> = [1.0, 5.0]
            .into_iter()
            .map(|time| Snapshot {
                time,
                replicate: "train".into(),
                observation: SampleObservation::Bulk(vec![Observation::Array {
                    beta: probability(0.0, 0.1, 0.15, time),
                    sd: 0.05,
                }]),
            })
            .collect();
        let fit = fit_particles(
            ModelFamily::Independent,
            std::slice::from_ref(&s),
            &samples,
            1000,
            0.05,
            7,
        )
        .unwrap();
        let pred = particle_predict(
            &fit.family,
            &fit.parameters,
            &[s],
            &samples[1],
            20000,
            0.05,
            111,
        )
        .unwrap();
        let expected = probability(0.0, 0.1, 0.15, 5.0);
        assert!((pred.survivor_probabilities[0] - expected).abs() < 0.08);
    }
}
