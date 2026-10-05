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

// Revision 3: context/exposure-only molecular transitions, with state fitness
// entering demography. The older Site/Mechanism APIs above are comparators.
#[derive(Clone, Debug)]
pub struct MolecularSite {
    /// Effective rates calibrated from molecular context, never functional cost.
    pub gain: f64,
    pub loss: f64,
}
#[derive(Clone, Debug)]
pub struct ExposureEpoch {
    pub end: f64,
    pub gain_multiplier: f64,
    pub loss_multiplier: f64,
}
#[derive(Clone, Debug)]
pub struct PopulationModel {
    pub sites: Vec<MolecularSite>,
    pub epochs: Vec<ExposureEpoch>,
    /// State-indexed rates; binary bit i is the state of site i.
    pub birth: Vec<f64>,
    pub death: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct PopulationDecomposition {
    pub methylation: Vec<f64>,
    pub intrinsic_change: Vec<f64>,
    pub selection_change: Vec<f64>,
    pub mean_net_growth: f64,
}
#[derive(Clone, Debug)]
pub struct ClonePopulation {
    /// Rows retain the supplied independent clone labels; columns are states.
    pub counts: Vec<Vec<u64>>,
    pub events: usize,
    pub extinction_time: Option<f64>,
    pub end_time: f64,
}
/// Signed costs are external effects of departures from the specified fitness
/// reference. Negative costs allow beneficial drift; this is not identity loss.
pub fn state_fitness(
    reference: &[bool],
    costs: &[f64],
    baseline: f64,
    strength: f64,
) -> Result<Vec<f64>, String> {
    if reference.is_empty()
        || reference.len() > 12
        || reference.len() != costs.len()
        || costs.iter().any(|x| !x.is_finite())
        || !baseline.is_finite()
        || !strength.is_finite()
        || strength < 0.0
    {
        return Err("invalid external state-fitness definition".into());
    }
    let fitness: Vec<_> = (0..1 << reference.len())
        .map(|s| {
            baseline
                - strength
                    * costs
                        .iter()
                        .enumerate()
                        .map(|(i, q)| q * f64::from(((s >> i) & 1 == 1) != reference[i]))
                        .sum::<f64>()
        })
        .collect();
    if fitness.iter().any(|r| !r.is_finite()) {
        return Err("state fitness overflow".into());
    }
    Ok(fitness)
}
/// Net growth does not identify turnover. This explicit nonnegative decomposition
/// is a modeling choice; alternative birth/death pairs can share these means.
pub fn state_demography(fitness: &[f64], turnover: f64) -> Result<(Vec<f64>, Vec<f64>), String> {
    if fitness.is_empty()
        || fitness.iter().any(|r| !r.is_finite())
        || !turnover.is_finite()
        || turnover < 0.0
    {
        return Err("invalid state demography".into());
    }
    let birth: Vec<_> = fitness.iter().map(|r| turnover + r.max(0.0)).collect();
    let death: Vec<_> = fitness.iter().map(|r| turnover + (-r).max(0.0)).collect();
    if birth.iter().chain(&death).any(|r| !r.is_finite()) {
        return Err("demographic rate overflow".into());
    }
    Ok((birth, death))
}
fn validate_population(
    model: &PopulationModel,
    counts: &[f64],
    horizon: f64,
) -> Result<(), String> {
    if model.sites.is_empty()
        || model.sites.len() > 12
        || !horizon.is_finite()
        || horizon < 0.0
        || counts.len() != 1 << model.sites.len()
        || counts.iter().any(|n| !n.is_finite() || *n < 0.0)
        || model.birth.len() != counts.len()
        || model.death.len() != counts.len()
        || model
            .birth
            .iter()
            .chain(&model.death)
            .any(|r| !r.is_finite() || *r < 0.0)
        || model
            .sites
            .iter()
            .any(|s| !s.gain.is_finite() || s.gain < 0.0 || !s.loss.is_finite() || s.loss < 0.0)
        || model.epochs.is_empty()
    {
        return Err("invalid finite-state population model".into());
    }
    let mut start = 0.0;
    for e in &model.epochs {
        if !e.end.is_finite()
            || e.end <= start
            || !e.gain_multiplier.is_finite()
            || e.gain_multiplier < 0.0
            || !e.loss_multiplier.is_finite()
            || e.loss_multiplier < 0.0
            || model.sites.iter().any(|s| {
                !(s.gain * e.gain_multiplier).is_finite()
                    || !(s.loss * e.loss_multiplier).is_finite()
                    || !(s.gain * e.gain_multiplier + s.loss * e.loss_multiplier).is_finite()
            })
        {
            return Err("invalid molecular exposure schedule".into());
        }
        start = e.end;
    }
    if start < horizon {
        return Err("exposure schedule ends before prediction horizon".into());
    }
    Ok(())
}
fn population_epoch(model: &PopulationModel, time: f64) -> &ExposureEpoch {
    model
        .epochs
        .iter()
        .find(|e| time < e.end)
        .unwrap_or_else(|| model.epochs.last().unwrap())
}
/// Derivative of expected absolute counts under transitions plus birth/death.
pub fn population_rhs(
    model: &PopulationModel,
    counts: &[f64],
    time: f64,
) -> Result<Vec<f64>, String> {
    validate_population(model, counts, time)?;
    let epoch = population_epoch(model, time);
    let mut derivative: Vec<_> = counts
        .iter()
        .enumerate()
        .map(|(s, n)| n * (model.birth[s] - model.death[s]))
        .collect();
    for (s, n) in counts.iter().enumerate() {
        for (i, site) in model.sites.iter().enumerate() {
            let rate = if (s >> i) & 1 == 0 {
                site.gain * epoch.gain_multiplier
            } else {
                site.loss * epoch.loss_multiplier
            };
            derivative[s] -= n * rate;
            derivative[s ^ (1 << i)] += n * rate;
        }
    }
    if derivative.iter().any(|v| !v.is_finite()) {
        return Err("population derivative overflow".into());
    }
    Ok(derivative)
}
/// Exact instantaneous Price decomposition for normalized expected composition.
/// It is not the expected composition of a finite population conditional on survival.
pub fn population_decomposition(
    model: &PopulationModel,
    counts: &[f64],
    time: f64,
) -> Result<PopulationDecomposition, String> {
    validate_population(model, counts, time)?;
    let total: f64 = counts.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err("composition undefined for an empty population".into());
    }
    let epoch = population_epoch(model, time);
    let p: Vec<_> = counts.iter().map(|n| n / total).collect();
    let growth: Vec<_> = model
        .birth
        .iter()
        .zip(&model.death)
        .map(|(b, d)| b - d)
        .collect();
    let mean_net_growth = mean_weighted(&p, &growth);
    let mut result = PopulationDecomposition {
        methylation: vec![],
        intrinsic_change: vec![],
        selection_change: vec![],
        mean_net_growth,
    };
    for (i, site) in model.sites.iter().enumerate() {
        let m: f64 = p
            .iter()
            .enumerate()
            .map(|(s, p)| p * f64::from((s >> i) & 1 == 1))
            .sum();
        let covariance: f64 = p
            .iter()
            .enumerate()
            .map(|(s, p)| p * (f64::from((s >> i) & 1 == 1) - m) * (growth[s] - mean_net_growth))
            .sum();
        result.methylation.push(m);
        result.intrinsic_change.push(
            site.gain * epoch.gain_multiplier * (1.0 - m) - site.loss * epoch.loss_multiplier * m,
        );
        result.selection_change.push(covariance);
    }
    Ok(result)
}
fn mean_weighted(weights: &[f64], values: &[f64]) -> f64 {
    weights.iter().zip(values).map(|(p, v)| p * v).sum()
}
/// Positivity-preserving Strang splitting of expected counts. Molecular updates
/// are exact conditional independent-site CTMC transitions within each epoch;
/// selection/transition coupling has O(dt^2) global error. Refine dt to verify.
pub fn population_expectation(
    model: &PopulationModel,
    initial: &[f64],
    horizon: f64,
    dt: f64,
) -> Result<Vec<f64>, String> {
    validate_population(model, initial, horizon)?;
    if !dt.is_finite() || dt <= 0.0 {
        return Err("invalid population solver step".into());
    }
    let mut counts = initial.to_vec();
    let mut time = 0.0;
    while time < horizon {
        let epoch = population_epoch(model, time);
        let boundary = epoch.end.min(horizon);
        let step = dt.min(boundary - time);
        if step <= 0.0 || time + step == time {
            return Err("population solver time step underflow".into());
        }
        let growth: Vec<_> = model
            .birth
            .iter()
            .zip(&model.death)
            .map(|(b, d)| ((b - d) * step / 2.0).exp())
            .collect();
        for (n, g) in counts.iter_mut().zip(&growth) {
            *n *= g;
        }
        for (i, site) in model.sites.iter().enumerate() {
            let gain = site.gain * epoch.gain_multiplier;
            let loss = site.loss * epoch.loss_multiplier;
            let from0 = probability(0.0, gain, loss, step);
            let from1 = probability(1.0, gain, loss, step);
            for s in 0..counts.len() {
                if (s >> i) & 1 == 0 {
                    let other = s | (1 << i);
                    let unmethylated = counts[s];
                    let methylated = counts[other];
                    counts[s] = unmethylated * (1.0 - from0) + methylated * (1.0 - from1);
                    counts[other] = unmethylated * from0 + methylated * from1;
                }
            }
        }
        for (n, g) in counts.iter_mut().zip(&growth) {
            *n *= g;
        }
        if counts.iter().any(|n| !n.is_finite()) {
            return Err("expected population overflow".into());
        }
        time = if step == boundary - time {
            boundary
        } else {
            time + step
        };
    }
    Ok(counts)
}
/// Exact branching CTMC under piecewise-constant exogenous molecular exposures.
/// Birth preserves parental state and clone label; a molecular jump can change
/// descendants later. Replication-linked strand-error inheritance is not modeled.
/// Limits return errors rather than silently reporting truncated paths as predictions.
pub fn simulate_population(
    model: &PopulationModel,
    initial: &[Vec<u64>],
    horizon: f64,
    seed: u64,
    event_limit: usize,
    population_limit: u64,
) -> Result<ClonePopulation, String> {
    if initial.is_empty() || event_limit == 0 || population_limit == 0 {
        return Err("invalid branching simulation limits".into());
    }
    let states = model.birth.len();
    if initial.iter().any(|r| r.len() != states) {
        return Err("invalid clone state count matrix".into());
    }
    let mut population = initial
        .iter()
        .flatten()
        .try_fold(0_u64, |a, b| a.checked_add(*b))
        .ok_or("initial population overflow")?;
    if population > population_limit {
        return Err("initial population exceeds limit".into());
    }
    let aggregate: Vec<_> = (0..states)
        .map(|s| initial.iter().map(|r| r[s] as f64).sum())
        .collect();
    validate_population(model, &aggregate, horizon)?;
    let mut rng = StdRng::seed_from_u64(seed);
    let mut counts = initial.to_vec();
    let mut time = 0.0;
    let mut events = 0;
    let mut extinction_time = if population == 0 { Some(0.0) } else { None };
    while time < horizon && population > 0 {
        let epoch = population_epoch(model, time);
        let boundary = epoch.end.min(horizon);
        let mut total = 0.0;
        for row in &counts {
            for (s, n) in row.iter().enumerate() {
                let mutations: f64 = model
                    .sites
                    .iter()
                    .enumerate()
                    .map(|(i, site)| {
                        if (s >> i) & 1 == 0 {
                            site.gain * epoch.gain_multiplier
                        } else {
                            site.loss * epoch.loss_multiplier
                        }
                    })
                    .sum();
                total += *n as f64 * (mutations + model.birth[s] + model.death[s]);
            }
        }
        if !total.is_finite() {
            return Err("branching event rate overflow".into());
        }
        let next = time + exponential(&mut rng, total);
        if next >= boundary {
            time = boundary;
            continue;
        }
        if next <= time {
            return Err("branching event time underflow".into());
        }
        time = next;
        if events == event_limit {
            return Err("branching event limit reached".into());
        }
        let mut draw = rng.gen::<f64>() * total;
        let mut chosen = None;
        'find: for (clone, row) in counts.iter().enumerate() {
            for (s, n) in row.iter().enumerate() {
                if *n == 0 {
                    continue;
                }
                for event in 0..model.sites.len() + 2 {
                    let rate = if event < model.sites.len() {
                        let site = &model.sites[event];
                        if (s >> event) & 1 == 0 {
                            site.gain * epoch.gain_multiplier
                        } else {
                            site.loss * epoch.loss_multiplier
                        }
                    } else if event == model.sites.len() {
                        model.birth[s]
                    } else {
                        model.death[s]
                    };
                    let mass = *n as f64 * rate;
                    if draw < mass {
                        chosen = Some((clone, s, event));
                        break 'find;
                    }
                    draw -= mass;
                }
            }
        }
        let (clone, state, event) = chosen.ok_or("failed branching event sampling")?;
        if event < model.sites.len() {
            counts[clone][state] -= 1;
            let target = state ^ (1 << event);
            counts[clone][target] = counts[clone][target]
                .checked_add(1)
                .ok_or("clone count overflow")?;
        } else if event == model.sites.len() {
            if population == population_limit {
                return Err("branching population limit reached".into());
            }
            counts[clone][state] = counts[clone][state]
                .checked_add(1)
                .ok_or("clone count overflow")?;
            population += 1;
        } else {
            counts[clone][state] -= 1;
            population -= 1;
            if population == 0 {
                extinction_time = Some(time);
            }
        }
        events += 1;
    }
    Ok(ClonePopulation {
        counts,
        events,
        extinction_time,
        end_time: horizon,
    })
}

#[cfg(test)]
mod population_tests {
    use super::*;
    fn model(gain: f64, loss: f64, birth: Vec<f64>, death: Vec<f64>) -> PopulationModel {
        PopulationModel {
            sites: vec![MolecularSite { gain, loss }],
            epochs: vec![ExposureEpoch {
                end: 10.0,
                gain_multiplier: 1.0,
                loss_multiplier: 1.0,
            }],
            birth,
            death,
        }
    }
    #[test]
    fn context_and_exposure_reduce_to_exact_ctmc_with_common_growth() {
        let mut m = model(0.2, 0.1, vec![0.3; 2], vec![0.1; 2]);
        m.epochs = vec![
            ExposureEpoch {
                end: 1.0,
                gain_multiplier: 0.5,
                loss_multiplier: 2.0,
            },
            ExposureEpoch {
                end: 4.0,
                gain_multiplier: 2.0,
                loss_multiplier: 0.5,
            },
        ];
        let counts = population_expectation(&m, &[80.0, 20.0], 3.0, 0.07).unwrap();
        let expected = probability(probability(0.2, 0.1, 0.2, 1.0), 0.4, 0.05, 2.0);
        assert!((counts.iter().sum::<f64>() - 100.0 * 0.6_f64.exp()).abs() < 1e-9);
        assert!((counts[1] / counts.iter().sum::<f64>() - expected).abs() < 1e-12);
    }
    #[test]
    fn selection_changes_bulk_without_molecular_repair_and_obeys_price() {
        let fitness = state_fitness(&[true], &[1.0], 0.2, 0.3).unwrap();
        let (birth, death) = state_demography(&fitness, 0.1).unwrap();
        let m = model(0.0, 0.0, birth, death);
        let decomposition = population_decomposition(&m, &[100.0, 100.0], 0.0).unwrap();
        assert_eq!(decomposition.intrinsic_change, vec![0.0]);
        assert!((decomposition.selection_change[0] - 0.075).abs() < 1e-12);
        let counts = population_expectation(&m, &[100.0, 100.0], 5.0, 0.05).unwrap();
        let expected = 1.0 / (1.0 + (-1.5_f64).exp());
        assert!((counts[1] / counts.iter().sum::<f64>() - expected).abs() < 1e-12);
        let rhs = population_rhs(&m, &[100.0, 100.0], 0.0).unwrap();
        let derivative = rhs[1] / 200.0 - 0.5 * rhs.iter().sum::<f64>() / 200.0;
        assert!((derivative - decomposition.selection_change[0]).abs() < 1e-12);
        // Fitness can favor drift; biological identity value is not clone fitness.
        let beneficial = state_fitness(&[true], &[-1.0], 0.2, 0.3).unwrap();
        assert!(beneficial[0] > beneficial[1]);
    }
    #[test]
    fn fitness_does_not_enter_intrinsic_transition_equation() {
        let a = model(0.2, 0.1, vec![0.1; 2], vec![0.0; 2]);
        let b = model(0.2, 0.1, vec![0.8, 0.1], vec![0.1, 0.3]);
        let x = population_decomposition(&a, &[20.0, 80.0], 0.0).unwrap();
        let y = population_decomposition(&b, &[20.0, 80.0], 0.0).unwrap();
        assert_eq!(x.intrinsic_change, y.intrinsic_change);
        assert_ne!(x.selection_change, y.selection_change);
    }
    #[test]
    fn expected_count_solver_converges_to_exact_selected_chain() {
        let m = model(0.3, 0.1, vec![0.2; 2], vec![0.05, 0.3]);
        let exact = survivor_probability(0.25, 0.3, 0.1, 0.05, 0.3, 3.0)
            .unwrap()
            .0;
        let error = |dt| {
            let n = population_expectation(&m, &[75.0, 25.0], 3.0, dt).unwrap();
            (n[1] / n.iter().sum::<f64>() - exact).abs()
        };
        assert!(error(0.05) < error(0.1) / 3.5);
        assert!(error(0.01) < 1e-6);
    }
    #[test]
    fn branching_counts_match_expectation_and_preserve_clone_labels() {
        let m = model(0.2, 0.1, vec![0.15; 2], vec![0.05; 2]);
        let expected = population_expectation(&m, &[50.0, 50.0], 2.0, 0.01).unwrap();
        let mut totals = [0.0; 2];
        for seed in 0..300 {
            let path =
                simulate_population(&m, &[vec![50, 0], vec![0, 50]], 2.0, seed, 10000, 10000)
                    .unwrap();
            assert_eq!(path.counts.len(), 2);
            for (s, total) in totals.iter_mut().enumerate() {
                *total += path.counts.iter().map(|r| r[s] as f64).sum::<f64>();
            }
        }
        for (observed, expected) in totals.iter().zip(expected) {
            assert!((observed / 300.0 - expected).abs() < 2.0);
        }
        let no_transition = model(0.0, 0.0, vec![0.15; 2], vec![0.05; 2]);
        let path = simulate_population(
            &no_transition,
            &[vec![20, 0], vec![0, 20]],
            2.0,
            9,
            10000,
            10000,
        )
        .unwrap();
        assert_eq!(path.counts[0][1], 0);
        assert_eq!(path.counts[1][0], 0);
    }
    #[test]
    fn extinction_and_simulation_limits_are_explicit() {
        let death = model(0.0, 0.0, vec![0.0; 2], vec![4.0; 2]);
        let path = simulate_population(&death, &[vec![1, 0]], 5.0, 1, 100, 100).unwrap();
        assert_eq!(path.counts, vec![vec![0, 0]]);
        assert!(path.extinction_time.is_some());
        let birth = model(0.0, 0.0, vec![100.0; 2], vec![0.0; 2]);
        assert!(simulate_population(&birth, &[vec![1, 0]], 5.0, 1, 100, 1).is_err());
        assert!(population_expectation(&birth, &[1.0, 0.0], 11.0, 0.1).is_err());
        assert!(population_decomposition(&death, &[0.0, 0.0], 0.0).is_err());
    }
}
