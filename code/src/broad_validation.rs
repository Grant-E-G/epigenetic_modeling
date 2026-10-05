//! Independent editing consequences and intervention-specific held-out kinetics.
use super::*;
use crate::validation::quantile;

#[derive(Clone, Default)]
struct Region {
    gene: String,
    chrom: String,
    start: u64,
    end: u64,
    density: f64,
    // One entry per genomic CpG, preserving the early-defined measurement set.
    counts: BTreeMap<u64, Vec<(f64, f64)>>,
}
#[derive(Clone)]
struct Well {
    plate: String,
    batch: String,
    cell: String,
    construct: String,
    gene: String,
    control: bool,
    signal: f64,
}
fn screen(out: &str) -> Result<BTreeMap<String, f64>> {
    let mut wells = Vec::new();
    for line in reader("data/derived/broad_screen.tsv")?.lines().skip(1) {
        let f = fields(&line?, '\t');
        let signal: f64 = f[8].parse()?;
        if signal <= 0.0 {
            return Err("nonpositive viability signal".into());
        }
        wells.push(Well {
            plate: f[0].clone(),
            batch: f[2].clone(),
            cell: f[4].clone(),
            construct: f[5].clone(),
            gene: f[7].clone(),
            control: f[3] == "neg",
            signal,
        });
    }
    let mut controls: BTreeMap<(&str, &str, &str), Vec<f64>> = BTreeMap::new();
    for w in &wells {
        if w.control {
            controls
                .entry((&w.cell, &w.construct, &w.plate))
                .or_default()
                .push(w.signal.ln());
        }
    }
    let mut targets: BTreeMap<(&str, &str, &str, &str), Vec<f64>> = BTreeMap::new();
    for w in &wells {
        if !w.control {
            let reference = controls
                .get(&(&*w.cell, &*w.construct, &*w.plate))
                .ok_or("missing plate controls")?;
            targets
                .entry((&w.cell, &w.gene, &w.construct, &w.batch))
                .or_default()
                .push(w.signal.ln() - mean(reference));
        }
    }
    let cells = ["HCT116", "DLD1"];
    let genes: BTreeSet<_> = wells
        .iter()
        .filter(|w| !w.control)
        .map(|w| w.gene.as_str())
        .collect();
    let mut summary = writer(&format!("{out}/broad_functional_scores.csv"))?;
    writeln!(summary,"cell,gene,batches,signed_cost,absolute_cost,positive_batches,vp64_log_effect,active_log_effect,inactive_log_effect")?;
    let mut batch_output = writer(&format!("{out}/broad_screen_batches.csv"))?;
    writeln!(batch_output,"cell,gene,batch,active_log_effect,inactive_log_effect,vp64_log_effect,signed_cost,active_wells,inactive_wells")?;
    let mut scores = BTreeMap::new();
    for cell in cells {
        for gene in &genes {
            let mut costs = vec![];
            let mut effects = [vec![], vec![], vec![]];
            for batch in ["1", "2", "3"] {
                let rows: Vec<_> = ["TET1", "TET1_IM", "VP64"]
                    .iter()
                    .map(|construct| targets.get(&(cell, *gene, *construct, batch)))
                    .collect();
                if rows.iter().any(|r| r.is_none()) {
                    continue;
                }
                let values: Vec<_> = rows.iter().map(|r| mean(r.unwrap())).collect();
                let cost = values[1] - values[0];
                costs.push(cost);
                for i in 0..3 {
                    effects[i].push(values[i]);
                }
                writeln!(
                    batch_output,
                    "{cell},{gene},{batch},{},{},{},{cost},{},{}",
                    values[0],
                    values[1],
                    values[2],
                    rows[0].unwrap().len(),
                    rows[1].unwrap().len()
                )?;
            }
            if costs.is_empty() {
                return Err("target without matched construct batches".into());
            }
            let cost = mean(&costs);
            writeln!(
                summary,
                "{cell},{gene},{},{cost},{},{},{},{},{}",
                costs.len(),
                cost.abs(),
                costs.iter().filter(|x| **x > 0.0).count(),
                mean(&effects[2]),
                mean(&effects[0]),
                mean(&effects[1])
            )?;
            if cell == "HCT116" {
                scores.insert(gene.to_string(), cost);
            }
        }
    }
    Ok(scores)
}
fn region_beta(region: &Region, sites: &[u64], sample: usize, coverage: f64) -> Option<f64> {
    let values: Vec<_> = sites
        .iter()
        .filter_map(|site| {
            let (m, n) = region.counts[site][sample];
            (n >= coverage).then_some(m / n)
        })
        .collect();
    // Fixed early sites: later coverage failure is missing, never unmethylated.
    (values.len() >= 3 && values.len() as f64 >= 0.8 * sites.len() as f64).then(|| mean(&values))
}
#[derive(Clone)]
struct Recovery {
    gene: String,
    features: Vec<f64>,
    cost: f64,
    outcomes: [Option<f64>; 2],
}
fn recovery_predictions(
    rows: &[Recovery],
    out: &str,
    label: &str,
    days: [i32; 2],
) -> Result<String> {
    let mut errors: BTreeMap<(usize, usize, usize), Vec<f64>> = BTreeMap::new();
    // Keep the primary forecasts reviewable; regenerate exploratory/gate tables
    // under ignored data rather than fragmenting the versioned results directory.
    let directory = if label == "primary" {
        out
    } else {
        "data/derived"
    };
    let mut w = writer(&format!(
        "{directory}/broad_recovery_predictions_{label}.csv"
    ))?;
    writeln!(w, "gene,day,penalty,model,observed,predicted,squared_error")?;
    if rows.len() < 12 {
        return Ok(format!("{label}: insufficient regions for forecasting.\n"));
    }
    for (day, _) in days.iter().enumerate() {
        let available: Vec<_> = rows.iter().filter(|r| r.outcomes[day].is_some()).collect();
        if available.len() < 12 {
            continue;
        }
        for (penalty_index, penalty) in [0.001, 0.01, 0.1].into_iter().enumerate() {
            for held in &available {
                for model in 0..3 {
                    let p = held.features.len() + usize::from(model > 0);
                    let x = |r: &Recovery| {
                        let mut v = r.features.clone();
                        if model > 0 {
                            v.push(if model == 1 { r.cost } else { r.cost.abs() });
                        }
                        v
                    };
                    let train: Vec<_> = available.iter().filter(|r| r.gene != held.gene).collect();
                    let center: Vec<_> = (0..p)
                        .map(|j| {
                            if j == 0 {
                                0.0
                            } else {
                                mean(&train.iter().map(|r| x(r)[j]).collect::<Vec<_>>())
                            }
                        })
                        .collect();
                    let scale: Vec<_> = (0..p)
                        .map(|j| {
                            if j == 0 {
                                1.0
                            } else {
                                mean(
                                    &train
                                        .iter()
                                        .map(|r| (x(r)[j] - center[j]).powi(2))
                                        .collect::<Vec<_>>(),
                                )
                                .sqrt()
                                .max(1e-8)
                            }
                        })
                        .collect();
                    let standardized = |r: &Recovery| {
                        x(r).iter()
                            .enumerate()
                            .map(|(j, v)| (v - center[j]) / scale[j])
                            .collect::<Vec<_>>()
                    };
                    let mut normal = NormalEquations::new(p);
                    for r in train {
                        normal.add_xy(&standardized(r), r.outcomes[day].unwrap(), 1.0);
                    }
                    let (b, _) = normal.fit(p, penalty)?;
                    let predicted = standardized(held)
                        .iter()
                        .zip(b)
                        .map(|(a, b)| a * b)
                        .sum::<f64>();
                    let observed = held.outcomes[day].unwrap();
                    errors
                        .entry((day, penalty_index, model))
                        .or_default()
                        .push((observed - predicted).powi(2));
                    writeln!(
                        w,
                        "{},{},{penalty},{model},{observed},{predicted},{}",
                        held.gene,
                        days[day],
                        (observed - predicted).powi(2)
                    )?;
                }
            }
        }
    }
    let mut summary = String::new();
    for (day, _) in days.iter().enumerate() {
        for (penalty_index, penalty) in [0.001, 0.01, 0.1].into_iter().enumerate() {
            if let Some(baseline) = errors.get(&(day, penalty_index, 0)) {
                let mse = mean(baseline);
                let signed = mean(&errors[&(day, penalty_index, 1)]);
                let absolute = mean(&errors[&(day, penalty_index, 2)]);
                summary.push_str(&format!(
                    "| {label} | {} | {penalty} | {} | {:.5} | {:.2}% | {:.2}% |\n",
                    days[day],
                    baseline.len(),
                    mse.sqrt(),
                    100.0 * (signed / mse - 1.0),
                    100.0 * (absolute / mse - 1.0)
                ));
            }
        }
    }
    Ok(summary)
}
fn human_recovery(root: &str, out: &str, scores: &BTreeMap<String, f64>) -> Result<String> {
    let root = format!("{root}/broad_validation");
    let files: Vec<_> = fs::read_dir(&root)?
        .filter_map(|p| p.ok())
        .map(|p| p.file_name().to_string_lossy().into_owned())
        .filter(|s| s.ends_with(".cov.gz"))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    if files.len() != 27 {
        return Err("expected 27 recovery RRBS files".into());
    }
    let mut regions = vec![];
    let mut unresolved = 0;
    for line in reader("data/derived/broad_regions.tsv")?.lines().skip(1) {
        let f = fields(&line?, '\t');
        if !scores.contains_key(&f[0]) {
            continue;
        }
        if f[7] != "unique" {
            unresolved += 1;
            continue;
        }
        regions.push(Region {
            gene: f[0].clone(),
            chrom: f[1].clone(),
            start: f[2].parse()?,
            end: f[3].parse()?,
            density: f[8].parse()?,
            counts: BTreeMap::new(),
        });
    }
    let mut by_chrom: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, r) in regions.iter().enumerate() {
        by_chrom.entry(r.chrom.clone()).or_default().push(i);
    }
    for indices in by_chrom.values_mut() {
        indices.sort_by_key(|i| regions[*i].start);
    }
    for (sample, file) in files.iter().enumerate() {
        eprintln!("Regional recovery {}/{}: {file}", sample + 1, files.len());
        for line in reader(&format!("{root}/{file}"))?.lines() {
            let line = line?;
            let mut parts = line.split('\t');
            let chrom = parts.next().ok_or("RRBS chromosome missing")?;
            let Some(indices) = by_chrom.get(chrom) else {
                continue;
            };
            let position: u64 = parts.next().ok_or("RRBS position missing")?.parse()?;
            // Author-processed files are BED chr/start/end/total/methylated,
            // not the six-column Bismark format. Coordinates are already 0-based.
            let candidates = indices.partition_point(|i| regions[*i].start <= position);
            for index in &indices[..candidates] {
                if position >= regions[*index].end {
                    continue;
                }
                let mut tail = parts.clone();
                tail.next();
                let n: f64 = tail.next().ok_or("RRBS total count missing")?.parse()?;
                let m: f64 = tail
                    .next()
                    .ok_or("RRBS methylated count missing")?
                    .parse()?;
                if m < 0.0 || n < m {
                    return Err("invalid RRBS counts".into());
                }
                let counts = regions[*index]
                    .counts
                    .entry(position)
                    .or_insert_with(|| vec![(0.0, 0.0); files.len()]);
                counts[sample].0 += m;
                counts[sample].1 += n;
            }
        }
    }
    let find = |tag: &str| -> Result<usize> {
        files
            .iter()
            .position(|s| s.contains(tag))
            .ok_or_else(|| format!("missing recovery sample {tag}").into())
    };
    let baseline = [find("HCT116_exp_1_")?, find("HCT116_exp_2_")?];
    let early = find("HCT116_aza_d3_")?;
    let time_indices: Vec<_> = [6, 13, 16, 22, 28, 40]
        .iter()
        .map(|day| find(&format!("HCT116_aza_d{day}_")))
        .collect::<Result<_>>()?;
    let mut raw = writer("data/derived/broad_region_counts.csv")?;
    writeln!(raw, "gene,chrom,position_1based,sample,methylated,total")?;
    for r in &regions {
        for (site, counts) in &r.counts {
            for (sample, (m, n)) in counts.iter().enumerate() {
                writeln!(
                    raw,
                    "{},{},{},{},{m},{n}",
                    r.gene,
                    r.chrom,
                    site + 1,
                    files[sample]
                )?;
            }
        }
    }
    let mut w = writer(&format!("{out}/broad_recovery_regions.csv"))?;
    writeln!(w,"coverage,min_baseline,min_loss,gene,early_shared_cpgs,baseline,day3,loss,eligible,day28_fraction,day40_fraction,day40_vehicle_fraction,dnmt3bko_day40_fraction,dnmt1ko_day40_fraction")?;
    let gates = [
        (10.0, 0.5, 0.15, "primary"),
        (5.0, 0.5, 0.15, "coverage5"),
        (20.0, 0.5, 0.15, "coverage20"),
        (10.0, 0.3, 0.15, "baseline03"),
        (10.0, 0.7, 0.15, "baseline07"),
        (10.0, 0.5, 0.1, "loss010"),
        (10.0, 0.5, 0.2, "loss020"),
    ];
    let mut report = String::new();
    let mut prediction_summary=String::from("\n| Gate / score | Day | Penalty | Regions | Context RMSE | Signed cost MSE change | Absolute cost MSE change |\n|---|---:|---:|---:|---:|---:|---:|\n");
    for (coverage, min_baseline, min_loss, label) in gates {
        let mut rows = vec![];
        let mut speed_rows = vec![];
        let mut diagnostics = [vec![], vec![], vec![], vec![], vec![]];
        let mut no_early = 0;
        let mut losses = vec![];
        for r in &regions {
            let sites: Vec<_> = r
                .counts
                .iter()
                .filter(|(_, v)| {
                    v[baseline[0]].1 >= coverage
                        && v[baseline[1]].1 >= coverage
                        && v[early].1 >= coverage
                })
                .map(|(p, _)| *p)
                .collect();
            let bases = baseline
                .iter()
                .filter_map(|sample| region_beta(r, &sites, *sample, coverage))
                .collect::<Vec<_>>();
            let Some(d3) = region_beta(r, &sites, early, coverage) else {
                no_early += 1;
                continue;
            };
            if bases.len() != 2 {
                no_early += 1;
                continue;
            }
            let base = mean(&bases);
            let loss = base - d3;
            losses.push(loss);
            let eligible = base >= min_baseline && loss >= min_loss;
            let fraction =
                |sample| region_beta(r, &sites, sample, coverage).map(|x| (x - d3) / loss);
            let later = [fraction(time_indices[4]), fraction(time_indices[5])];
            let control = fraction(find("HCT116_cont_d40_")?);
            // KO fractions use their own early baseline/loss on the same fixed WT sites.
            let knockout = |prefix: &str| -> Result<Option<f64>> {
                let b = region_beta(r, &sites, find(&format!("{prefix}_RRBS"))?, coverage);
                let e = region_beta(r, &sites, find(&format!("{prefix}_aza_d3_"))?, coverage);
                let l = region_beta(r, &sites, find(&format!("{prefix}_aza_d40_"))?, coverage);
                Ok(match (b, e, l) {
                    (Some(b), Some(e), Some(l)) if b - e >= min_loss && b >= min_baseline => {
                        Some((l - e) / (b - e))
                    }
                    _ => None,
                })
            };
            let ko3 = knockout("HCT116_3BKO")?;
            let ko1 = knockout("HCT116_1KO")?;
            let display = |x: Option<f64>| x.map_or_else(|| "NA".to_string(), |v| v.to_string());
            writeln!(w,"{coverage},{min_baseline},{min_loss},{},{},{base},{d3},{loss},{eligible},{},{},{},{},{}",r.gene,sites.len(),display(later[0]),display(later[1]),display(control),display(ko3),display(ko1))?;
            if eligible {
                for (i, value) in [later[0], later[1], control, ko3, ko1].iter().enumerate() {
                    if let Some(value) = value {
                        diagnostics[i].push(*value);
                    }
                }
                let e6 = fraction(time_indices[0]);
                let e13 = fraction(time_indices[1]);
                let e22 = fraction(time_indices[3]);
                speed_rows.push(Recovery {
                    gene: r.gene.clone(),
                    features: vec![1.0, base, loss, r.density],
                    cost: scores[&r.gene],
                    outcomes: [e6, e13],
                });
                if let (Some(e6), Some(e13), Some(e22)) = (e6, e13, e22) {
                    rows.push(Recovery {
                        gene: r.gene.clone(),
                        features: vec![1.0, base, loss, r.density, e6, e13, e22],
                        cost: scores[&r.gene],
                        outcomes: later,
                    });
                }
            }
        }
        prediction_summary.push_str(&recovery_predictions(&rows, out, label, [28, 40])?);
        if label == "primary" {
            prediction_summary.push_str(&recovery_predictions(
                &speed_rows,
                out,
                "exploratory_speed",
                [6, 13],
            )?);
            for (i, name) in [
                "WT drug day28",
                "WT drug day40",
                "WT vehicle day40",
                "DNMT3B KO day40",
                "DNMT1 KO day40",
            ]
            .iter()
            .enumerate()
            {
                report.push_str(&format!("- {name}: {} eligible regions; median recovery fraction {:.4}, quartiles {:.4} / {:.4}.\n",diagnostics[i].len(),quantile(diagnostics[i].clone(),0.5),quantile(diagnostics[i].clone(),0.25),quantile(diagnostics[i].clone(),0.75)));
            }
            let simpler: Vec<_> = rows
                .iter()
                .cloned()
                .map(|mut r| {
                    r.features.truncate(4);
                    r
                })
                .collect();
            prediction_summary.push_str(&recovery_predictions(
                &simpler,
                out,
                "exploratory_baseline_context",
                [28, 40],
            )?);
            for (name, cell, column, sign) in [
                ("exploratory_dld1", "DLD1", 3, 1.0),
                ("exploratory_vp64", "HCT116", 6, -1.0),
                ("exploratory_inactive", "HCT116", 8, -1.0),
            ] {
                let mut external = BTreeMap::new();
                for line in reader(&format!("{out}/broad_functional_scores.csv"))?
                    .lines()
                    .skip(1)
                {
                    let f = fields(&line?, ',');
                    if f[0] == cell {
                        external.insert(f[1].clone(), sign * f[column].parse::<f64>()?);
                    }
                }
                let sensitivity: Vec<_> = rows
                    .iter()
                    .cloned()
                    .map(|mut r| {
                        r.cost = external[&r.gene];
                        r
                    })
                    .collect();
                prediction_summary.push_str(&recovery_predictions(
                    &sensitivity,
                    out,
                    name,
                    [28, 40],
                )?);
            }
        }
        report.push_str(&format!("- {label}: {} sequence-joined regions, {no_early} fail early coverage; {} have sufficient perturbation and training coverage, {} / {} have day-28 / day-40 outcomes. Median early methylation loss {:.4}.\n",regions.len(),rows.len(),rows.iter().filter(|r|r.outcomes[0].is_some()).count(),rows.iter().filter(|r|r.outcomes[1].is_some()).count(),quantile(losses,0.5)));
    }
    Ok(format!("## Independent local function versus human recovery\n\n56 screen targets; {unresolved} fail the conservative two-guide/unique-assembly gate. Mean log viability is normalized against negative controls within construct and plate, then averaged equally across available batches. Each target belongs to one supplied batch with eight wells; the three batch IDs distribute targets rather than replicate each target. No biological significance test is inferred from wells. Signed cost is inactive minus active; absolute cost is an explicit sensitivity. VP64 and DLD1 are reported separately.\n\n{report}\nPositive MSE changes mean worse held-out prediction. Primary scores and gate sensitivities were frozen before outcomes. Early recovery-speed, baseline-only context, DLD1 score, VP64 and inactive-construct comparisons below were added as exploratory robustness checks after the first run; none replaces the primary test.\n{prediction_summary}\nA coverage-qualified late regional mean requires at least three CpGs and 80% of the frozen early CpG set, each with the specified read coverage. CpGs receive equal weight. Recovery fractions may overshoot; none are clipped. Late vehicle uses the same WT reference and diagnoses untreated culture drift. KO fractions use genotype-specific early reference on the same WT-defined sites. There is one recovery culture per genotype/time: CpGs are not biological replicates.\n\nPrediction CSVs contain whole-region-excluded ridge predictions at three fixed penalties, with training-only feature scaling. Other regions' outcomes at the evaluated day calibrate the regression: this is cross-region generalization, not a simultaneous global future-time holdout or a new-condition transfer test; models 0/1/2 are early recovery/context, plus signed cost, plus absolute cost. Empty prediction tables mean the frozen >=12-region feasibility threshold failed. H3K36me3, expression, division and clone measurements are unavailable in this regional join, so any surviving association would remain mechanistically confounded. Active-versus-inactive editing also lacks a measured per-target editing-efficiency control.\n"))
}

// m' = a exp(-clearance*t)(1-m) - d*m. All rate parameters are /day.
// A midpoint transition with frozen rates is a positivity-preserving integration.
fn kinetic_curve(initial: f64, a: f64, d: f64, clearance: f64) -> [f64; 7] {
    let mut values = [initial; 7];
    let times: [f64; 7] = [0.0, 4.0, 8.0, 10.0, 13.0, 17.0, 29.0];
    let mut t = 0.0;
    let mut m = initial;
    for i in 1..times.len() {
        while t < times[i] - 1e-9 {
            let dt: f64 = (times[i] - t).min(0.1);
            let on = a * (-clearance * (t + dt / 2.0)).exp();
            let rate = on + d;
            if rate > 0.0 {
                let equilibrium = on / rate;
                m = equilibrium + (m - equilibrium) * (-rate * dt).exp();
            }
            t += dt;
        }
        values[i] = m;
    }
    values
}
#[derive(Clone)]
struct Amplicon {
    chrom: String,
    position: u64,
    tag: String,
    methylation: [[f64; 7]; 3],
    mock: [[Option<f64>; 6]; 3],
}
fn load_amplicons(root: &str, genotype: &str) -> Result<Vec<Amplicon>> {
    let mut lines = reader(&format!(
        "{root}/broad_validation/kinetics_{genotype}.txt.gz"
    ))?
    .lines();
    let header = fields(&lines.next().ok_or("empty kinetics")??, '\t');
    let locate = |name: &str| {
        header
            .iter()
            .position(|x| x == name)
            .ok_or_else(|| format!("missing kinetics column {name}"))
    };
    let prefix = if genotype == "WT" { "WT" } else { "TTKO" };
    let days = [0, 4, 8, 10, 13, 17, 29];
    let mut rows = vec![];
    let tag = locate("tag")?;
    for line in lines {
        let f = fields(&line?, '\t');
        let mut methylation = [[0.0; 7]; 3];
        let mut mock = [[None; 6]; 3];
        let mut eligible = true;
        for replicate in 1..=3 {
            for (time, day) in days.iter().enumerate() {
                let stem = if time == 0 {
                    format!("{prefix}_{replicate}_d0")
                } else {
                    format!("{prefix}_cre{replicate}_d{day}")
                };
                let n = numeric(&f[locate(&format!("{stem}.coverage"))?]);
                let m = numeric(&f[locate(&format!("{stem}.methRatio"))?]);
                match (n, m) {
                    (Some(n), Some(m)) if n >= 10.0 && (0.0..=1.0).contains(&m) => {
                        methylation[replicate - 1][time] = m
                    }
                    _ => eligible = false,
                }
                if time > 0 {
                    let stem = format!("{prefix}_m{replicate}_d{day}");
                    let n = numeric(&f[locate(&format!("{stem}.coverage"))?]);
                    let m = numeric(&f[locate(&format!("{stem}.methRatio"))?]);
                    if n.is_some_and(|n| n >= 10.0) {
                        mock[replicate - 1][time - 1] = m.filter(|m| (0.0..=1.0).contains(m));
                    }
                }
            }
        }
        if eligible {
            rows.push(Amplicon {
                chrom: f[0].clone(),
                position: f[1].parse()?,
                tag: f[tag].clone(),
                methylation,
                mock,
            });
        }
    }
    Ok(rows)
}
// Empirical two-parameter loss curve anchored at day 4, after the early
// transduction transient. The fitted amplitude is independent of held-out culture.
fn delayed_loss(observed: &[[f64; 7]]) -> [f64; 7] {
    let times: [f64; 7] = [0.0, 4.0, 8.0, 10.0, 13.0, 17.0, 29.0];
    let mut best = (f64::INFINITY, 0.0, 0.0);
    for d in std::iter::once(0.0).chain((0..=40).map(|i| 10.0_f64.powf(-4.0 + f64::from(i) * 0.1)))
    {
        let x: [f64; 7] = std::array::from_fn(|t| (-d * (times[t] - 4.0)).exp());
        let numerator = observed
            .iter()
            .map(|y| (1..=4).map(|t| y[t] * x[t]).sum::<f64>())
            .sum::<f64>();
        let denominator = observed.len() as f64 * (1..=4).map(|t| x[t] * x[t]).sum::<f64>();
        let amplitude = (numerator / denominator).clamp(0.0, 1.0);
        let error = observed
            .iter()
            .map(|y| {
                (1..=4)
                    .map(|t| (amplitude * x[t] - y[t]).powi(2))
                    .sum::<f64>()
            })
            .sum::<f64>();
        if error < best.0 {
            best = (error, amplitude, d);
        }
    }
    std::array::from_fn(|t| (best.1 * (-best.2 * (times[t] - 4.0)).exp()).clamp(0.0, 1.0))
}

type KineticCurve = (f64, f64, [f64; 7], [f64; 7]);
type KineticGrid = (&'static str, Vec<KineticCurve>);

fn genome_kinetics(root: &str, out: &str, grids: &[KineticGrid]) -> Result<String> {
    let mut errors: BTreeMap<(&str, String), Vec<f64>> = BTreeMap::new();
    let mut predictions = writer("data/derived/broad_genome_predictions.csv")?;
    writeln!(
        predictions,
        "chrom,position,day,model,observed,predicted,squared_error"
    )?;
    let mut sampled = 0;
    let mut retained = 0;
    for (index, line) in reader(&format!("{root}/broad_validation/kinetics_genome.txt.gz"))?
        .lines()
        .skip(1)
        .enumerate()
    {
        let line = line?;
        // Fixed input-index thinning, before coverage, labels or methylation.
        if index % 256 != 0 {
            continue;
        }
        sampled += 1;
        let f = fields(&line, '\t');
        if f.len() != 24 {
            return Err("unexpected genome kinetics format".into());
        }
        let mut observed = [0.0; 7];
        let mut eligible = true;
        for t in 0..7 {
            match (numeric(&f[8 + t]), numeric(&f[16 + t])) {
                (Some(m), Some(n)) if n >= 50.0 && (0.0..=100.0).contains(&m) => {
                    observed[t] = m / 100.0
                }
                _ => eligible = false,
            }
        }
        if !eligible {
            continue;
        }
        retained += 1;
        let mut models = vec![
            ("initial_constant", [observed[0]; 7]),
            ("last_training", [observed[4]; 7]),
            ("delayed_empirical_loss", delayed_loss(&[observed])),
        ];
        for (name, curves) in grids {
            let best = curves
                .iter()
                .min_by(|a, b| {
                    let loss = |(_, _, zero, one): &KineticCurve| {
                        (1..=4)
                            .map(|t| {
                                let p = zero[t] + observed[0] * (one[t] - zero[t]);
                                (p - observed[t]).powi(2)
                            })
                            .sum::<f64>()
                    };
                    loss(a).total_cmp(&loss(b))
                })
                .ok_or("empty rate grid")?;
            let curve = std::array::from_fn(|t| best.2[t] + observed[0] * (best.3[t] - best.2[t]));
            models.push((*name, curve));
        }
        for (model, curve) in models {
            for t in [5, 6] {
                let error = (curve[t] - observed[t]).powi(2);
                errors.entry((model, f[0].clone())).or_default().push(error);
                writeln!(
                    predictions,
                    "{},{},{},{model},{},{},{error}",
                    f[0],
                    f[1],
                    [0, 4, 8, 10, 13, 17, 29][t],
                    observed[t],
                    curve[t]
                )?;
            }
        }
    }
    let mut blocks = writer(&format!("{out}/broad_genome_summary.csv"))?;
    writeln!(blocks, "chrom,model,forecasts,mse")?;
    let mut totals: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
    for ((model, chrom), v) in errors {
        writeln!(blocks, "{chrom},{model},{},{}", v.len(), mean(&v))?;
        totals.entry(model).or_default().extend(v);
    }
    let mut report=format!("\n## Exploratory genome-wide enzyme benchmark\n\nFixed every-256th input-row sampling retains {retained} / {sampled} sampled CpGs with >=50 reads at all seven times. The processed mouse TTKO table aggregates cultures, so this is a temporal holdout, not an additional biological replication study. Methylation values are percentages converted to fractions. Published fitted rates, TET labels and confidence classifications are ignored. The processed source already requires coverage across all times and reflects the authors' data filtering; this is not an untouched raw-data evaluation. Training uses days 4/8/10/13; evaluation uses days 17/29. This broader benchmark was added after the amplicon results and is exploratory.\n\n| Model | Future RMSE |\n|---|---:|\n");
    for (name, v) in totals {
        report.push_str(&format!("| {name} | {:.5} |\n", mean(&v).sqrt()));
    }
    Ok(report)
}

fn mouse_kinetics(root: &str, out: &str) -> Result<String> {
    let mut forecasts = writer("data/derived/broad_kinetic_predictions.csv")?;
    writeln!(forecasts,"genotype,chrom,position,amplicon,held_culture,day,model,a_per_day,d_per_day,observed,predicted,squared_error,boundary")?;
    let mut mocks = writer(&format!("{out}/broad_mock_controls.csv"))?;
    writeln!(
        mocks,
        "genotype,amplicon,culture,day,cpgs,mean_change_from_baseline"
    )?;
    let rates: Vec<_> = std::iter::once(0.0)
        .chain((0..=40).map(|i| 10.0_f64.powf(-4.0 + f64::from(i) * 0.1)))
        .collect();
    let models = [
        ("constant_ctmc", 0.0),
        ("clearance025", 0.25),
        ("clearance050", 0.5),
        ("clearance100", 1.0),
    ];
    let grids: Vec<KineticGrid> = models
        .iter()
        .map(|(name, clearance)| {
            let mut curves = Vec::new();
            for a in &rates {
                for d in &rates {
                    curves.push((
                        *a,
                        *d,
                        kinetic_curve(0.0, *a, *d, *clearance),
                        kinetic_curve(1.0, *a, *d, *clearance),
                    ));
                }
            }
            (*name, curves)
        })
        .collect();
    let mut report = String::new();
    for genotype in ["WT", "TTKO"] {
        let rows = load_amplicons(root, genotype)?;
        let tags: BTreeSet<_> = rows.iter().map(|r| r.tag.clone()).collect();
        let mut errors: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
        let mut blocks: BTreeMap<(&str, String), Vec<f64>> = BTreeMap::new();
        let mut bounds: BTreeMap<&str, usize> = BTreeMap::new();
        for (site, row) in rows.iter().enumerate() {
            if site % 100 == 0 {
                eprintln!("{genotype}: kinetic CpG {site}/{}", rows.len());
            }
            for held in 0..3 {
                let train: Vec<_> = (0..3).filter(|r| *r != held).collect();
                let observed = &row.methylation[held];
                let last = mean(
                    &train
                        .iter()
                        .map(|r| row.methylation[*r][4])
                        .collect::<Vec<_>>(),
                );
                let mut predictions = vec![
                    ("initial_constant", [observed[0]; 7], 0.0, 0.0, false),
                    ("last_training", [last; 7], 0.0, 0.0, false),
                ];
                predictions.push((
                    "delayed_empirical_loss",
                    delayed_loss(
                        &train
                            .iter()
                            .map(|rep| row.methylation[*rep])
                            .collect::<Vec<_>>(),
                    ),
                    0.0,
                    0.0,
                    false,
                ));
                for (model, curves) in &grids {
                    let mut best = (f64::INFINITY, 0.0, 0.0);
                    for (a, d, zero, one) in curves {
                        let error = train
                            .iter()
                            .map(|rep| {
                                (1..=4)
                                    .map(|t| {
                                        let predicted =
                                            zero[t] + row.methylation[*rep][0] * (one[t] - zero[t]);
                                        (predicted - row.methylation[*rep][t]).powi(2)
                                    })
                                    .sum::<f64>()
                            })
                            .sum::<f64>();
                        if error < best.0 {
                            best = (error, *a, *d);
                        }
                    }
                    let clearance = models.iter().find(|(name, _)| name == model).unwrap().1;
                    let boundary = best.1 == 1.0 || best.2 == 1.0;
                    predictions.push((
                        *model,
                        kinetic_curve(observed[0], best.1, best.2, clearance),
                        best.1,
                        best.2,
                        boundary,
                    ));
                }
                for (model, curve, a, d, boundary) in predictions {
                    for t in [5, 6] {
                        let error = (curve[t] - observed[t]).powi(2);
                        errors.entry(model).or_default().push(error);
                        blocks
                            .entry((model, row.tag.clone()))
                            .or_default()
                            .push(error);
                        *bounds.entry(model).or_default() += usize::from(boundary);
                        writeln!(
                            forecasts,
                            "{genotype},{},{},{},{},{},{model},{a},{d},{},{},{error},{boundary}",
                            row.chrom,
                            row.position,
                            row.tag,
                            held + 1,
                            [0, 4, 8, 10, 13, 17, 29][t],
                            observed[t],
                            curve[t]
                        )?;
                    }
                }
            }
        }
        report.push_str(&format!("\n### {genotype}: {} CpGs in {} amplicons\n\n| Model | Future RMSE | Equal-amplicon MSE | Upper-rate-bound forecasts |\n|---|---:|---:|---:|\n",rows.len(),tags.len()));
        for (model, error) in &errors {
            let equal = blocks
                .iter()
                .filter(|((m, _), _)| m == model)
                .map(|(_, v)| mean(v))
                .collect::<Vec<_>>();
            report.push_str(&format!(
                "| {model} | {:.5} | {:.6} | {} / {} |\n",
                mean(error).sqrt(),
                mean(&equal),
                bounds[model],
                error.len()
            ));
        }
        for tag in tags {
            for replicate in 0..3 {
                for t in 0..6 {
                    let changes: Vec<_> = rows
                        .iter()
                        .filter(|r| r.tag == tag)
                        .filter_map(|r| {
                            r.mock[replicate][t].map(|m| m - r.methylation[replicate][0])
                        })
                        .collect();
                    if !changes.is_empty() {
                        writeln!(
                            mocks,
                            "{genotype},{tag},{},{},{},{}",
                            replicate + 1,
                            [4, 8, 10, 13, 17, 29][t],
                            changes.len(),
                            mean(&changes)
                        )?;
                    }
                }
            }
        }
    }
    report.push_str(&genome_kinetics(root, out, &grids)?);
    Ok(format!("## Replicated mouse enzyme-deletion forecasts\n\nRates fit only Cre days 4/8/10/13 in two cultures; predictions exclude the third culture and forecast days 17/29 from its measured day-0 baseline. Three rotations share training data and are not three independent experiments. Equal CpG/time/culture squared error is the fitting criterion; read counts are coverage filters, not independent biological samples. No published confidence label or fitted rate selects our subset. The source table itself was filtered by its authors using all time points and excludes one high-variability amplicon; this limits the independence of the retrospective validation. An exploratory two-parameter empirical comparator fits an amplitude and exponential loss rate from day 4 onward, without using the held-out baseline. It was added after the initial amplicon results and does not identify enzyme clearance. WT has TET activity; TTKO lacks TET enzymes. Both undergo de-novo DNMT deletion. Constant CTMC and decaying de-novo activity have two site-specific rates; fixed clearance 0.5/day is the published enzyme-clearance mechanism, 0.25/1.0 are prespecified sensitivities. Grid rates include zero and log10 -4..0 in 0.1 steps, maximum 1/day.\n{report}\nMock controls are summarized separately by amplicon/culture/time. Clearance fits test an already published intervention mechanism, not functional protection. Amplicons are correlated CpG blocks and were deliberately selected by the original authors; hundreds of CpGs do not constitute hundreds of independent biological replicates. No uncertainty interval assumes otherwise.\n"))
}

pub fn run(root: &str, out: &str) -> Result<()> {
    fs::create_dir_all(out)?;
    let scores = screen(out)?;
    let human = human_recovery(root, out, &scores)?;
    let mouse = mouse_kinetics(root, out)?;
    let mut control_groups: BTreeMap<(String, String), Vec<f64>> = BTreeMap::new();
    for line in reader(&format!("{out}/broad_mock_controls.csv"))?
        .lines()
        .skip(1)
    {
        let f = fields(&line?, ',');
        control_groups
            .entry((f[0].clone(), f[3].clone()))
            .or_default()
            .push(f[5].parse()?);
    }
    let mut controls = String::from("\n## Mock-control diagnostic\n\nEqual amplicon/culture mean changes from day 0, in methylation fractions:\n\n| Genotype | Day | Mean change |\n|---|---:|---:|\n");
    for ((genotype, day), values) in control_groups {
        controls.push_str(&format!("| {genotype} | {day} | {:.4} |\n", mean(&values)));
    }
    controls.push_str("\nWT mock methylation rises during the early transduction period. Primary predictions use uncorrected Cre fractions; clearance and empirical transient fits can absorb this background. The near-competitive empirical loss curve and mock drift prevent treating forecast gains as unique evidence for the exact clearance equation or as isolated TET effects. TTKO controls drift much less. No genotype-transfer rule or functional-restoration law is established here.\n");
    let mut report = writer(&format!("{out}/broad_validation.md"))?;
    writeln!(report,"# Public-data stopping tests (2026-10-04)\n\nThe [stopping assessment](../../notes/experiment_brief.md#computational-stopping-assessment-2026-10-04) distinguishes negative functional predictions from generic kinetic plausibility. Protocol frozen in notes/research_plan.md before new recovery outcomes. Inputs verified against broad_data_manifest.csv and broad_reference_manifest.csv. Python handles acquisition, verification, XLSX and reference coordinate formats only; all estimates and predictions are Rust.\n\n{human}\n{mouse}\n{controls}")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn clearance_solver_preserves_bounds_and_constant_ctmc() {
        for initial in [0.0, 0.3, 1.0] {
            let result = kinetic_curve(initial, 0.2, 0.1, 0.0);
            let exact = 2.0 / 3.0 + (initial - 2.0 / 3.0) * (-0.3_f64 * 29.0).exp();
            assert!((result[6] - exact).abs() < 1e-12);
            for clearance in [0.25, 0.5, 1.0] {
                assert!(kinetic_curve(initial, 1.0, 0.1, clearance)
                    .iter()
                    .all(|p| (0.0..=1.0).contains(p)));
            }
        }
        assert_eq!(kinetic_curve(0.4, 0.0, 0.0, 0.5), [0.4; 7]);
        let initial = 0.3;
        let exact = 1.0 - (1.0 - initial) * (-(1.0 - (-29.0_f64).exp())).exp();
        assert!((kinetic_curve(initial, 1.0, 0.0, 1.0)[6] - exact).abs() < 2e-4);
    }
    #[test]
    fn empirical_loss_forecasts_known_decay_from_training_only() {
        let times: [f64; 7] = [0.0, 4.0, 8.0, 10.0, 13.0, 17.0, 29.0];
        let mut observed = std::array::from_fn(|t| 0.8 * (-0.1 * (times[t] - 4.0)).exp());
        let expected = observed;
        observed[5] = 0.99;
        observed[6] = 0.01;
        let prediction = delayed_loss(&[observed]);
        assert!((prediction[5] - expected[5]).abs() < 1e-12);
        assert!((prediction[6] - expected[6]).abs() < 1e-12);
    }
    #[test]
    fn missing_late_cpgs_do_not_become_zero() {
        let mut region = Region::default();
        for site in 0..5 {
            region.counts.insert(site, vec![(8.0, 10.0), (0.0, 0.0)]);
        }
        let sites: Vec<_> = (0..5).collect();
        assert_eq!(region_beta(&region, &sites, 1, 10.0), None);
        assert_eq!(region_beta(&region, &sites, 0, 10.0), Some(0.8));
    }
}
