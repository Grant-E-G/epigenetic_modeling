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
fn export_audit_design(
    rows: &[Recovery],
    out: &str,
    label: &str,
    output: &mut BufWriter<File>,
) -> Result<()> {
    let mut batches = BTreeMap::new();
    for line in reader(&format!("{out}/broad_screen_batches.csv"))?
        .lines()
        .skip(1)
    {
        let f = fields(&line?, ',');
        if f[0] == "HCT116" {
            batches.insert(f[1].clone(), f[2].clone());
        }
    }
    for r in rows {
        let mut values = vec![
            label.to_string(),
            r.gene.clone(),
            batches[&r.gene].clone(),
            r.cost.to_string(),
        ];
        values.extend(
            r.outcomes
                .iter()
                .map(|v| v.map_or_else(|| "NA".into(), |v| v.to_string())),
        );
        values.extend((0..7).map(|i| {
            r.features
                .get(i)
                .map_or_else(|| "NA".into(), ToString::to_string)
        }));
        writeln!(output, "{}", values.join(","))?;
    }
    Ok(())
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
    let mut audit = writer("data/derived/broad_power_design.csv")?;
    writeln!(
        audit,
        "analysis,gene,batch,functional_cost,y0,y1,x0,x1,x2,x3,x4,x5,x6"
    )?;
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
            export_audit_design(&rows, out, "primary", &mut audit)?;
            export_audit_design(&speed_rows, out, "speed", &mut audit)?;
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

// Linear weights reproduce the frozen ridge pipeline exactly while avoiding
// repeated matrix factorization during thousands of design-conditional audits.
struct AuditFold {
    held: usize,
    train: Vec<usize>,
    x: Vec<Vec<f64>>,
    held_x: Vec<f64>,
    inverse: Vec<Vec<f64>>,
    baseline: Vec<f64>,
    penalty: f64,
}
fn audit_standardize(x: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let p = x[0].len();
    let center: Vec<_> = (0..p)
        .map(|j| {
            if j == 0 {
                0.0
            } else {
                mean(&x.iter().map(|r| r[j]).collect::<Vec<_>>())
            }
        })
        .collect();
    let scale: Vec<_> = (0..p)
        .map(|j| {
            if j == 0 {
                1.0
            } else {
                mean(
                    &x.iter()
                        .map(|r| (r[j] - center[j]).powi(2))
                        .collect::<Vec<_>>(),
                )
                .sqrt()
                .max(1e-8)
            }
        })
        .collect();
    let z = x
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .map(|(j, v)| (v - center[j]) / scale[j])
                .collect()
        })
        .collect();
    (center, scale, z)
}
fn audit_inverse(x: &[Vec<f64>], penalty: f64) -> Result<Vec<Vec<f64>>> {
    let p = x[0].len();
    let mut normal = NormalEquations::new(p);
    for r in x {
        normal.add_xy(r, 0.0, 1.0);
    }
    let mut inverse = vec![vec![0.0; p]; p];
    for j in 0..p {
        normal.rhs.fill(0.0);
        normal.rhs[j] = 1.0;
        let (column, _) = normal.fit(p, penalty)?;
        for i in 0..p {
            inverse[i][j] = column[i];
        }
    }
    Ok(inverse)
}
fn audit_dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
fn audit_product(inverse: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    inverse.iter().map(|r| audit_dot(r, rhs)).collect()
}
impl AuditFold {
    fn new(features: &[Vec<f64>], held: usize, penalty: f64) -> Result<Self> {
        let train: Vec<_> = (0..features.len()).filter(|i| *i != held).collect();
        let raw: Vec<_> = train.iter().map(|i| features[*i].clone()).collect();
        let (center, scale, x) = audit_standardize(&raw);
        let held_x: Vec<_> = features[held]
            .iter()
            .enumerate()
            .map(|(j, v)| (v - center[j]) / scale[j])
            .collect();
        let inverse = audit_inverse(&x, penalty)?;
        let coefficient = audit_product(&inverse, &held_x);
        let baseline = x.iter().map(|r| audit_dot(r, &coefficient)).collect();
        Ok(Self {
            held,
            train,
            x,
            held_x,
            inverse,
            baseline,
            penalty,
        })
    }
    fn weights(&self, q: &[f64]) -> Vec<f64> {
        let mut weights = vec![0.0; q.len()];
        let values: Vec<_> = self.train.iter().map(|i| q[*i]).collect();
        let center = mean(&values);
        let scale = mean(
            &values
                .iter()
                .map(|q| (q - center).powi(2))
                .collect::<Vec<_>>(),
        )
        .sqrt()
        .max(1e-8);
        let z: Vec<_> = values.iter().map(|q| (q - center) / scale).collect();
        let rhs: Vec<_> = (0..self.held_x.len())
            .map(|j| self.x.iter().zip(&z).map(|(x, q)| x[j] * q).sum())
            .collect();
        let projection = audit_product(&self.inverse, &rhs);
        let denominator = audit_dot(&z, &z) + self.penalty * self.train.len() as f64
            - audit_dot(&rhs, &projection);
        let leverage =
            ((q[self.held] - center) / scale - audit_dot(&self.held_x, &projection)) / denominator;
        for (j, i) in self.train.iter().enumerate() {
            weights[*i] = self.baseline[j] + leverage * (z[j] - audit_dot(&self.x[j], &projection));
        }
        weights
    }
    fn baseline_weights(&self, n: usize) -> Vec<f64> {
        let mut w = vec![0.0; n];
        for (j, i) in self.train.iter().enumerate() {
            w[*i] = self.baseline[j];
        }
        w
    }
}
fn audit_gain(y: &[f64], base: &[Vec<f64>], model: &[Vec<f64>]) -> f64 {
    // Explicit loops keep observation vectors separate from scalar outcomes.
    let mse = |matrix: &[Vec<f64>]| {
        matrix
            .iter()
            .enumerate()
            .map(|(i, w)| (audit_dot(w, y) - y[i]).powi(2))
            .sum::<f64>()
    };
    1.0 - mse(model) / mse(base).max(1e-15)
}
fn audit_full_weights(x: &[Vec<f64>], penalty: f64) -> Result<Vec<Vec<f64>>> {
    let (_, _, z) = audit_standardize(x);
    let inverse = audit_inverse(&z, penalty)?;
    Ok(z.iter()
        .map(|r| {
            let coefficient = audit_product(&inverse, r);
            z.iter().map(|t| audit_dot(t, &coefficient)).collect()
        })
        .collect())
}
fn audit_fit_predict(
    rows: &[&Recovery],
    held: &Recovery,
    day: usize,
    penalty: f64,
    with_q: bool,
) -> Result<f64> {
    let feature = |r: &Recovery| {
        let mut x = r.features.clone();
        if with_q {
            x.push(r.cost);
        }
        x
    };
    let raw: Vec<_> = rows.iter().map(|r| feature(r)).collect();
    let (center, scale, x) = audit_standardize(&raw);
    let mut normal = NormalEquations::new(x[0].len());
    for (r, x) in rows.iter().zip(&x) {
        normal.add_xy(x, r.outcomes[day].unwrap(), 1.0);
    }
    let (b, _) = normal.fit(x[0].len(), penalty)?;
    let held = feature(held)
        .iter()
        .enumerate()
        .map(|(j, v)| (v - center[j]) / scale[j])
        .collect::<Vec<_>>();
    Ok(audit_dot(&held, &b))
}
fn audit_bootstrap(
    rows: &[Recovery],
    day: usize,
    penalty: f64,
    rng: &mut rand::rngs::StdRng,
) -> Result<f64> {
    use rand::Rng;
    let draw: Vec<_> = (0..rows.len())
        .map(|_| &rows[rng.gen_range(0..rows.len())])
        .collect();
    let mut errors = [0.0; 2];
    for held in &draw {
        let train: Vec<_> = draw
            .iter()
            .copied()
            .filter(|r| r.gene != held.gene)
            .collect();
        for (model, error) in errors.iter_mut().enumerate() {
            let p = audit_fit_predict(&train, held, day, penalty, model == 1)?;
            *error += (p - held.outcomes[day].unwrap()).powi(2);
        }
    }
    Ok(1.0 - errors[1] / errors[0].max(1e-15))
}
fn audit_coefficient_weights(features: &[Vec<f64>], q: &[f64], penalty: f64) -> Result<Vec<f64>> {
    let x: Vec<_> = features
        .iter()
        .zip(q)
        .map(|(x, q)| {
            let mut v = x.clone();
            v.push(*q);
            v
        })
        .collect();
    let (_, scale, z) = audit_standardize(&x);
    let inverse = audit_inverse(&z, penalty)?;
    Ok(z.iter()
        .map(|r| audit_dot(inverse.last().unwrap(), r) / scale.last().unwrap())
        .collect())
}
fn audit_bound(value: f64, row: &Recovery) -> f64 {
    let loss = row.features[2];
    let day3 = row.features[1] - loss;
    value.clamp(-day3 / loss, (1.0 - day3) / loss)
}
pub fn power_audit(out: &str) -> Result<()> {
    use rand::{seq::SliceRandom, Rng, SeedableRng};
    let mut designs: BTreeMap<String, Vec<(Recovery, String)>> = BTreeMap::new();
    for line in reader("data/derived/broad_power_design.csv")?
        .lines()
        .skip(1)
    {
        let f = fields(&line?, ',');
        let p = if f[0] == "primary" { 7 } else { 4 };
        let features = f[6..6 + p]
            .iter()
            .map(|x| x.parse())
            .collect::<std::result::Result<Vec<f64>, _>>()?;
        designs.entry(f[0].clone()).or_default().push((
            Recovery {
                gene: f[1].clone(),
                cost: f[3].parse()?,
                features,
                outcomes: [numeric(&f[4]), numeric(&f[5])],
            },
            f[2].clone(),
        ));
    }
    let mut summary = writer(&format!("{out}/functional_power_summary.csv"))?;
    writeln!(summary,"analysis,day,penalty,regions,observed_gain,global_permutation_p,stratified_permutation_p,global_null_gain_p05,global_null_gain_p50,global_null_gain_p95,stratified_null_gain_p05,stratified_null_gain_p50,stratified_null_gain_p95,bootstrap_gain_p025,bootstrap_gain_p975")?;
    let mut power = writer(&format!("{out}/functional_power_injection.csv"))?;
    writeln!(power,"analysis,day,penalty,regions,mode,target_oracle_gain,achieved_oracle_signal_fraction,beta_per_score_unit,effect_per_score_sd,noise_sd,coefficient_threshold_gain,any_gain_probability,calibrated_detection_probability,mc_standard_error,mean_realized_gain,bound_violation_fraction,replicates")?;
    let mut draws = writer("data/derived/functional_power_draws.csv")?;
    writeln!(
        draws,
        "analysis,day,penalty,procedure,mode,target,replicate,gain"
    )?;
    let mut report=String::from("# Functional-effect detectability audit\n\nPrimary audit: signed HCT116 score, day 40, ridge penalty 0.01. Other days/penalties are sensitivities. Positive gain means lower held-out MSE. The linear-weight implementation exactly reproduces the frozen training-only scaling, ridge penalties and whole-region exclusion. The design is regenerated from checksum-verified source inputs; neither injected outcomes nor resampling alter coverage/eligibility or missingness. These estimates condition on one observed recovery culture and noisy observed functional scores; they are not biological-culture power estimates.\n\n## Observed effect and empirical null\n\n| Analysis | Day | Penalty | Regions | Observed gain | Global permutation p | Within-batch permutation p | Bootstrap 95% range |\n|---|---:|---:|---:|---:|---:|---:|---:|\n");
    let mut primary_power=String::from("\n## Primary injection results\n\n| Generator | Target oracle reduction | Achieved signal fraction | Effect per score SD | Any MSE gain | Calibrated detection | MC SE |\n|---|---:|---:|---:|---:|---:|---:|\n");
    let mut diagnostics = String::from("\n## Fixed-design diagnostics\n\n| Analysis | Day | Regions | Score SD | Residual score SD | Score variance explained by context | Largest region share of residual score energy |\n|---|---:|---:|---:|---:|---:|---:|\n");
    for (analysis, design) in designs {
        for day in 0..2 {
            let available: Vec<_> = design
                .iter()
                .filter(|(r, _)| r.outcomes[day].is_some())
                .collect();
            let rows: Vec<_> = available.iter().map(|(r, _)| r.clone()).collect();
            let batches: Vec<_> = available.iter().map(|(_, b)| b.clone()).collect();
            let n = rows.len();
            let features: Vec<_> = rows.iter().map(|r| r.features.clone()).collect();
            let q: Vec<_> = rows.iter().map(|r| r.cost).collect();
            let y: Vec<_> = rows.iter().map(|r| r.outcomes[day].unwrap()).collect();
            let actual_day = if analysis == "primary" {
                [28, 40][day]
            } else {
                [6, 13][day]
            };
            // Context residualization defines an incremental available signal,
            // not a claim that functional effects change a biological rate.
            let projection = audit_full_weights(&features, 0.0)?;
            let q_res: Vec<_> = projection
                .iter()
                .enumerate()
                .map(|(i, w)| q[i] - audit_dot(w, &q))
                .collect();
            let qvar = mean(&q_res.iter().map(|q| q * q).collect::<Vec<_>>());
            let total_score_variance =
                mean(&q.iter().map(|v| (v - mean(&q)).powi(2)).collect::<Vec<_>>());
            let largest = q_res.iter().map(|q| q * q).fold(0.0_f64, f64::max) / (n as f64 * qvar);
            diagnostics.push_str(&format!(
                "| {analysis} | {actual_day} | {n} | {:.5} | {:.5} | {:.2}% | {:.2}% |\n",
                total_score_variance.sqrt(),
                qvar.sqrt(),
                100.0 * (1.0 - qvar / total_score_variance),
                100.0 * largest
            ));
            if qvar < 1e-15 {
                return Err("no residual functional-score variation".into());
            }
            for penalty in [0.001, 0.01, 0.1] {
                eprintln!("Power audit {analysis} day{actual_day} penalty{penalty}, n={n}");
                let seed = 20261004 + actual_day as u64 * 10000 + (penalty * 1000.0) as u64;
                let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
                let folds: Vec<_> = (0..n)
                    .map(|i| AuditFold::new(&features, i, penalty))
                    .collect::<Result<_>>()?;
                let base: Vec<_> = folds.iter().map(|f| f.baseline_weights(n)).collect();
                let model: Vec<_> = folds.iter().map(|f| f.weights(&q)).collect();
                let observed = audit_gain(&y, &base, &model);
                let mut permutations = [vec![], vec![]];
                let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
                for (i, b) in batches.iter().enumerate() {
                    groups.entry(b.clone()).or_default().push(i);
                }
                for (kind, distribution) in permutations.iter_mut().enumerate() {
                    for replicate in 0..999 {
                        let mut perm = q.clone();
                        if kind == 0 {
                            perm.shuffle(&mut rng);
                        } else {
                            for indices in groups.values() {
                                let mut shuffled: Vec<_> = indices.iter().map(|i| q[*i]).collect();
                                shuffled.shuffle(&mut rng);
                                for (i, value) in indices.iter().zip(shuffled) {
                                    perm[*i] = value;
                                }
                            }
                        }
                        let weights: Vec<_> = folds.iter().map(|f| f.weights(&perm)).collect();
                        let gain = audit_gain(&y, &base, &weights);
                        distribution.push(gain);
                        writeln!(draws,"{analysis},{actual_day},{penalty},permutation{kind},NA,NA,{replicate},{gain}")?;
                    }
                }
                let pvalue = |v: &[f64]| {
                    (1 + v.iter().filter(|g| **g >= observed).count()) as f64 / (v.len() + 1) as f64
                };
                let mut bootstrap = vec![];
                for replicate in 0..500 {
                    let gain = audit_bootstrap(&rows, day, penalty, &mut rng)?;
                    bootstrap.push(gain);
                    writeln!(
                        draws,
                        "{analysis},{actual_day},{penalty},bootstrap,NA,NA,{replicate},{gain}"
                    )?;
                }
                let interval = [
                    quantile(bootstrap.clone(), 0.025),
                    quantile(bootstrap, 0.975),
                ];
                writeln!(summary,"{analysis},{actual_day},{penalty},{n},{observed},{},{},{},{},{},{},{},{},{},{}",pvalue(&permutations[0]),pvalue(&permutations[1]),quantile(permutations[0].clone(),0.05),quantile(permutations[0].clone(),0.5),quantile(permutations[0].clone(),0.95),quantile(permutations[1].clone(),0.05),quantile(permutations[1].clone(),0.5),quantile(permutations[1].clone(),0.95),interval[0],interval[1])?;
                report.push_str(&format!("| {analysis} | {actual_day} | {penalty} | {n} | {:.2}% | {:.3} | {:.3} | {:.2}% .. {:.2}% |\n",100.0*observed,pvalue(&permutations[0]),pvalue(&permutations[1]),100.0*interval[0],100.0*interval[1]));
                let full = audit_full_weights(&features, penalty)?;
                let mu: Vec<_> = full.iter().map(|w| audit_dot(w, &y)).collect();
                let raw_res: Vec<_> = base
                    .iter()
                    .enumerate()
                    .map(|(i, w)| y[i] - audit_dot(w, &y))
                    .collect();
                let average = mean(&raw_res);
                let residual: Vec<_> = raw_res.iter().map(|e| e - average).collect();
                let noise = mean(&residual.iter().map(|e| e * e).collect::<Vec<_>>());
                let coefficients = audit_coefficient_weights(&features, &q, penalty)?;
                for bounded in [false, true] {
                    let mode = if bounded {
                        "physical_bounds"
                    } else {
                        "unbounded"
                    };
                    let clip =
                        |v: f64, i: usize| if bounded { audit_bound(v, &rows[i]) } else { v };
                    // Calibrate detection using a separate null bank; assessment
                    // repeats do not choose the rejection threshold.
                    let mut null_bank = vec![];
                    for _ in 0..2000 {
                        let sample: Vec<_> = (0..n)
                            .map(|i| {
                                clip(
                                    mu[i]
                                        + if rng.gen_bool(0.5) {
                                            residual[i]
                                        } else {
                                            -residual[i]
                                        },
                                    i,
                                )
                            })
                            .collect();
                        let g = audit_gain(&sample, &base, &model);
                        null_bank.push(if audit_dot(&coefficients, &sample) > 0.0 {
                            g
                        } else {
                            f64::NEG_INFINITY
                        });
                    }
                    // Quantile helper drops non-finite values; use a finite floor
                    // so negative-coefficient draws remain in the null bank.
                    for g in &mut null_bank {
                        if !g.is_finite() {
                            *g = -1e6;
                        }
                    }
                    let cutoff = quantile(null_bank, 0.95).max(0.0);
                    for target in [0.0, 0.03, 0.05, 0.1, 0.2] {
                        let beta = if target == 0.0 {
                            0.0
                        } else {
                            (target / (1.0 - target) * noise / qvar).sqrt()
                        };
                        let null_mean: Vec<_> = (0..n)
                            .map(|i| {
                                0.5 * (clip(mu[i] + residual[i], i) + clip(mu[i] - residual[i], i))
                            })
                            .collect();
                        let true_mean: Vec<_> = (0..n)
                            .map(|i| {
                                0.5 * (clip(mu[i] + residual[i] + beta * q_res[i], i)
                                    + clip(mu[i] - residual[i] + beta * q_res[i], i))
                            })
                            .collect();
                        let change: Vec<_> = true_mean
                            .iter()
                            .zip(null_mean)
                            .map(|(a, b)| a - b)
                            .collect();
                        let available_signal = mean(
                            &projection
                                .iter()
                                .enumerate()
                                .map(|(i, w)| (change[i] - audit_dot(w, &change)).powi(2))
                                .collect::<Vec<_>>(),
                        );
                        let actual_noise = mean(
                            &(0..n)
                                .map(|i| {
                                    0.25 * (clip(mu[i] + residual[i] + beta * q_res[i], i)
                                        - clip(mu[i] - residual[i] + beta * q_res[i], i))
                                    .powi(2)
                                })
                                .collect::<Vec<_>>(),
                        );
                        let achieved =
                            available_signal / (available_signal + actual_noise).max(1e-15);
                        let mut gains = vec![];
                        let mut any = 0;
                        let mut detected = 0;
                        let mut violations = 0;
                        for replicate in 0..2000 {
                            let sample: Vec<_> = (0..n)
                                .map(|i| {
                                    let v = mu[i]
                                        + if rng.gen_bool(0.5) {
                                            residual[i]
                                        } else {
                                            -residual[i]
                                        }
                                        + beta * q_res[i];
                                    if (audit_bound(v, &rows[i]) - v).abs() > 1e-12 {
                                        violations += 1;
                                    }
                                    clip(v, i)
                                })
                                .collect();
                            let gain = audit_gain(&sample, &base, &model);
                            any += usize::from(gain > 0.0);
                            detected += usize::from(
                                gain > cutoff && audit_dot(&coefficients, &sample) > 0.0,
                            );
                            gains.push(gain);
                            writeln!(draws,"{analysis},{actual_day},{penalty},injection,{mode},{target},{replicate},{gain}")?;
                        }
                        let probability = detected as f64 / 2000.0;
                        let se = (probability * (1.0 - probability) / 2000.0).sqrt();
                        let effect_sd = beta * variance(&q).sqrt();
                        writeln!(power,"{analysis},{actual_day},{penalty},{n},{mode},{target},{achieved},{beta},{effect_sd},{},{cutoff},{},{probability},{se},{},{},2000",noise.sqrt(),any as f64/2000.0,mean(&gains),violations as f64/(2000*n) as f64)?;
                        if analysis == "primary" && actual_day == 40 && penalty == 0.01 {
                            primary_power.push_str(&format!(
                                "| {mode} | {:.0}% | {:.2}% | {:.5} | {:.1}% | {:.1}% | {:.2}% |\n",
                                100.0 * target,
                                100.0 * achieved,
                                effect_sd,
                                100.0 * any as f64 / 2000.0,
                                100.0 * probability,
                                100.0 * se
                            ));
                        }
                    }
                }
            }
        }
    }
    report.push_str(&primary_power);
    report.push_str(&diagnostics);
    report.push_str("\n## Interpretation and limits\n\nAn injected beta acts on recovery fraction conditional on the real covariates: y*=context mean + beta*q_res + a signed region residual. It is not a DNMT rate change. Target oracle reductions 3/5/10/20% specify available fixed-design signal relative to the observed context leave-one-out residual variance. They do not assert that ridge achieves those reductions. q_res is the OLS projection residual of the actual observed score; score reliability is treated optimistically as perfect. Context means are ridge fits to this culture. Independent regional wild signs are an assumed error generator, not observed replicate biology. Read-depth-dependent residual magnitudes and actual missingness/filters remain fixed.\n\nPhysical-bound sensitivity clips only simulated methylation to [0,1], using each region's observed day3 and induced loss, then converts back to recovery fraction. It does not clip recovery fractions to [0,1]. The achieved signal fraction measures incremental simulated mean variation after context projection relative to simulated noise; clipping can alter the nominal oracle signal. This ratio is an operational signal/noise diagnostic, not a known true biological effect. Calibrated detection requires positive global fitted score coefficient and MSE gain above a separate zero-injection 95th percentile (floored at zero); zero-injection assessment reports realized false positives. MC standard errors describe simulations only.\n\nGlobal score permutations are not generally exchangeable conditional on context; within-screen-batch permutations address supplied batch grouping but not all biological confounding. Their p-values are empirical diagnostics, not causal significance. Bootstrap intervals resample complete region records and refit, excluding every copy of the held-out gene; they describe regional heterogeneity conditional on this culture. This resampling reduces the unique training-region count, so intervals are not a formally calibrated independent-culture confidence interval. All three penalties and four days are reported without selecting a successful variant.\n\nA low powered audit cannot permanently falsify functional restoration. Even high conditional detectability would only downgrade the tested measured-score/region/regression relationship, not all kinetic effects, functional definitions, tissues or selection. A null protection test does not validate selection.\n");
    let mut w = writer(&format!("{out}/functional_power_audit.md"))?;
    write!(w, "{report}")?;
    Ok(())
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
    fn audit_linear_weights_reproduce_frozen_ridge_and_holdout() {
        let rows: Vec<_> = (0..15)
            .map(|i| Recovery {
                gene: format!("gene{i}"),
                features: vec![1.0, i as f64 / 15.0, (i as f64).sin()],
                cost: (i as f64 * 0.7).cos(),
                outcomes: [Some(0.2 + 0.03 * i as f64), None],
            })
            .collect();
        let features: Vec<_> = rows.iter().map(|r| r.features.clone()).collect();
        let q: Vec<_> = rows.iter().map(|r| r.cost).collect();
        let y: Vec<_> = rows.iter().map(|r| r.outcomes[0].unwrap()).collect();
        for penalty in [0.001, 0.01, 0.1] {
            for held in 0..rows.len() {
                let fold = AuditFold::new(&features, held, penalty).unwrap();
                let train: Vec<_> = rows
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != held)
                    .map(|(_, r)| r)
                    .collect();
                let weights = fold.weights(&q);
                assert_eq!(weights[held], 0.0);
                assert!(
                    (audit_dot(&weights, &y)
                        - audit_fit_predict(&train, &rows[held], 0, penalty, true).unwrap())
                    .abs()
                        < 1e-10
                );
                assert!(
                    (audit_dot(&fold.baseline_weights(rows.len()), &y)
                        - audit_fit_predict(&train, &rows[held], 0, penalty, false).unwrap())
                    .abs()
                        < 1e-10
                );
            }
        }
    }
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
