//! Public-data falsification tests. Measurement units, biological replicates and
//! external feature definitions remain explicit; no age cutoff is called identity loss.
use super::*;
use rand::{seq::SliceRandom, Rng, SeedableRng};

fn soft_samples(path: &str) -> Result<Vec<BTreeMap<String, String>>> {
    let mut samples = vec![];
    let mut current = BTreeMap::new();
    for line in reader(path)?.lines() {
        let line = line?;
        if let Some(id) = line.strip_prefix("^SAMPLE = ") {
            if !current.is_empty() {
                samples.push(std::mem::take(&mut current));
            }
            current.insert("accession".into(), id.into());
        } else if let Some(title) = line.strip_prefix("!Sample_title = ") {
            current.insert("title".into(), title.into());
        } else if let Some(a) = line.strip_prefix("!Sample_characteristics_ch1 = ") {
            if let Some((k, v)) = a.split_once(": ") {
                current.insert(k.to_lowercase(), v.into());
            }
        }
    }
    if !current.is_empty() {
        samples.push(current);
    }
    Ok(samples)
}
pub(crate) fn quantile(mut values: Vec<f64>, q: f64) -> f64 {
    values.retain(|x| x.is_finite());
    values.sort_by(f64::total_cmp);
    if values.is_empty() {
        return f64::NAN;
    }
    values[((values.len() - 1) as f64 * q).round() as usize]
}
fn interval(values: &[f64], seed: u64) -> (f64, f64) {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let draws: Vec<_> = (0..4000)
        .map(|_| {
            (0..values.len())
                .map(|_| values[rng.gen_range(0..values.len())])
                .sum::<f64>()
                / values.len() as f64
        })
        .collect();
    (quantile(draws.clone(), 0.025), quantile(draws, 0.975))
}
fn correlation(a: &[f64], b: &[f64]) -> f64 {
    if a.len() != b.len() || a.len() < 3 {
        return f64::NAN;
    }
    let (ma, mb) = (mean(a), mean(b));
    let cov: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let va: f64 = a.iter().map(|x| (x - ma).powi(2)).sum();
    let vb: f64 = b.iter().map(|x| (x - mb).powi(2)).sum();
    cov / (va * vb).sqrt()
}
fn pathway_genes(root: &str) -> Result<BTreeSet<String>> {
    let output = std::process::Command::new("unzip")
        .args(["-p", &format!("{root}/ReactomePathways.gmt.zip")])
        .output()?;
    if !output.status.success() {
        return Err("Reactome archive extraction failed".into());
    }
    let text = String::from_utf8(output.stdout)?;
    let row = text
        .lines()
        .find(|l| l.starts_with("Extracellular matrix organization\tR-HSA-1474244\t"))
        .ok_or("external ECM pathway missing")?;
    Ok(row.split('\t').skip(2).map(String::from).collect())
}
fn promoter_map(root: &str) -> Result<BTreeMap<String, Vec<String>>> {
    let mut table = false;
    let mut indices = None;
    let mut map = BTreeMap::new();
    for line in reader(&format!("{root}/GSE73115_family.soft.gz"))?.lines() {
        let line = line?;
        if line == "!platform_table_begin" {
            table = true;
            continue;
        }
        if line == "!platform_table_end" {
            break;
        }
        if !table {
            continue;
        }
        let a = fields(&line, '\t');
        if indices.is_none() {
            indices = Some((
                a.iter()
                    .position(|s| s == "UCSC_RefGene_Name")
                    .ok_or("missing gene names")?,
                a.iter()
                    .position(|s| s == "UCSC_RefGene_Group")
                    .ok_or("missing regulatory groups")?,
            ));
            continue;
        }
        let (names, groups) = indices.unwrap();
        let genes: Vec<_> = a[names]
            .split(';')
            .zip(a[groups].split(';'))
            .filter(|(_, g)| g.starts_with("TSS") || *g == "5'UTR")
            .map(|(g, _)| g.to_string())
            .filter(|g| !g.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if !genes.is_empty() {
            map.insert(a[0].clone(), genes);
        }
    }
    Ok(map)
}
fn feature_stratum(f: &Feature) -> Option<String> {
    let sequence = f.sequence?;
    if f.island.is_empty() || f.island == "Unknown" {
        return None;
    }
    Some(format!(
        "{}:{}:{}:{}",
        f.chrom,
        f.island,
        (f.density * 50.0).floor(),
        (sequence * 5.0).floor()
    ))
}

fn prepare_functional(root: &str, out: &str) -> Result<()> {
    let genes = pathway_genes(root)?;
    let promoters = promoter_map(root)?;
    let features = load_features()?;
    let mut gene_counts: BTreeMap<String, usize> = BTreeMap::new();
    for (cpg, names) in &promoters {
        if features
            .get(cpg)
            .is_none_or(|f| feature_stratum(f).is_none())
        {
            continue;
        }
        for name in names.iter().filter(|g| genes.contains(*g)) {
            *gene_counts.entry(name.clone()).or_default() += 1;
        }
    }
    let mut functional = BTreeMap::new();
    for (cpg, names) in &promoters {
        if features
            .get(cpg)
            .is_none_or(|f| feature_stratum(f).is_none())
        {
            continue;
        }
        let weight: f64 = names
            .iter()
            .filter_map(|g| gene_counts.get(g))
            .map(|n| 1.0 / *n as f64)
            .sum();
        if weight > 0.0 {
            functional.insert(cpg.clone(), weight);
        }
    }
    let mut pools: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for cpg in promoters.keys().filter(|c| !functional.contains_key(*c)) {
        if let Some(key) = features.get(cpg).and_then(feature_stratum) {
            pools.entry(key).or_default().push(cpg.clone());
        }
    }
    let mut rng = rand::rngs::StdRng::seed_from_u64(20261003);
    let mut weights: BTreeMap<String, Vec<(usize, f64)>> = BTreeMap::new();
    let mut supported = 0;
    let mut gene_weight = 0.0;
    // One fixed matched control for every functional probe, per random mask. Sampling
    // with replacement preserves the target's external strata and assigned gene weight.
    for (cpg, weight) in &functional {
        let key = feature_stratum(&features[cpg]).unwrap();
        let Some(pool) = pools.get(&key).filter(|p| !p.is_empty()) else {
            continue;
        };
        supported += 1;
        gene_weight += weight;
        weights.entry(cpg.clone()).or_default().push((1, *weight));
        for mask in 0..20 {
            let control = pool.choose(&mut rng).unwrap();
            weights
                .entry(control.clone())
                .or_default()
                .push((mask + 3, *weight));
        }
    }
    for (cpg, f) in &features {
        if f.esl && feature_stratum(f).is_some() {
            weights.entry(cpg.clone()).or_default().push((2, 1.0));
        }
    }
    let mut manifest = writer(&format!("{out}/functional_mask_manifest.csv"))?;
    writeln!(manifest,"source,pathway,genes_in_source,annotated_genes,eligible_promoter_probes,matched_supported_probes,gene_weight,random_masks,status")?;
    writeln!(manifest,"Reactome_frozen_download,R-HSA-1474244,{},{},{},{supported},{gene_weight},20,ECM_functional_program_not_identity_threshold",genes.len(),gene_counts.len(),functional.len())?;
    let metadata = fibro_metadata(root)?;
    let mut r = reader(&format!("{root}/GSE179847_matrix.csv.gz"))?;
    let mut line = String::new();
    r.read_line(&mut line)?;
    let header = fields(&line, ',');
    let mut columns = vec![];
    let mut samples = vec![];
    for (i, name) in header
        .iter()
        .enumerate()
        .filter(|(_, s)| s.ends_with(" beta"))
    {
        let s = metadata
            .get(name.trim_end_matches(" beta"))
            .ok_or("DNA sample metadata absent")?;
        if s.treatment == "Control" {
            columns.push(i);
            samples.push(s.clone());
        }
    }
    let mut initial: BTreeMap<String, usize> = BTreeMap::new();
    for (i, s) in samples.iter().enumerate() {
        if !initial.contains_key(&s.lineage) || s.time < samples[initial[&s.lineage]].time {
            initial.insert(s.lineage.clone(), i);
        }
    }
    let mut sums = vec![vec![(0.0, 0.0, 0.0); 23]; samples.len()];
    let mut total = [0.0; 23];
    let mut index = 0;
    let mut trajectory = writer("data/derived/validation_trajectories.tsv")?;
    writeln!(trajectory, "cpg\tsample\tdonor\tlineage\ttime\tbeta")?;
    loop {
        line.clear();
        if r.read_line(&mut line)? == 0 {
            break;
        }
        let a = fields(&line, ',');
        if a.len() != header.len() {
            return Err("DNA matrix width differs".into());
        }
        let mut membership = weights.get(&a[0]).cloned().unwrap_or_default();
        if index % 256 == 0 {
            membership.push((0, 1.0));
        }
        index += 1;
        if membership.is_empty() {
            continue;
        }
        let values: Vec<_> = columns
            .iter()
            .map(|j| {
                numeric(&a[*j]).filter(|x| {
                    (0.0..=1.0).contains(x) && numeric(&a[*j + 1]).is_some_and(|p| p <= 0.01)
                })
            })
            .collect();
        for (mask, weight) in &membership {
            total[*mask] += weight;
        }
        for (j, s) in samples.iter().enumerate() {
            if membership.iter().any(|m| m.0 == 0) {
                if let Some(beta) = values[j] {
                    writeln!(
                        trajectory,
                        "{}\t{}\t{}\t{}\t{}\t{beta}",
                        a[0], s.key, s.donor, s.lineage, s.time
                    )?;
                }
            }
            if let (Some(beta), Some(start)) = (values[j], values[initial[&s.lineage]]) {
                for (mask, weight) in &membership {
                    let sum = &mut sums[j][*mask];
                    sum.0 += weight * (beta - start).abs();
                    sum.1 += weight;
                    sum.2 += weight * start;
                }
            }
        }
    }
    let mut w = writer("data/derived/functional_distances.tsv")?;
    write!(w, "sample\tdonor\tlineage\ttime\tinitial_time")?;
    for mask in 0..23 {
        write!(w, "\tdistance_{mask}\tcoverage_{mask}\tinitial_beta_{mask}")?;
    }
    writeln!(w)?;
    for (j, s) in samples.iter().enumerate() {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            s.key, s.donor, s.lineage, s.time, samples[initial[&s.lineage]].time
        )?;
        for (mask, (sum, observed, beta)) in sums[j].iter().enumerate() {
            write!(
                w,
                "\t{}\t{}\t{}",
                if *observed >= 0.9 * total[mask] && *observed > 0.0 {
                    sum / observed
                } else {
                    f64::NAN
                },
                observed / total[mask],
                beta / observed
            )?;
        }
        writeln!(w)?;
    }
    println!("Functional masks: {} external genes, {supported} matched promoter probes; {} control DNA arrays",gene_counts.len(),samples.len());
    Ok(())
}

pub fn run(root: &str, out: &str, command: &str) -> Result<()> {
    match command {
        "broad" => crate::broad_validation::run(root, out)?,
        "power" => crate::broad_validation::power_audit(out)?,
        "prepare"=>prepare_functional(root,out)?,
        "rna"=>{functional_rna(root,out,8.0)?;functional_rna(root,out,32.0)?;},
        "kinetics"=>maintenance_kinetics(out)?,
        "clones"=>{clones(root,out,0.95)?;clones(root,out,0.8)?;},
        "perturbations"=>perturbations(root,out)?,
        "distributions"=>trajectory_distributions(out)?,
        "identifiability"=>identifiability(out)?,
        "recovery"=>{recovery_pilot(root,out,false)?;recovery_pilot(root,out,true)?;},
        "report"=>extended_report(out)?,
        "all"=>{
            prepare_functional(root,out)?;
            functional_rna(root,out,8.0)?;functional_rna(root,out,32.0)?;
            maintenance_kinetics(out)?;
            clones(root,out,0.95)?;clones(root,out,0.8)?;
            perturbations(root,out)?;
            trajectory_distributions(out)?;
            identifiability(out)?;
            extended_report(out)?;
        }
        _=>return Err("validation subcommand: prepare|rna|kinetics|clones|perturbations|distributions|identifiability|recovery|report|all".into()),
    }
    Ok(())
}

#[derive(Clone)]
struct FunctionalRow {
    donor: String,
    x: Vec<f64>,
    distances: Vec<f64>,
    outcome: f64,
}
fn cv_functional(rows: &[FunctionalRow], out: &str, label: &str) -> Result<()> {
    let donors: BTreeSet<_> = rows.iter().map(|r| r.donor.clone()).collect();
    let mut w = writer(&format!("{out}/{label}_prediction.csv"))?;
    writeln!(
        w,
        "outcome,held_out_donor,n,penalty,model,baseline_mse,model_mse,mse_gain_vs_baseline"
    )?;
    // Analyze all fixed penalties and controls; no outcome-guided mask selection.
    for donor in &donors {
        for penalty in [0.001, 0.01, 0.1] {
            for model in 0..24 {
                let dimensions = rows[0].x.len() + usize::from(model > 0);
                let mut normal = NormalEquations::new(dimensions);
                let design = |r: &FunctionalRow| {
                    let mut x = r.x.clone();
                    if model > 0 {
                        x.push(r.distances[model - 1] * 10.0);
                    }
                    x
                };
                let mut baseline = NormalEquations::new(rows[0].x.len());
                for r in rows.iter().filter(|r| &r.donor != donor) {
                    // A donor contributes equal total weight, independent of culture/timepoint counts.
                    let n = rows.iter().filter(|s| s.donor == r.donor).count();
                    normal.add_xy(&design(r), r.outcome, 1.0 / n as f64);
                    baseline.add_xy(&r.x, r.outcome, 1.0 / n as f64);
                }
                let (b, _) = normal.fit(dimensions, penalty)?;
                let (b0, _) = baseline.fit(rows[0].x.len(), penalty)?;
                let (mut sum, mut sum0, mut n) = (0.0, 0.0, 0);
                for r in rows.iter().filter(|r| &r.donor == donor) {
                    let mu = b.iter().zip(design(r)).map(|(a, b)| a * b).sum::<f64>();
                    let mu0 = b0.iter().zip(&r.x).map(|(a, b)| a * b).sum::<f64>();
                    sum += (r.outcome - mu).powi(2);
                    sum0 += (r.outcome - mu0).powi(2);
                    n += 1;
                }
                let name = match model {
                    0 => "baseline".into(),
                    1 => "unweighted".into(),
                    2 => "functional_ECM".into(),
                    3 => "ESL".into(),
                    _ => format!("matched_random_{}", model - 4),
                };
                writeln!(
                    w,
                    "{label},{donor},{n},{penalty},{name},{},{},{}",
                    sum0 / n as f64,
                    sum / n as f64,
                    (sum0 - sum) / n as f64
                )?;
            }
        }
    }
    Ok(())
}
fn functional_rna(root: &str, out: &str, tolerance: f64) -> Result<()> {
    let genes = pathway_genes(root)?;
    let metadata = soft_samples(&format!("{root}/GSE179848_family.soft.gz"))?;
    let mut samples = BTreeMap::new();
    for m in metadata {
        if m.get("treatments").is_none_or(|s| s != "Control") {
            continue;
        }
        let (Some(id), Some(time), Some(lineage), Some(donor)) = (
            m.get("rnaseq_sampleid"),
            m.get("days_grown_udays").and_then(|s| numeric(s)),
            m.get("cell_line_group"),
            m.get("cell_line"),
        ) else {
            continue;
        };
        samples.insert(
            format!("Sample_{id}"),
            (donor.clone(), lineage.clone(), time),
        );
    }
    let mut r = reader(&format!("{root}/GSE179848_expression.csv.gz"))?;
    let mut line = String::new();
    r.read_line(&mut line)?;
    let header = fields(&line, ',');
    let columns: Vec<_> = header
        .iter()
        .enumerate()
        .filter(|(_, s)| samples.contains_key(*s))
        .map(|(i, s)| (i, samples[s].clone()))
        .collect();
    let mut ecm = vec![0.0; columns.len()];
    let mut ecm_expression = vec![vec![]; columns.len()];
    let mut gene_n = 0;
    let mut all_n = 0;
    let mut expression = vec![vec![]; columns.len()];
    while {
        line.clear();
        r.read_line(&mut line)? > 0
    } {
        let a = fields(&line, ',');
        if a.len() != header.len() {
            return Err("RNA matrix width differs".into());
        }
        let values: Vec<_> = columns.iter().map(|(i, _)| numeric(&a[*i])).collect();
        if values.iter().any(Option::is_none) {
            continue;
        }
        if genes.contains(&a[0]) {
            for (sum, value) in ecm.iter_mut().zip(&values) {
                *sum += value.unwrap();
            }
            for (v, x) in ecm_expression.iter_mut().zip(&values) {
                v.push(x.unwrap());
            }
            gene_n += 1;
        }
        // Outcome-independent systematic transcript subset for global remodelling.
        if all_n % 16 == 0 {
            for (v, x) in expression.iter_mut().zip(&values) {
                v.push(x.unwrap());
            }
        }
        all_n += 1;
    }
    if gene_n < 100 {
        return Err("insufficient external ECM gene RNA coverage".into());
    }
    for x in &mut ecm {
        *x /= gene_n as f64;
    }
    let mut by_line: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (j, (_, (_, lineage, _))) in columns.iter().enumerate() {
        by_line.entry(lineage.clone()).or_default().push(j);
    }
    for js in by_line.values_mut() {
        js.sort_by(|a, b| columns[*a].1 .2.total_cmp(&columns[*b].1 .2));
    }
    let mut distances = vec![];
    let mut lines = reader("data/derived/functional_distances.tsv")?.lines();
    lines.next();
    for line in lines {
        let a = fields(&line?, '\t');
        let d: Vec<_> = (0..23)
            .map(|j| a[5 + 3 * j].parse::<f64>())
            .collect::<std::result::Result<_, _>>()?;
        if d.iter().any(|x| !x.is_finite()) {
            continue;
        }
        distances.push((
            a[1].clone(),
            a[2].clone(),
            a[3].parse::<f64>()?,
            a[4].parse::<f64>()?,
            d,
            a[7].parse::<f64>()?,
            a[10].parse::<f64>()?,
        ));
    }
    let mut ecm_rows = vec![];
    let mut global_rows = vec![];
    let mut ecm_absolute_rows = vec![];
    let mut w = writer(&format!("{out}/rna_pairs_{tolerance}.csv"))?;
    writeln!(w,"donor,lineage,DNA_time,RNA_time,forecast_horizon,initial_RNA_time,external_ECM_genes,ECM_expression_change,global_expression_absolute_change,unweighted_distance,functional_distance,ESL_distance")?;
    for (donor, lineage, time, initial_time, d, baseline_beta, functional_beta) in distances {
        let Some(js) = by_line.get(&lineage) else {
            continue;
        };
        let Some(&start) = js.iter().min_by(|a, b| {
            (columns[**a].1 .2 - initial_time)
                .abs()
                .total_cmp(&(columns[**b].1 .2 - initial_time).abs())
        }) else {
            continue;
        };
        if (columns[start].1 .2 - initial_time).abs() > tolerance
            || time - initial_time < 16.0
            || columns[start].1 .2 > time - 8.0
        {
            continue;
        }
        // Each DNA sample predicts the next RNA sample, 1..16 days later.
        let Some(&next) = js
            .iter()
            .find(|j| columns[**j].1 .2 > time + 1.0 && columns[**j].1 .2 <= time + 16.0)
        else {
            continue;
        };
        let horizon = columns[next].1 .2 - time;
        let outcome = ecm[next] - ecm[start];
        let global = mean(
            &expression[next]
                .iter()
                .zip(&expression[start])
                .map(|(a, b)| (a - b).abs())
                .collect::<Vec<_>>(),
        );
        let x = vec![
            1.0,
            (time - initial_time) / 100.0,
            horizon / 16.0,
            f64::from(donor.starts_with("SURF1")),
            f64::from(lineage.ends_with("ox3")),
            ecm[start] / 10.0,
            baseline_beta,
            functional_beta,
        ];
        let base = FunctionalRow {
            donor: donor.clone(),
            x,
            distances: d.clone(),
            outcome,
        };
        ecm_rows.push(base.clone());
        let mut absolute = base.clone();
        absolute.outcome = mean(
            &ecm_expression[next]
                .iter()
                .zip(&ecm_expression[start])
                .map(|(a, b)| (a - b).abs())
                .collect::<Vec<_>>(),
        );
        ecm_absolute_rows.push(absolute);
        let mut g = base;
        g.outcome = global;
        global_rows.push(g);
        writeln!(
            w,
            "{donor},{lineage},{time},{},{horizon},{},{gene_n},{outcome},{global},{},{},{}",
            columns[next].1 .2, columns[start].1 .2, d[0], d[1], d[2]
        )?;
    }
    if ecm_rows.len() < 20 {
        return Err(format!(
            "only {} future RNA pairs: insufficient for planned CV",
            ecm_rows.len()
        )
        .into());
    }
    cv_functional(&ecm_rows, out, &format!("RNA_ECM_{tolerance}"))?;
    cv_functional(
        &ecm_absolute_rows,
        out,
        &format!("RNA_ECM_absolute_{tolerance}"),
    )?;
    cv_functional(&global_rows, out, &format!("RNA_global_{tolerance}"))?;
    println!(
        "RNA validation: {} forward pairs, {} donors, {gene_n} external ECM genes",
        ecm_rows.len(),
        ecm_rows
            .iter()
            .map(|r| &r.donor)
            .collect::<BTreeSet<_>>()
            .len()
    );
    Ok(())
}

fn binomial_score(m: f64, u: f64, p: f64) -> f64 {
    let p = p.clamp(1e-9, 1.0 - 1e-9);
    m * p.ln() + u * (1.0 - p).ln()
}
fn pulse_probability(rate: f64, fraction: f64, time: f64, uniform: bool) -> f64 {
    if uniform {
        fraction * (1.0 - (-rate * time).exp() * (-(-rate).exp_m1()) / rate)
    } else {
        fraction * (-(-rate * (time + 0.5)).exp_m1())
    }
}
fn maintenance_kinetics(out: &str) -> Result<()> {
    let mut w = writer(&format!("{out}/maintenance_holdout.csv"))?;
    writeln!(w,"chrom,position,model,held_out_hour,methylated,total,predicted_beta,observed_beta,log_score_per_read,score_gain_vs_early_constant,rate,fraction,boundary_fit")?;
    let mut lines = reader("data/derived/maintenance_counts.tsv")?.lines();
    lines.next();
    let mut n = 0;
    for line in lines {
        let a = fields(&line?, '\t');
        let counts: Vec<_> = a[2..]
            .iter()
            .map(|s| s.parse::<f64>())
            .collect::<std::result::Result<_, _>>()?;
        // Coverage selection concerns counts only, and is reported; no response screening.
        if (0..4).any(|j| counts[2 * j] + counts[2 * j + 1] < 5.0) || counts[0] + counts[1] < 10.0 {
            continue;
        }
        for heldout in [2, 3] {
            let times = [0.0, 1.0, 4.0, 16.0];
            let train: Vec<_> = (0..4).filter(|j| *j != heldout).collect();
            let early = counts[0] / (counts[0] + counts[1]);
            let sum_m: f64 = train.iter().map(|j| counts[2 * j]).sum();
            let sum_n: f64 = train
                .iter()
                .map(|j| counts[2 * j] + counts[2 * j + 1])
                .sum();
            let constant = sum_m / sum_n;
            for uniform in [false, true] {
                let objective = |p: &[f64]| {
                    -train
                        .iter()
                        .map(|j| {
                            binomial_score(
                                counts[2 * j],
                                counts[2 * j + 1],
                                pulse_probability(p[0].exp(), p[1], times[*j], uniform),
                            )
                        })
                        .sum::<f64>()
                };
                let (_, fit) = optimize(
                    objective,
                    &[vec![-1.0, 0.8], vec![1.0, 0.5]],
                    &[(-6.0, 4.0), (0.001, 0.999)],
                );
                let predicted = pulse_probability(fit[0].exp(), fit[1], times[heldout], uniform);
                let m = counts[2 * heldout];
                let u = counts[2 * heldout + 1];
                let total = m + u;
                let score = binomial_score(m, u, predicted) / total;
                let model = if uniform {
                    "pulse_uniform_CTMC_equivalent"
                } else {
                    "pulse_midpoint_CTMC_equivalent"
                };
                writeln!(
                    w,
                    "{},{},{model},{},{m},{total},{predicted},{},{score},{},{},{},{}",
                    a[0],
                    a[1],
                    times[heldout],
                    m / total,
                    score - binomial_score(m, u, early) / total,
                    fit[0].exp(),
                    fit[1],
                    fit[0] < -5.95 || fit[0] > 3.95 || fit[1] < 0.006 || fit[1] > 0.994
                )?;
            }
            let m = counts[2 * heldout];
            let u = counts[2 * heldout + 1];
            let total = m + u;
            for (model, predicted) in [("early_constant", early), ("pooled_constant", constant)] {
                let score = binomial_score(m, u, predicted) / total;
                writeln!(
                    w,
                    "{},{},{model},{},{m},{total},{predicted},{},{score},{},NaN,NaN,false",
                    a[0],
                    a[1],
                    times[heldout],
                    m / total,
                    score - binomial_score(m, u, early) / total
                )?;
            }
        }
        n += 1;
    }
    println!("Direct maintenance: {n} covered CpGs on chr1/22, independent 4h/16h holdouts");
    Ok(())
}

#[derive(Clone)]
struct PerturbationSite {
    chrom: u32,
    initial: f64,
    counts: Vec<(f64, f64)>,
}
fn perturbations(root: &str, out: &str) -> Result<()> {
    let mut r = reader(&format!("{root}/GSE145698_counts.txt.gz"))?;
    let mut line = String::new();
    r.read_line(&mut line)?;
    let header = fields(&line, '\t');
    if header.len() != 37 {
        return Err("unexpected perturbation matrix width".into());
    }
    let mut sites = vec![];
    let mut index = 0;
    let mut total = 0;
    while {
        line.clear();
        r.read_line(&mut line)? > 0
    } {
        let take = index % 32 == 0;
        index += 1;
        if !take {
            continue;
        }
        let a = fields(&line, '\t');
        if a.len() != header.len() {
            return Err("perturbation row width differs".into());
        }
        let Some((chrom, _)) = a[0].split_once(':') else {
            continue;
        };
        let Ok(chrom) = chrom.trim_start_matches("chr").parse::<u32>() else {
            continue;
        };
        if !(1..=19).contains(&chrom) {
            continue;
        }
        let values: Vec<_> = a[1..].iter().map(|s| numeric(s)).collect();
        if values.iter().any(Option::is_none) {
            continue;
        }
        let counts: Vec<_> = values
            .chunks_exact(2)
            .map(|v| (v[0].unwrap(), v[1].unwrap()))
            .collect();
        if counts.iter().any(|(m, n)| *m < 0.0 || m > n || *n < 0.0) {
            return Err("invalid methylated/total perturbation counts".into());
        }
        let m: f64 = counts[..3].iter().map(|v| v.0).sum();
        let n: f64 = counts[..3].iter().map(|v| v.1).sum();
        if n < 15.0 || counts[..3].iter().filter(|v| v.1 >= 5.0).count() < 2 {
            continue;
        }
        total += 1;
        sites.push(PerturbationSite {
            chrom,
            initial: (m + 0.5) / (n + 1.0),
            counts: counts[3..].to_vec(),
        });
    }
    let mut w = writer(&format!("{out}/perturbation_prediction.csv"))?;
    writeln!(w,"mutation,replicate,held_out_chromosome_group,model,sites,mean_log_score_per_read,score_gain_vs_unchanged,MSE,mean_observed_minus_WT,parameter_a,parameter_b")?;
    let mut effects = writer(&format!("{out}/perturbation_effects.csv"))?;
    writeln!(effects,"mutation,replicate,sites,mean_WT_beta,mean_mutant_beta,mean_change,fraction_decreased,correlation_WT_mutant")?;
    for mutation in 0..5 {
        let name = header[1 + 2 * (3 + 3 * mutation)].trim_end_matches("-rep1.mC");
        for genome in 0..2 {
            let train: Vec<_> = sites.iter().filter(|s| s.chrom % 2 != genome).collect();
            let objective = |p: &[f64], affine: bool| {
                let mut sum = 0.0;
                let mut n = 0;
                for s in &train {
                    let predicted = if affine {
                        p[0] + p[1] * s.initial
                    } else {
                        s.initial * (-p[0]).exp()
                    };
                    for (m, total) in s.counts[mutation * 3..mutation * 3 + 3]
                        .iter()
                        .filter(|v| v.1 >= 5.0)
                    {
                        sum -= binomial_score(*m, total - m, predicted) / total;
                        n += 1;
                    }
                }
                sum / n as f64
            };
            let loss = optimize(
                |p| objective(p, false),
                &[vec![0.1], vec![1.0]],
                &[(0.0, 6.0)],
            )
            .1;
            let affine = optimize(
                |p| objective(p, true),
                &[vec![0.0, 0.8], vec![0.1, 0.5]],
                &[(-0.5, 0.5), (0.0, 1.5)],
            )
            .1;
            let components = match mutation {
                2 => Some((1, 0)),
                4 => Some((3, 0)),
                _ => None,
            };
            let combination_loss = components.map(|(left, right)| {
                [left, right]
                    .iter()
                    .map(|component| {
                        let objective = |p: &[f64]| {
                            let mut sum = 0.0;
                            let mut n = 0;
                            for s in &train {
                                for (m, total) in s.counts[component * 3..component * 3 + 3]
                                    .iter()
                                    .filter(|v| v.1 >= 5.0)
                                {
                                    sum -= binomial_score(*m, total - m, s.initial * (-p[0]).exp())
                                        / total;
                                    n += 1;
                                }
                            }
                            sum / n as f64
                        };
                        optimize(objective, &[vec![0.1], vec![1.0]], &[(0.0, 6.0)]).1[0]
                    })
                    .sum::<f64>()
            });
            for replicate in 0..3 {
                for model in [
                    "unchanged",
                    "CTMC_effective_loss",
                    "affine_empirical",
                    "unseen_combination_independent_effects",
                ] {
                    if model == "unseen_combination_independent_effects"
                        && combination_loss.is_none()
                    {
                        continue;
                    }

                    let (mut sum, mut baseline, mut mse, mut effect, mut n) =
                        (0.0, 0.0, 0.0, 0.0, 0);
                    for s in sites.iter().filter(|s| s.chrom % 2 == genome) {
                        let (m, total) = s.counts[mutation * 3 + replicate];
                        if total < 5.0 {
                            continue;
                        }
                        let p = match model {
                            "unchanged" => s.initial,
                            "CTMC_effective_loss" => s.initial * (-loss[0]).exp(),
                            "unseen_combination_independent_effects" => {
                                s.initial * (-combination_loss.unwrap()).exp()
                            }
                            _ => affine[0] + affine[1] * s.initial,
                        }
                        .clamp(1e-9, 1.0 - 1e-9);
                        sum += binomial_score(m, total - m, p) / total;
                        baseline += binomial_score(m, total - m, s.initial) / total;
                        mse += (m / total - p).powi(2);
                        effect += m / total - s.initial;
                        n += 1;
                    }
                    let (a, b) = match model {
                        "unchanged" => (0.0, 1.0),
                        "CTMC_effective_loss" => (loss[0], f64::NAN),
                        "unseen_combination_independent_effects" => {
                            (combination_loss.unwrap(), f64::NAN)
                        }
                        _ => (affine[0], affine[1]),
                    };
                    writeln!(
                        w,
                        "{name},{replicate},{genome},{model},{n},{},{},{},{},{a},{b}",
                        sum / n as f64,
                        (sum - baseline) / n as f64,
                        mse / n as f64,
                        effect / n as f64
                    )?;
                }
            }
        }
        for replicate in 0..3 {
            let paired: Vec<_> = sites
                .iter()
                .filter_map(|s| {
                    let (m, n) = s.counts[mutation * 3 + replicate];
                    (n >= 5.0).then_some((s.initial, m / n))
                })
                .collect();
            let x: Vec<_> = paired.iter().map(|p| p.0).collect();
            let y: Vec<_> = paired.iter().map(|p| p.1).collect();
            writeln!(
                effects,
                "{name},{replicate},{},{},{},{},{},{}",
                x.len(),
                mean(&x),
                mean(&y),
                mean(&y) - mean(&x),
                paired.iter().filter(|(a, b)| b < a).count() as f64 / x.len() as f64,
                correlation(&x, &y)
            )?;
        }
    }
    println!("Perturbation validation: {total} fixed-sampled WT-covered CpGs, five mutations x three biological cultures");
    Ok(())
}

#[derive(Clone)]
struct CloneCell {
    mouse: String,
    clone: String,
    x: Vec<f64>,
    y: Vec<f64>,
}
fn clones(root: &str, out: &str, quality_threshold: f64) -> Result<()> {
    let mut metadata = BTreeMap::new();
    let mut lines = reader("data/derived/clone_metadata.tsv")?.lines();
    let header = fields(&lines.next().ok_or("no clone metadata")??, '\t');
    let col = |name: &str| {
        header
            .iter()
            .position(|x| x == name)
            .ok_or("clone metadata field absent")
    };
    let (batch, larry, quality, experiment) = (
        col("ProcessingBatch")?,
        col("LARRY")?,
        col("PerformanceNonHhaI")?,
        col("Experiment")?,
    );
    for line in lines {
        let a = fields(&line?, '\t');
        if a[experiment] != "LARRY main experiment"
            || !["LARRY_mouse3", "LARRY_mouse4"].contains(&a[batch].as_str())
            || a[larry] == "NA"
            || numeric(&a[quality]).is_none_or(|q| q < quality_threshold)
        {
            continue;
        }
        metadata.insert(
            a[0].clone(),
            (
                a[batch].clone(),
                a[larry].clone(),
                numeric(&a[quality]).unwrap(),
            ),
        );
    }
    // Build protein covariates from independent UMI counts. In stained libraries,
    // absent marker rows are zero detected molecules; require a nonempty library.
    // This avoids selecting only rare cells positive for every marker.
    let mut raw: BTreeMap<String, (Vec<f64>, f64)> = metadata
        .keys()
        .map(|id| (id.clone(), (vec![0.0; 6], 0.0)))
        .collect();
    let markers = ["SCA1", "cKIT", "CD48", "CD150", "CD135", "CD201"];
    for (mouse, suffix) in [("mouse3", "larry3"), ("mouse4", "larry4")] {
        for line in reader(&format!("{root}/EPIClone_{mouse}_umi.tsv.gz"))?
            .lines()
            .skip(1)
        {
            let a = fields(&line?, '\t');
            let id = format!("{}_{}", a[0], suffix);
            let Some((values, total)) = raw.get_mut(&id) else {
                continue;
            };
            let count = a[2].parse::<f64>()?;
            *total += count;
            if let Some(j) = markers
                .iter()
                .position(|m| a[1].trim_end_matches("_totalseqb") == *m)
            {
                values[j] += count;
            }
        }
    }
    let proteins: BTreeMap<_, _> = raw
        .into_iter()
        .filter(|(_, (_, total))| *total > 0.0)
        .map(|(id, (values, total))| {
            let mut x = vec![1.0];
            x.extend(values.iter().map(|v| (1.0 + 10000.0 * v / total).ln()));
            x.push(total.ln() / 10.0);
            x.push(metadata[&id].2);
            (id, x)
        })
        .collect();
    let mut types = BTreeMap::new();
    for line in reader(&format!("{root}/EPIClone_cpg_selection.csv"))?
        .lines()
        .skip(1)
    {
        let a = fields(&line?, ',');
        types.insert(a[1].clone(), a[5].clone());
    }
    // Probe filtering uses the external undigested-library amplification control.
    let mut lines = reader(&format!("{root}/EPIClone_panel.tsv"))?.lines();
    let header = fields(&lines.next().ok_or("panel header absent")??, '\t');
    let dropout = header
        .iter()
        .position(|s| s == "DropoutUncutLKs")
        .ok_or("uncut amplification control missing")?;
    let mut reliable = BTreeSet::new();
    for line in lines {
        let a = fields(&line?, '\t');
        if numeric(&a[dropout + 1]).is_some_and(|p| p >= 0.95) {
            reliable.insert(a[0].clone());
        }
    }
    // The panel file has R row names preceding its named columns; detect offset below.
    if reliable.is_empty() {
        return Err("no assay probes pass undigested control".into());
    }
    let mut lines = reader("data/derived/clone_methylation.tsv")?.lines();
    let header = fields(
        &lines.next().ok_or("clone methylation header absent")??,
        '\t',
    );
    let indices: Vec<_> = header
        .iter()
        .enumerate()
        .skip(1)
        .filter(|(_, id)| reliable.contains(*id))
        .map(|(i, id)| (i, id.clone()))
        .collect();
    let mut cells = vec![];
    for line in lines {
        let a = fields(&line?, '\t');
        let Some(x) = proteins.get(&a[0]) else {
            continue;
        };
        let (mouse, clone, _) = &metadata[&a[0]];
        let y: Vec<_> = indices
            .iter()
            .map(|(i, _)| numeric(&a[*i]).unwrap_or(f64::NAN))
            .collect();
        cells.push(CloneCell {
            mouse: mouse.clone(),
            clone: clone.clone(),
            x: x.clone(),
            y,
        });
    }
    let mut w = writer(&format!("{out}/clone_conditioning_{quality_threshold}.csv"))?;
    writeln!(w,"mouse,probe,class,heldout_cells,training_clones,protein_model_log_score,protein_plus_clone_log_score,log_score_gain,cell_weighted_beta,equal_clone_beta,absolute_reweighting_change")?;
    let mut decomposition = writer(&format!(
        "{out}/clone_representation_{quality_threshold}.csv"
    ))?;
    writeln!(decomposition,"probe,class,shared_clones,cells_mouse3,cells_mouse4,bulk_difference,symmetric_within_clone_component,symmetric_representation_component,decomposition_error,status")?;
    for (probe, (_, id)) in indices.iter().enumerate() {
        let mut means: BTreeMap<String, [(f64, usize); 2]> = BTreeMap::new();
        for c in &cells {
            if c.y[probe].is_finite() {
                let mouse = usize::from(c.mouse == "LARRY_mouse4");
                let g = &mut means.entry(c.clone.clone()).or_default()[mouse];
                g.0 += c.y[probe];
                g.1 += 1;
            }
        }
        let shared: Vec<_> = means
            .values()
            .filter(|g| g[0].1 >= 5 && g[1].1 >= 5)
            .collect();
        if shared.len() < 3 {
            continue;
        }
        let n0: usize = shared.iter().map(|g| g[0].1).sum();
        let n1: usize = shared.iter().map(|g| g[1].1).sum();
        let (mut total, mut within, mut representation) = (0.0, 0.0, 0.0);
        for g in &shared {
            let a = g[0].1 as f64 / n0 as f64;
            let b = g[1].1 as f64 / n1 as f64;
            let m0 = g[0].0 / g[0].1 as f64;
            let m1 = g[1].0 / g[1].1 as f64;
            total += b * m1 - a * m0;
            within += 0.5 * (a + b) * (m1 - m0);
            representation += 0.5 * (m0 + m1) * (b - a);
        }
        writeln!(decomposition,"{id},{},{},{n0},{n1},{total},{within},{representation},{},different_recipients_not_longitudinal_selection",types.get(id).map(String::as_str).unwrap_or("Unknown"),shared.len(),total-within-representation)?;
    }
    for mouse in ["LARRY_mouse3", "LARRY_mouse4"] {
        let mouse_cells: Vec<_> = cells.iter().filter(|c| c.mouse == mouse).collect();
        // Both folds test different cells within independently barcode-labelled clones.
        for (probe, (_, id)) in indices.iter().enumerate() {
            let mut sum = 0.0;
            let mut sum_clone = 0.0;
            let mut n = 0;
            let mut train_clones = 0;
            let mut groups: BTreeMap<String, (f64, usize)> = BTreeMap::new();
            for c in &mouse_cells {
                if c.y[probe].is_finite() {
                    let g = groups.entry(c.clone.clone()).or_default();
                    g.0 += c.y[probe];
                    g.1 += 1;
                }
            }
            let beta = groups
                .values()
                .filter(|g| g.1 >= 10)
                .map(|g| g.0)
                .sum::<f64>()
                / groups
                    .values()
                    .filter(|g| g.1 >= 10)
                    .map(|g| g.1)
                    .sum::<usize>() as f64;
            let equal = mean(
                &groups
                    .values()
                    .filter(|g| g.1 >= 10)
                    .map(|g| g.0 / g.1 as f64)
                    .collect::<Vec<_>>(),
            );
            for fold in 0..2 {
                let mut normal = NormalEquations::new(cells[0].x.len());
                for (j, c) in mouse_cells
                    .iter()
                    .enumerate()
                    .filter(|(j, c)| j % 2 != fold && c.y[probe].is_finite())
                {
                    let _ = j;
                    normal.add_xy(&c.x, c.y[probe], 1.0);
                }
                let (fit, _) = normal.fit(cells[0].x.len(), 0.01)?;
                let predict = |c: &CloneCell| {
                    fit.iter()
                        .zip(&c.x)
                        .map(|(a, b)| a * b)
                        .sum::<f64>()
                        .clamp(0.01, 0.99)
                };
                let mut residuals: BTreeMap<String, (f64, usize)> = BTreeMap::new();
                for (_, c) in mouse_cells
                    .iter()
                    .enumerate()
                    .filter(|(j, c)| j % 2 != fold && c.y[probe].is_finite())
                {
                    let s = residuals.entry(c.clone.clone()).or_default();
                    s.0 += c.y[probe] - predict(c);
                    s.1 += 1;
                }
                train_clones = residuals.len();
                for (_, c) in mouse_cells
                    .iter()
                    .enumerate()
                    .filter(|(j, c)| j % 2 == fold && c.y[probe].is_finite())
                {
                    let p = predict(c);
                    let shift = residuals
                        .get(&c.clone)
                        .map(|(sum, n)| sum / (*n as f64 + 20.0))
                        .unwrap_or(0.0);
                    sum += binomial_score(c.y[probe], 1.0 - c.y[probe], p);
                    sum_clone +=
                        binomial_score(c.y[probe], 1.0 - c.y[probe], (p + shift).clamp(0.01, 0.99));
                    n += 1;
                }
            }
            writeln!(
                w,
                "{mouse},{id},{},{n},{train_clones},{},{},{},{beta},{equal},{}",
                types.get(id).map(String::as_str).unwrap_or("Unknown"),
                sum / n as f64,
                sum_clone / n as f64,
                (sum_clone - sum) / n as f64,
                (beta - equal).abs()
            )?;
        }
    }
    println!("Clone conditioning: {} high-quality cells, {} reliable probes, independent protein covariates and LARRY labels",cells.len(),indices.len());
    Ok(())
}

fn trajectory_distributions(out: &str) -> Result<()> {
    type TrajectoryKey = (String, String, String);
    let mut data: BTreeMap<TrajectoryKey, Vec<(String, f64, f64)>> = BTreeMap::new();
    let mut lines = reader("data/derived/validation_trajectories.tsv")?.lines();
    lines.next();
    for line in lines {
        let a = fields(&line?, '\t');
        data.entry((a[0].clone(), a[2].clone(), a[3].clone()))
            .or_default()
            .push((a[1].clone(), a[4].parse()?, a[5].parse()?));
    }
    type ResidualRecords = BTreeMap<(String, String, String), Vec<(String, f64, f64, f64)>>;
    let mut residuals: ResidualRecords = BTreeMap::new();
    let mut count = 0;
    for ((cpg, donor, _), rows) in &mut data {
        rows.sort_by(|a, b| a.1.total_cmp(&b.1));
        if rows.len() < 6 {
            continue;
        }
        let train_n = (rows.len() * 3 / 4).min(rows.len() - 1);
        let train = &rows[..train_n];
        let test = &rows[train_n..];
        let initial = train[0].2;
        let origin = train[0].1;
        let objective = |p: &[f64]| {
            train
                .iter()
                .skip(1)
                .map(|(_, time, beta)| {
                    (beta - probability(initial, p[0].exp(), p[1].exp(), (time - origin) / 100.0))
                        .powi(2)
                })
                .sum::<f64>()
        };
        let fit = optimize(
            objective,
            &[vec![-3.0, -3.0], vec![-5.0, -1.0]],
            &[(-10.0, 3.0); 2],
        )
        .1;
        let times: Vec<_> = train.iter().map(|r| (r.1 - origin) / 100.0).collect();
        let betas: Vec<_> = train.iter().map(|r| r.2).collect();
        let (mt, mb) = (mean(&times), mean(&betas));
        let slope = times
            .iter()
            .zip(&betas)
            .map(|(x, y)| (x - mt) * (y - mb))
            .sum::<f64>()
            / times
                .iter()
                .map(|x| (x - mt).powi(2))
                .sum::<f64>()
                .max(1e-12);
        for model in ["early_constant", "independent_CTMC", "linear_empirical"] {
            let predict = |time: f64| match model {
                "early_constant" => initial,
                "independent_CTMC" => {
                    probability(initial, fit[0].exp(), fit[1].exp(), (time - origin) / 100.0)
                }
                _ => (mb + slope * ((time - origin) / 100.0 - mt)).clamp(0.0, 1.0),
            };
            let v = (train
                .iter()
                .skip(1)
                .map(|(_, t, b)| (b - predict(*t)).powi(2))
                .sum::<f64>()
                / (train_n - 1) as f64)
                .max(0.02_f64.powi(2));
            for (sample, time, beta) in test {
                let error = beta - predict(*time);
                residuals
                    .entry((donor.clone(), sample.clone(), model.into()))
                    .or_default()
                    .push((cpg.clone(), error, error / v.sqrt(), v));
            }
        }
        count += 1;
    }
    let mut w = writer(&format!("{out}/trajectory_distribution_checks.csv"))?;
    writeln!(w,"donor,held_out_array,model,n_probes,RMSE,mean_error,residual_SD,mean_log_predictive_density,z_SD,fraction_abs_z_gt_1_96,normal_expected_fraction,fraction_abs_z_gt_3,normal_expected_3,abs_error_95")?;
    for ((donor, sample, model), rows) in &residuals {
        let errors: Vec<_> = rows.iter().map(|r| r.1).collect();
        let z: Vec<_> = rows.iter().map(|r| r.2).collect();
        writeln!(
            w,
            "{donor},{sample},{model},{},{},{},{},{},{},{},0.05,{},0.0026998,{}",
            rows.len(),
            mean(&errors.iter().map(|e| e * e).collect::<Vec<_>>()).sqrt(),
            mean(&errors),
            variance(&errors).sqrt(),
            mean(
                &rows
                    .iter()
                    .map(|r| normal_logpdf(r.1, 0.0, r.3))
                    .collect::<Vec<_>>()
            ),
            variance(&z).sqrt(),
            z.iter().filter(|x| x.abs() > 1.96).count() as f64 / z.len() as f64,
            z.iter().filter(|x| x.abs() > 3.0).count() as f64 / z.len() as f64,
            quantile(errors.iter().map(|x| x.abs()).collect(), 0.95)
        )?;
    }
    // Residual correlation uses strictly future arrays, centres each probe, and has
    // no biological interpretation without proliferation/batch/composition adjustment.
    let arrays: Vec<_> = residuals
        .keys()
        .filter(|k| k.2 == "independent_CTMC")
        .map(|k| (k.0.clone(), k.1.clone()))
        .collect();
    let mut by_probe: BTreeMap<String, Vec<Option<f64>>> = BTreeMap::new();
    for (j, (donor, sample)) in arrays.iter().enumerate() {
        for (cpg, error, _, _) in
            &residuals[&(donor.clone(), sample.clone(), "independent_CTMC".into())]
        {
            by_probe
                .entry(cpg.clone())
                .or_insert_with(|| vec![None; arrays.len()])[j] = Some(*error);
        }
    }
    let matrix: Vec<_> = by_probe
        .values()
        .filter(|v| v.iter().filter(|x| x.is_some()).count() >= arrays.len() * 4 / 5)
        .take(128)
        .map(|v| {
            let m = mean(&v.iter().flatten().copied().collect::<Vec<_>>());
            v.iter().map(|x| x.unwrap_or(m)).collect::<Vec<_>>()
        })
        .collect();
    if matrix.len() >= 10 && arrays.len() >= 6 {
        let transposed: Vec<_> = (0..arrays.len())
            .map(|i| matrix.iter().map(|v| v[i]).collect())
            .collect();
        let mut null_rng = rand::rngs::StdRng::seed_from_u64(20261003);
        let leading = |x: &[Vec<f64>]| -> f64 {
            let centred: Vec<_> = (0..x[0].len())
                .map(|j| {
                    let v: Vec<_> = x.iter().map(|r| r[j]).collect();
                    let m = mean(&v);
                    v.iter().map(|a| a - m).collect::<Vec<_>>()
                })
                .collect();
            let mut v = vec![1.0 / (centred.len() as f64).sqrt(); centred.len()];
            for _ in 0..50 {
                let mut next = vec![0.0; v.len()];
                for i in 0..v.len() {
                    for j in 0..v.len() {
                        next[i] += centred[i]
                            .iter()
                            .zip(&centred[j])
                            .map(|(a, b)| a * b)
                            .sum::<f64>()
                            * v[j];
                    }
                }
                let norm = next.iter().map(|x| x * x).sum::<f64>().sqrt();
                if norm == 0.0 {
                    return 0.0;
                }
                v = next.iter().map(|x| x / norm).collect();
            }
            let eigen: f64 = (0..v.len())
                .map(|i| {
                    v[i] * (0..v.len())
                        .map(|j| {
                            centred[i]
                                .iter()
                                .zip(&centred[j])
                                .map(|(a, b)| a * b)
                                .sum::<f64>()
                                * v[j]
                        })
                        .sum::<f64>()
                })
                .sum();
            eigen / centred.iter().flatten().map(|x| x * x).sum::<f64>()
        };
        let actual = leading(&transposed);
        let mut nulls = vec![];
        // Independently shuffle residuals within donor for each probe. Retains donor
        // offsets; tests shared within-donor array structure rather than donor biology.
        let mut permutations = transposed.clone();
        for _ in 0..100 {
            for j in 0..matrix.len() {
                for donor in arrays.iter().map(|r| &r.0).collect::<BTreeSet<_>>() {
                    let ids: Vec<_> = arrays
                        .iter()
                        .enumerate()
                        .filter(|(_, r)| &r.0 == donor)
                        .map(|(i, _)| i)
                        .collect();
                    let mut values: Vec<_> = ids.iter().map(|i| transposed[*i][j]).collect();
                    values.shuffle(&mut null_rng);
                    for (i, value) in ids.iter().zip(values) {
                        permutations[*i][j] = value;
                    }
                }
            }
            nulls.push(leading(&permutations));
        }
        let mut w = writer(&format!("{out}/residual_dependence.csv"))?;
        writeln!(w,"held_out_arrays,probes,leading_variance_fraction,within_donor_shuffle_lower,within_donor_shuffle_upper,permutation_p,status")?;
        writeln!(
            w,
            "{},{},{actual},{},{},{},descriptive_residual_structure_not_validated_M5",
            arrays.len(),
            matrix.len(),
            quantile(nulls.clone(), 0.025),
            quantile(nulls.clone(), 0.975),
            (1 + nulls.iter().filter(|x| **x >= actual).count()) as f64 / 101.0
        )?;
    }
    println!(
        "Distribution checks: {count} CpG/culture trajectories; {} future array/model summaries",
        residuals.len()
    );
    Ok(())
}

// Exact finite-state killed chain by uniformization. Three loci retain shared
// state-dependent killing, so selection is not replaced by independent locus removal.
fn killed_states(parameters: &[f64], time: f64) -> (Vec<f64>, f64) {
    let gain = parameters[0].exp();
    let loss = parameters[1].exp();
    let alpha = parameters[2];
    let gamma = parameters[3];
    let q = [0.0, 0.5, 1.0];
    let mut generator = vec![vec![0.0; 8]; 8];
    let mut lambda: f64 = 0.0;
    for (state, row) in generator.iter_mut().enumerate() {
        let deviation = (0..3)
            .filter(|i| state & (1 << i) != 0)
            .map(|i| q[i])
            .sum::<f64>()
            / 1.5;
        let mut total = 0.05 * (gamma * deviation).exp();
        for (site, qi) in q.iter().enumerate() {
            let rate = if state & (1 << site) == 0 {
                gain * (-alpha * qi).exp()
            } else {
                loss
            };
            row[state ^ (1 << site)] = rate;
            total += rate;
        }
        row[state] = -total;
        lambda = lambda.max(total);
    }
    let steps = (lambda * time / 10.0).ceil().max(1.0) as usize;
    let dt = time / steps as f64;
    let mut p = vec![0.0; 8];
    p[0] = 1.0;
    for _ in 0..steps {
        let mut term = p.clone();
        let mut weight = (-lambda * dt).exp();
        let mut accumulated: Vec<_> = term.iter().map(|x| x * weight).collect();
        for k in 1..150 {
            let mut next = vec![0.0; 8];
            for i in 0..8 {
                for (j, value) in next.iter_mut().enumerate() {
                    *value += term[i] * (f64::from(i == j) + generator[i][j] / lambda);
                }
            }
            term = next;
            weight *= lambda * dt / k as f64;
            for (p, x) in accumulated.iter_mut().zip(&term) {
                *p += weight * x;
            }
            if weight < 1e-14 && k as f64 > lambda * dt {
                break;
            }
        }
        p = accumulated;
    }
    let survival = p.iter().sum::<f64>();
    for x in &mut p {
        *x /= survival;
    }
    (p, survival)
}
fn observed_pattern(p: &[f64], pattern: usize, mask: usize) -> f64 {
    p.iter()
        .enumerate()
        .map(|(state, p)| {
            p * (0..3)
                .filter(|i| mask & (1 << i) != 0)
                .map(|i| {
                    if (state ^ (pattern)) & (1 << i) == 0 {
                        0.98
                    } else {
                        0.02
                    }
                })
                .product::<f64>()
        })
        .sum::<f64>()
}
type SimObservation = (f64, BTreeMap<(usize, usize), usize>, usize);
fn identifiability(out: &str) -> Result<()> {
    let mut w = writer(&format!("{out}/multisite_identifiability.csv"))?;
    writeln!(w,"truth,architecture,seed,candidate,training_negative_log_likelihood,test_log_predictive_density,estimated_gain,estimated_loss,estimated_protection,estimated_selection,estimated_survival,truth_survival")?;
    let mut rng = rand::rngs::StdRng::seed_from_u64(20261003);
    for (truth, alpha, gamma) in [
        ("protection", 1.5, 0.0),
        ("selection", 0.0, 3.0),
        ("both", 1.0, 2.0),
    ] {
        let actual = [0.15_f64.ln(), 0.3_f64.ln(), alpha, gamma];
        for architecture in [
            "endpoint_bulk",
            "longitudinal_cells",
            "longitudinal_cells_and_survival",
        ] {
            let times = if architecture == "endpoint_bulk" {
                vec![4.0]
            } else {
                vec![1.0, 2.0, 4.0]
            };
            for seed in 0..8 {
                let mut train = vec![];
                let mut test = vec![];
                for time in &times {
                    let (states, survival) = killed_states(&actual, *time);
                    for split in 0..2 {
                        let mut observations: BTreeMap<(usize, usize), usize> = BTreeMap::new();
                        for _ in 0..400 {
                            let draw = rng.gen::<f64>();
                            let mut cumulative = 0.0;
                            let mut state = 7;
                            for (i, p) in states.iter().enumerate() {
                                cumulative += p;
                                if draw <= cumulative {
                                    state = i;
                                    break;
                                }
                            }
                            if architecture == "endpoint_bulk" {
                                // Bulk yields three marginal locus frequencies, not joint patterns.
                                for site in 0..3 {
                                    if rng.gen::<f64>() < 0.6 {
                                        let call = usize::from(
                                            (state & (1 << site) != 0) ^ rng.gen_bool(0.02),
                                        ) << site;
                                        *observations.entry((1 << site, call)).or_default() += 1;
                                    }
                                }
                            } else {
                                let (mut mask, mut pattern) = (0, 0);
                                for site in 0..3 {
                                    if rng.gen::<f64>() < 0.6 {
                                        mask |= 1 << site;
                                        if (state & (1 << site) != 0) ^ rng.gen_bool(0.02) {
                                            pattern |= 1 << site;
                                        }
                                    }
                                }
                                if mask > 0 {
                                    *observations.entry((mask, pattern)).or_default() += 1;
                                }
                            }
                        }
                        let survived = (0..500).filter(|_| rng.gen_bool(survival)).count();
                        let row = (*time, observations, survived);
                        if split == 0 {
                            train.push(row);
                        } else {
                            test.push(row);
                        }
                    }
                }
                for candidate in ["independent", "protection", "selection", "both"] {
                    let bounds = [
                        (-4.0, 0.0),
                        (-4.0, 1.0),
                        if ["protection", "both"].contains(&candidate) {
                            (0.0, 4.0)
                        } else {
                            (0.0, 0.0)
                        },
                        if ["selection", "both"].contains(&candidate) {
                            (0.0, 5.0)
                        } else {
                            (0.0, 0.0)
                        },
                    ];
                    let score = |parameters: &[f64], rows: &[SimObservation]| {
                        rows.iter()
                            .map(|(time, observations, survived)| {
                                let (p, survival) = killed_states(parameters, *time);
                                let mut ll = observations
                                    .iter()
                                    .map(|((mask, pattern), n)| {
                                        *n as f64
                                            * observed_pattern(&p, *pattern, *mask).max(1e-12).ln()
                                    })
                                    .sum::<f64>();
                                if architecture.ends_with("and_survival") {
                                    ll += binomial_score(
                                        *survived as f64,
                                        (500 - survived) as f64,
                                        survival,
                                    );
                                }
                                -ll
                            })
                            .sum::<f64>()
                    };
                    let (objective, fit) = optimize(
                        |p| score(p, &train),
                        &[
                            vec![-2.0, -1.0, 0.0, 0.0],
                            vec![
                                -1.5,
                                -0.5,
                                if bounds[2].1 > 0.0 { 1.0 } else { 0.0 },
                                if bounds[3].1 > 0.0 { 2.0 } else { 0.0 },
                            ],
                        ],
                        &bounds,
                    );
                    let n: usize = test
                        .iter()
                        .map(|r| {
                            r.1.values().sum::<usize>()
                                + if architecture.ends_with("and_survival") {
                                    500
                                } else {
                                    0
                                }
                        })
                        .sum();
                    writeln!(w,"{truth},{architecture},{seed},{candidate},{objective},{},{},{},{},{},{},{}",-score(&fit,&test)/n as f64,fit[0].exp(),fit[1].exp(),fit[2],fit[3],killed_states(&fit,4.0).1,killed_states(&actual,4.0).1)?;
                }
            }
        }
    }
    println!("Multisite identifiability: three truths x three observation architectures x eight synthetic replicates");
    Ok(())
}

fn csv_records(path: &str) -> Result<Vec<BTreeMap<String, String>>> {
    let mut lines = reader(path)?.lines();
    let header = fields(&lines.next().ok_or("empty CSV")??, ',');
    lines
        .map(|l| {
            let a = fields(&l?, ',');
            if a.len() != header.len() {
                return Err(format!("wrong CSV width: {path}").into());
            }
            Ok(header.iter().cloned().zip(a).collect())
        })
        .collect()
}
fn extended_report(out: &str) -> Result<()> {
    type Summaries = BTreeMap<String, BTreeMap<String, Vec<f64>>>;
    let mut summaries: Summaries = BTreeMap::new();
    let mut add = |key: String, unit: String, value: f64| {
        if value.is_finite() {
            summaries
                .entry(key)
                .or_default()
                .entry(unit)
                .or_default()
                .push(value);
        }
    };
    for tolerance in [8, 32] {
        for outcome in ["RNA_ECM", "RNA_ECM_absolute", "RNA_global"] {
            let records = csv_records(&format!("{out}/{outcome}_{tolerance}_prediction.csv"))?;
            let mut random: BTreeMap<(String, String), Vec<f64>> = BTreeMap::new();
            for r in &records {
                if r["model"].starts_with("matched_random") {
                    random
                        .entry((r["held_out_donor"].clone(), r["penalty"].clone()))
                        .or_default()
                        .push(r["mse_gain_vs_baseline"].parse()?);
                }
            }
            for r in records
                .iter()
                .filter(|r| !r["model"].starts_with("matched_random"))
            {
                let gain = r["mse_gain_vs_baseline"].parse()?;
                add(
                    format!(
                        "{outcome} tolerance={tolerance} penalty={} {} MSE gain",
                        r["penalty"], r["model"]
                    ),
                    r["held_out_donor"].clone(),
                    gain,
                );
                if r["model"] == "functional_ECM" {
                    add(format!("{outcome} tolerance={tolerance} penalty={} functional minus matched random MSE gain",r["penalty"]),r["held_out_donor"].clone(),gain-mean(&random[&(r["held_out_donor"].clone(),r["penalty"].clone())]));
                }
            }
        }
    }
    let kinetic = csv_records(&format!("{out}/maintenance_holdout.csv"))?;
    let pooled: BTreeMap<_, f64> = kinetic
        .iter()
        .filter(|r| r["model"] == "pooled_constant")
        .map(|r| {
            (
                (
                    r["chrom"].clone(),
                    r["position"].clone(),
                    r["held_out_hour"].clone(),
                ),
                r["log_score_per_read"].parse().unwrap(),
            )
        })
        .collect();
    for r in kinetic.iter().filter(|r| r["model"].starts_with("pulse_")) {
        let key = (
            r["chrom"].clone(),
            r["position"].clone(),
            r["held_out_hour"].clone(),
        );
        let unit = format!("{}:{}", r["chrom"], r["position"].parse::<u64>()? / 100000);
        add(
            format!(
                "kinetics {} hour={} log score gain vs pooled constant",
                r["model"], r["held_out_hour"]
            ),
            unit,
            r["log_score_per_read"].parse::<f64>()? - pooled[&key],
        );
    }
    for quality in [0.95, 0.8] {
        for r in csv_records(&format!("{out}/clone_conditioning_{quality}.csv"))? {
            add(
                format!(
                    "clone quality={quality} class={} log score gain",
                    r["class"]
                ),
                r["mouse"].clone(),
                r["log_score_gain"].parse()?,
            );
            add(
                format!(
                    "clone quality={quality} class={} equal-clone reweighting",
                    r["class"]
                ),
                r["mouse"].clone(),
                r["absolute_reweighting_change"].parse()?,
            );
        }
        for r in csv_records(&format!("{out}/clone_representation_{quality}.csv"))? {
            add(format!("clone representation quality={quality} class={} absolute composition component",r["class"]),r["probe"].clone(),r["symmetric_representation_component"].parse::<f64>()?.abs());
            add(format!("clone representation quality={quality} class={} absolute within-clone component",r["class"]),r["probe"].clone(),r["symmetric_within_clone_component"].parse::<f64>()?.abs());
        }
    }
    let perturb = csv_records(&format!("{out}/perturbation_prediction.csv"))?;
    for r in &perturb {
        add(
            format!(
                "perturbation {} {} log score gain vs unchanged",
                r["mutation"], r["model"]
            ),
            r["replicate"].clone(),
            r["score_gain_vs_unchanged"].parse()?,
        );
        if r["model"] == "unseen_combination_independent_effects" {
            let calibrated = perturb
                .iter()
                .find(|s| {
                    s["mutation"] == r["mutation"]
                        && s["replicate"] == r["replicate"]
                        && s["held_out_chromosome_group"] == r["held_out_chromosome_group"]
                        && s["model"] == "CTMC_effective_loss"
                })
                .unwrap();
            add(
                format!(
                    "perturbation {} independent combination minus calibrated log score",
                    r["mutation"]
                ),
                r["replicate"].clone(),
                r["mean_log_score_per_read"].parse::<f64>()?
                    - calibrated["mean_log_score_per_read"].parse::<f64>()?,
            );
        }
    }
    for r in csv_records(&format!("{out}/trajectory_distribution_checks.csv"))? {
        add(
            format!("trajectory {} RMSE", r["model"]),
            r["donor"].clone(),
            r["RMSE"].parse()?,
        );
        add(
            format!(
                "trajectory {} fraction outside nominal 95 percent",
                r["model"]
            ),
            r["donor"].clone(),
            r["fraction_abs_z_gt_1_96"].parse()?,
        );
    }
    let ident = csv_records(&format!("{out}/multisite_identifiability.csv"))?;
    type CandidateFits<'a> = BTreeMap<(String, String, String), Vec<&'a BTreeMap<String, String>>>;
    let mut fits: CandidateFits<'_> = BTreeMap::new();
    for r in &ident {
        fits.entry((
            r["truth"].clone(),
            r["architecture"].clone(),
            r["seed"].clone(),
        ))
        .or_default()
        .push(r);
    }
    let mut winners = writer(&format!("{out}/identifiability_winners.csv"))?;
    writeln!(winners,"truth,architecture,synthetic_replicate,best_heldout_candidate,heldout_log_score,includes_true_mechanism,competitive_candidates_at_0_001,competitive_candidate_count")?;
    for ((truth, architecture, seed), models) in fits {
        let best = models
            .iter()
            .max_by(|a, b| {
                a["test_log_predictive_density"]
                    .parse::<f64>()
                    .unwrap()
                    .total_cmp(&b["test_log_predictive_density"].parse::<f64>().unwrap())
            })
            .unwrap();
        let includes =
            best["candidate"] == truth || (best["candidate"] == "both" && truth != "both");
        let best_score = best["test_log_predictive_density"].parse::<f64>()?;
        let competitive: Vec<_> = models
            .iter()
            .filter(|r| {
                best_score - r["test_log_predictive_density"].parse::<f64>().unwrap() <= 0.001
            })
            .map(|r| r["candidate"].clone())
            .collect();
        add(
            format!("identifiability {architecture} fraction with multiple competitive candidates"),
            format!("{truth}:{seed}"),
            f64::from(competitive.len() > 1),
        );
        writeln!(
            winners,
            "{truth},{architecture},{seed},{},{},{includes},{},{}",
            best["candidate"],
            best["test_log_predictive_density"],
            competitive.join(";"),
            competitive.len()
        )?;
        add(
            format!("identifiability truth={truth} {architecture} strict highest-score candidate frequency"),
            seed,
            f64::from(best["candidate"] == truth),
        );
    }
    let mut w = writer(&format!("{out}/six_test_summary.csv"))?;
    writeln!(
        w,
        "comparison,units,mean,bootstrap_lower,bootstrap_upper,unit_definition_and_limit"
    )?;
    let mut md = writer(&format!("{out}/six_test_validation.md"))?;
    writeln!(md,"# Six public-data validation tests (2026-10-03)\n\nThe extension finds useful predictive behavior and several failures of the simplest models. It **does not establish functional constraint or biological identity loss**. The functional RNA comparison is negative in the tested settings; kinetic curves do not outperform the pooled constant comparator; future-array intervals are undercalibrated. Clone structure and non-additive perturbation effects require explicit treatment. These results qualify the earlier stability-only evidence rather than negate its observed transfer.\n\nAll metrics and sensitivity settings are retained in six_test_summary.csv and source-specific CSVs. Positive prediction gains mean improvement; negative RNA MSE gains mean worse prediction. The table uses the middle fixed ridge penalty 0.01 for display only.\n\n| Comparison | Units | Mean | 95% resampling interval |\n|---|---:|---:|---:|")?;
    for (key, units) in &summaries {
        let values: Vec<_> = units.values().map(|v| mean(v)).collect();
        let (lo, hi) = interval(&values, 20261003);
        let status = if key.starts_with("kinetics") {
            "100kb genomic blocks; one stem-cell line; not biological-replicate confidence"
        } else if key.starts_with("clone representation") {
            "probes; two recipients; descriptive no biological confidence"
        } else if key.starts_with("clone") {
            "two recipients; published classes selected in same cohort"
        } else if key.starts_with("identifiability") {
            "synthetic replicates; known baseline hazard; not biological validation"
        } else {
            "donors or culture replicates; fixed-fold bootstrap omits refitting uncertainty"
        };
        writeln!(
            w,
            "{key},{},{},{lo},{hi},{status}",
            values.len(),
            mean(&values)
        )?;
        let display = key.contains("penalty=0.01 functional_ECM")
            || key.contains("penalty=0.01 functional minus matched")
            || key.contains("pulse_uniform")
            || key.contains("quality=0.95 class=Static log")
            || key.contains("independent combination minus")
            || key.starts_with("trajectory")
            || key.contains("strict highest-score candidate frequency")
            || key.contains("multiple competitive candidates");
        if display {
            writeln!(
                md,
                "| {key} | {} | {:.6} | [{lo:.6}, {hi:.6}] |",
                values.len(),
                mean(&values)
            )?;
        }
    }
    writeln!(md,"\n## 1. Independent functional outcomes\n\nA frozen external Reactome extracellular-matrix pathway (R-HSA-1474244) supplies 307 annotated genes and 2,698 supported promoter probes; RNA includes 321 pathway genes. Probe weights balance genes. Twenty control masks sample other gene promoters with replacement in the same chromosome/CGI/probe-density/sequence strata, carrying the same weights. No methylation or RNA outcome selects the mask. This is a functional-program proxy, not essentiality or an identity-loss boundary.\n\nControl fibroblast DNA at time t forecasts RNA at t+1..16 days. Both assays use their independent earlier culture references; no RNA outcome or held-out donor selects methylation weights. Starting expression, earlier beta, culture time, forecast horizon, SURF1 status and oxygen are baseline covariates. The primary baseline-match tolerance is 8 days (33 pairs/four donors); the exploratory 32-day sensitivity retains 63 pairs/seven donors and still enforces RNA baseline earlier than prediction time. HC5/HC6 lack corresponding RNA culture groups. Outcomes are signed mean ECM change, absolute ECM transcript displacement and absolute displacement of an input-order transcript subset. The absolute ECM outcome and wider tolerance were added during the audit, so are exploratory. Most functional-weight gains are negative; matched masks also often fail. This does not support the proposed functional-distance prediction in this cohort/definition. It does not disprove all functional weights. Age/sex/batch/proliferation adjustment and RNA future-increment forecasting remain incomplete.\n\n## 2. Direct post-replication kinetics\n\nThe corrected author commit ec46d95b68fad06072702686c248cad55754fc3b provides observed methylated/unmethylated counts at post-pulse 0/1/4/16 hours. Fixed every-32nd-site sampling precedes coverage filtering (at least ten reads at zero and five at each later time). There are 2,959 covered CpGs on chr1 and chr22. Each of 4h and 16h is excluded in turn. We fit f*(1-exp(-k*t)) with either half-hour timing correction or exact one-hour pulse averaging; this marginal is mathematically CTMC-equivalent but cannot separate a mixed plateau from reversible maintenance. A training-pooled constant has better aggregate held-out score than these kinetics. Thus this test gives no predictive advantage for the tested time dependence. Many rates reach boundaries and only one biological cell line is represented. Hour-scale nascent-strand rates are not the year/day-scale drift rates inferred from aging arrays.\n\n## 3. Clone and representation controls\n\nThe public EPI-Clone M.1–M.3 object has 28,782 cells. Tests use independent LARRY barcodes and raw antibody UMI counts from the two stained main-experiment recipients. Zero detected antibody molecules are retained only in nonempty stained libraries. Six surface-marker covariates plus antibody-library size and undigested-control quality adjust cell state/measurement. Probe amplification in independent undigested controls must be at least 95%. Cell quality thresholds 0.95 and 0.8 retain 1,903 and 3,285 cells respectively; 370 probes pass assay controls. Methylation-sensitive-digestion amplification is a proxy with dropout, not clean allele-specific bisulfite calls.\n\nHeld-out cells within independently labelled clones benefit from a shrunken clone residual, especially at published static CpGs. Published static/dynamic classes use this cohort's proteins and are descriptive, not independent feature validation. Equal-clone and cell-weighted averages use the same qualifying clone universe. A symmetric exact decomposition across shared clones divides mouse3-to-mouse4 bulk differences into within-clone changes and clone-frequency changes; these are **different recipients**, not a longitudinal survival experiment. Clonal representation matters, but selection causing aging stability is not identified. Human raw reads are controlled access and were not used.\n\n## 4. Existing DNMT1 perturbations\n\nGSE145698 contains three WT cultures and three cultures for each of five DNMT1 mutants. Fixed every-32nd-row sampling and WT coverage yield 165,106 autosomal CpGs. WT counts define initial beta; mutant training chromosomes calibrate an effective CTMC loss exposure, then opposite-parity chromosomes are excluded for evaluation. An affine empirical comparator is virtually as predictive: calibrated scaling does not validate a unique mechanism. No elapsed culture exposure is given, so the fitted parameter is a rate-times-exposure product.\n\nAn additional genotype-held-out test predicts W465A/W796A from the fitted separate W465A and W796A exposures, and the triple mutant from W464A/W465A plus W796A. The simple independent-effects combination substantially underperforms mutation-specific calibration for W465A/W796A, consistent with domain interactions. This rejects that transfer rule, not every CTMC or the general drift framework. No mutation-specific data calibrate the held-out combination rule.\n\n## 5. Full predictive distribution\n\nControl fibroblast trajectories (fixed every-256th probe) fit their earlier three quarters and predict strictly later arrays. 54,083 covered CpG/culture trajectories produce 50 future array summaries. Independent CTMC means improve RMSE over initial-constant and linear forecasts, but Gaussian predictive intervals miss far more than the nominal 5%. Training residual variance plus a 0.02 array-noise floor omits parameter uncertainty and changing variance; this is a failure of the current predictive distribution, not proof that intrinsic CTMC dynamics are impossible. Covariance analysis uses a fixed 128-probe subset and independent within-donor residual permutations: shared residual structure remains. Batch, culture, oxygen and common proliferation changes can generate it. This is an exploratory rejection of unconditional independence, not a validated biological M5 rank. No threshold is called identity loss.\n\n## 6. What observations identify\n\nAn exact eight-state killed chain retains three loci and shared state-dependent killing. Three generating mechanisms, three observation architectures and eight synthetic repeats are compared with independent/protection/selection/both candidates and independent test simulations. There are 400 sampled survivors/timepoint, 40% missing locus calls and 2% call error. Bulk endpoints expose only locus marginals; repeated destructive cell snapshots expose joint patterns; adding independent survival counts supplies 500 at-risk observations/timepoint. Snapshot-only candidates often confuse protection and selection; survival helps recover selection and combined mechanisms. Strict highest-score winner frequencies are sensitive to nearly tied scores; identifiability_winners.csv also lists every candidate within 0.001 mean log score of the best. This tolerance was added during the numerical audit, not preregistered. Multiple candidates remain competitive in 22/24 bulk-endpoint runs, 21/24 repeated-cell runs and 20/24 survival-augmented runs. Survival improves some strict rankings but does not demonstrate robust mechanism identification at this sampling depth. Nested candidates can overfit pure mechanisms even with survival. Baseline killing (0.05), targets, error and missingness are known and no clone reproduction/site-context confounding is simulated. This measures conditional recoverability, not identifiability for unrestricted biology.\n\n## Reproduction and remaining gates\n\nAll source files live under ignored data/raw and are checksummed. A pinned Python helper only decodes MATLAB/R objects to TSV; all inference, matching, predictions, simulations and reports are Rust. `make validation` converts formats, rebuilds functional/trajectory summaries, runs all six tests and regenerates this report. Feature assembly is hg19 for human probes; mouse assays are never joined to human coordinates.\n\nIndependent functional weights beyond ECM promoters, RNA/identity thresholds, family-aware blood inference, clone reproduction and survival likelihoods for real aging cohorts, full observation uncertainty and intervention-specific molecular recruitment still require work. No broad causal biological-validation claim follows from these tests.")?;
    Ok(())
}

fn recovery_detected_probes(path: &str) -> Result<BTreeSet<String>> {
    let mut is_baseline = false;
    let mut table = false;
    let mut header = true;
    let mut samples: Vec<BTreeSet<String>> = vec![];
    for line in reader(path)?.lines() {
        let line = line?;
        if line.starts_with("^SAMPLE") {
            is_baseline = false;
            table = false;
        }
        if let Some(title) = line.strip_prefix("!Sample_title = ") {
            is_baseline =
                title.contains("Ctrl") && !title.contains("1KO") && !title.contains("3BKO");
        }
        if line == "!sample_table_begin" {
            table = is_baseline;
            header = true;
            if table {
                samples.push(BTreeSet::new());
            }
            continue;
        }
        if line == "!sample_table_end" {
            table = false;
            continue;
        }
        if !table {
            continue;
        }
        let a = fields(&line, '\t');
        if header {
            if a.len() != 3 || a[2] != "Detection Pval" {
                return Err("unexpected expression detection format".into());
            }
            header = false;
            continue;
        }
        if a.len() != 3 {
            return Err("bad baseline detection row".into());
        }
        if numeric(&a[2]).is_some_and(|p| (0.0..=0.01).contains(&p)) {
            samples.last_mut().unwrap().insert(a[0].clone());
        }
    }
    if samples.len() != 2 {
        return Err("expected two WT baseline detection tables".into());
    }
    Ok(samples[0].intersection(&samples[1]).cloned().collect())
}

// A narrow functional-constraint pilot; matching never reads recovery outcomes.
struct RecoverySite {
    id: String,
    gene: String,
    essential: bool,
    feature: Feature,
    beta: Vec<f64>,
    expression: Vec<f64>,
}
type RecoveryMatrix = (Vec<String>, BTreeMap<String, Vec<f64>>);
fn recovery_matrix(path: &str) -> Result<RecoveryMatrix> {
    let mut titles = vec![];
    let mut table = false;
    let mut rows = BTreeMap::new();
    for line in reader(path)?.lines() {
        let line = line?;
        if line.starts_with("!Sample_title\t") {
            titles = fields(&line, '\t').into_iter().skip(1).collect();
        } else if line == "!series_matrix_table_begin" {
            table = true;
        } else if line == "!series_matrix_table_end" {
            break;
        } else if table && !line.starts_with("\"ID_REF\"") {
            let a = fields(&line, '\t');
            if a.len() != titles.len() + 1 {
                return Err("recovery matrix/sample count mismatch".into());
            }
            rows.insert(
                a[0].clone(),
                a[1..]
                    .iter()
                    .map(|x| numeric(x).unwrap_or(f64::NAN))
                    .collect(),
            );
        }
    }
    if rows.is_empty() {
        return Err("empty recovery matrix".into());
    }
    Ok((titles, rows))
}
fn recovery_annotation(
    path: &str,
    gene_column: &str,
    body_only: bool,
) -> Result<BTreeMap<String, String>> {
    let mut table = false;
    let mut header = vec![];
    let mut map = BTreeMap::new();
    for line in reader(path)?.lines() {
        let line = line?;
        if line == "!platform_table_begin" {
            table = true;
            continue;
        }
        if line == "!platform_table_end" {
            break;
        }
        if !table {
            continue;
        }
        let a = fields(&line, '\t');
        if header.is_empty() {
            header = a;
            continue;
        }
        let names = header
            .iter()
            .position(|x| x == gene_column)
            .ok_or("missing recovery gene column")?;
        let genes: BTreeSet<_> = a[names].split(';').filter(|g| !g.is_empty()).collect();
        if genes.len() != 1 {
            continue;
        }
        if body_only {
            let groups = header
                .iter()
                .position(|x| x == "UCSC_RefGene_Group")
                .ok_or("missing BODY annotation")?;
            if !a[groups].split(';').all(|g| g == "Body") {
                continue;
            }
        }
        map.insert(a[0].clone(), genes.into_iter().next().unwrap().to_string());
    }
    if map.is_empty() {
        return Err("no recovery annotations".into());
    }
    Ok(map)
}
fn recovery_indices(titles: &[String], genotype: &str, days: &[u32]) -> Result<Vec<Vec<usize>>> {
    days.iter()
        .map(|day| {
            let found: Vec<_> = titles
                .iter()
                .enumerate()
                .filter(|(_, title)| {
                    let is_genotype = match genotype {
                        "WT" => !title.contains("1KO") && !title.contains("3BKO"),
                        other => title.contains(other),
                    };
                    let sample_day = if title.contains("Ctrl") {
                        Some(0)
                    } else {
                        title.split_whitespace().find_map(|part| {
                            part.strip_prefix('D')
                                .and_then(|d| d.split('_').next()?.parse::<u32>().ok())
                        })
                    };
                    is_genotype && sample_day == Some(*day)
                })
                .map(|(i, _)| i)
                .collect();
            if found.is_empty() {
                Err(format!("missing {genotype} day {day}").into())
            } else {
                Ok(found)
            }
        })
        .collect()
}
fn recovery_fraction(beta: &[f64], index: usize) -> f64 {
    (beta[index] - beta[1]) / (beta[0] - beta[1])
}
fn recovery_match(a: &RecoverySite, b: &RecoverySite, scales: [f64; 5]) -> Option<f64> {
    if a.feature.chrom != b.feature.chrom || a.feature.island != b.feature.island {
        return None;
    }
    let differences = [
        (a.beta[0] - b.beta[0]).abs() / 0.1,
        ((a.beta[0] - a.beta[1]) - (b.beta[0] - b.beta[1])).abs() / 0.1,
        (a.expression[0] - b.expression[0]).abs() / 2.0,
        (a.feature.sequence? - b.feature.sequence?).abs() / 0.25,
        (a.feature.density - b.feature.density).abs() / 0.02,
    ];
    differences
        .iter()
        .zip(scales)
        .all(|(d, scale)| *d <= scale)
        .then(|| differences.iter().map(|d| d * d).sum())
}
fn recovery_pilot(root: &str, out: &str, detected_only: bool) -> Result<()> {
    let prefix = if detected_only {
        "recovery_detected"
    } else {
        "recovery"
    };
    let reference = |name: &str| -> Result<BTreeSet<String>> {
        reader(&format!("{root}/{name}"))?
            .lines()
            .skip(1)
            .map(|l| Ok(l?.split('\t').next().unwrap_or("").to_string()))
            .collect()
    };
    let essential = reference("BAGEL_CEGv2.txt")?;
    let nonessential = reference("BAGEL_NEGv1.txt")?;
    if !essential.is_disjoint(&nonessential) {
        return Err("functional reference overlap".into());
    }
    let features = load_features()?;
    let genes = recovery_annotation(
        &format!("{root}/GSE73115_family.soft.gz"),
        "UCSC_RefGene_Name",
        true,
    )?;
    let rna_genes =
        recovery_annotation(&format!("{root}/GSE51811_family.soft.gz"), "Symbol", false)?;
    let (rna_titles, rna) = recovery_matrix(&format!("{root}/GSE51811_series_matrix.txt.gz"))?;
    let rna_indices = recovery_indices(&rna_titles, "WT", &[0, 5, 14, 24, 42])?;
    let detected = recovery_detected_probes(&format!("{root}/GSE51811_family.soft.gz"))?;
    let mut expression: BTreeMap<String, Vec<Vec<f64>>> = BTreeMap::new();
    for (id, values) in &rna {
        if detected_only && !detected.contains(id) {
            continue;
        }
        let Some(gene) = rna_genes.get(id) else {
            continue;
        };
        if !essential.contains(gene) && !nonessential.contains(gene) {
            continue;
        }
        if rna_indices[0]
            .iter()
            .any(|i| !values[*i].is_finite() || values[*i] < 0.0)
        {
            continue;
        }
        let x: Vec<_> = rna_indices
            .iter()
            .map(|indices| {
                mean(
                    &indices
                        .iter()
                        .map(|i| (values[*i] + 1.0).log2())
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        expression.entry(gene.clone()).or_default().push(x);
    }
    let expression: BTreeMap<_, Vec<_>> = expression
        .into_iter()
        .map(|(gene, rows)| {
            let average = (0..5)
                .map(|i| mean(&rows.iter().map(|r| r[i]).collect::<Vec<_>>()))
                .collect();
            (gene, average)
        })
        .collect();
    let (titles, matrix) = recovery_matrix(&format!("{root}/GSE51810_series_matrix.txt.gz"))?;
    let days = [0, 5, 14, 24, 42, 54, 68];
    let indices = recovery_indices(&titles, "WT", &days)?;
    if indices.iter().any(|i| i.len() != 1) {
        return Err("expected one methylation array per WT time".into());
    }
    let mutant_indices = recovery_indices(&titles, "3BKO", &days)?;
    let mut sites = vec![];
    for (id, values) in &matrix {
        let (Some(gene), Some(feature)) = (genes.get(id), features.get(id)) else {
            continue;
        };
        let Some(exp) = expression.get(gene) else {
            continue;
        };
        if feature.sequence.is_none() || feature.island == "Unknown" || feature.island.is_empty() {
            continue;
        }
        let beta: Vec<_> = indices.iter().map(|i| values[i[0]]).collect();
        // Eligibility uses only baseline and day 5; later finite values are checked after matching.
        if !beta[0].is_finite() || !beta[1].is_finite() || beta[0] < 0.7 || beta[0] - beta[1] < 0.15
        {
            continue;
        }
        sites.push(RecoverySite {
            id: id.clone(),
            gene: gene.clone(),
            essential: essential.contains(gene),
            feature: feature.clone(),
            beta,
            expression: exp.clone(),
        });
    }
    let mut q_genes: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, site) in sites.iter().enumerate().filter(|(_, s)| s.essential) {
        q_genes.entry(site.gene.clone()).or_default().push(i);
    }
    let controls: Vec<_> = sites
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.essential)
        .map(|(i, _)| i)
        .collect();
    let mut used_controls = BTreeSet::new();
    let mut pairs = vec![];
    for q_sites in q_genes.values() {
        let best = q_sites
            .iter()
            .flat_map(|i| {
                controls
                    .iter()
                    .filter(|j| !used_controls.contains(&sites[**j].gene))
                    .filter_map(|j| {
                        recovery_match(&sites[*i], &sites[*j], [1.0; 5]).map(|d| (*i, *j, d))
                    })
            })
            .min_by(|a, b| a.2.total_cmp(&b.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
        if let Some((i, j, d)) = best {
            used_controls.insert(sites[j].gene.clone());
            pairs.push((i, j, d));
        }
    }
    // Added after observing the primary result: exploratory confounding/support audits.
    let mut sensitivity = writer(&format!("{out}/{prefix}_matching_sensitivity.csv"))?;
    writeln!(sensitivity,"setting,pairs,mean_day42_contrast,mean_log2_expression_imbalance,mean_baseline_beta_imbalance,mean_induced_loss_imbalance")?;
    let mut audit_manifest = writer(&format!("{out}/{prefix}_sensitivity_loci.csv"))?;
    writeln!(
        audit_manifest,
        "setting,essential_gene,essential_cpg,control_gene,control_cpg,match_distance"
    )?;
    for (setting, scales, reverse, expression_floor) in [
        ("primary", [1.0; 5], false, 0.0),
        ("reverse_gene_order", [1.0; 5], true, 0.0),
        (
            "expression_half_log2",
            [1.0, 1.0, 0.25, 1.0, 1.0],
            false,
            0.0,
        ),
        ("all_calipers_halved", [0.5; 5], false, 0.0),
        ("both_baseline_expression_at_least_8", [1.0; 5], false, 8.0),
    ] {
        let mut used = BTreeSet::new();
        let mut audit_pairs = vec![];
        let ordered: Vec<_> = if reverse {
            q_genes.values().rev().collect()
        } else {
            q_genes.values().collect()
        };
        for q_sites in ordered {
            let best = q_sites
                .iter()
                .filter(|i| sites[**i].expression[0] >= expression_floor)
                .flat_map(|i| {
                    controls
                        .iter()
                        .filter(|j| {
                            !used.contains(&sites[**j].gene)
                                && sites[**j].expression[0] >= expression_floor
                        })
                        .filter_map(|j| {
                            recovery_match(&sites[*i], &sites[*j], scales).map(|d| (*i, *j, d))
                        })
                })
                .min_by(|a, b| a.2.total_cmp(&b.2).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
            if let Some((i, j, d)) = best {
                used.insert(sites[j].gene.clone());
                audit_pairs.push((i, j, d));
            }
        }
        for (i, j, d) in &audit_pairs {
            writeln!(
                audit_manifest,
                "{setting},{},{},{},{},{d}",
                sites[*i].gene, sites[*i].id, sites[*j].gene, sites[*j].id
            )?;
        }
        let diffs: Vec<_> = audit_pairs
            .iter()
            .map(|(i, j, _)| {
                recovery_fraction(&sites[*i].beta, 4) - recovery_fraction(&sites[*j].beta, 4)
            })
            .filter(|v| v.is_finite())
            .collect();
        let exp: Vec<_> = audit_pairs
            .iter()
            .map(|(i, j, _)| sites[*i].expression[0] - sites[*j].expression[0])
            .collect();
        let baseline: Vec<_> = audit_pairs
            .iter()
            .map(|(i, j, _)| sites[*i].beta[0] - sites[*j].beta[0])
            .collect();
        let loss: Vec<_> = audit_pairs
            .iter()
            .map(|(i, j, _)| {
                (sites[*i].beta[0] - sites[*i].beta[1]) - (sites[*j].beta[0] - sites[*j].beta[1])
            })
            .collect();
        writeln!(
            sensitivity,
            "{setting},{},{},{},{},{}",
            diffs.len(),
            mean(&diffs),
            mean(&exp),
            mean(&baseline),
            mean(&loss)
        )?;
    }
    let mut manifest = writer(&format!("{out}/{prefix}_matched_loci.csv"))?;
    writeln!(manifest,"pair,class,gene,cpg,chrom,cgi,baseline_beta,day5_beta,baseline_log2_expression,sequence_rank,density,match_distance,eligible_for_primary")?;
    let mut observations = writer(&format!("{out}/{prefix}_pair_outcomes.csv"))?;
    writeln!(
        observations,
        "pair,essential_gene,control_gene,day,essential_recovery,control_recovery,contrast"
    )?;
    let mut summary = writer(&format!("{out}/{prefix}_pilot_summary.csv"))?;
    writeln!(summary,"endpoint,units,essential_mean,control_mean,contrast,heterogeneity_bootstrap_low,heterogeneity_bootstrap_high")?;
    for (pair, (i, j, d)) in pairs.iter().enumerate() {
        for site in [&sites[*i], &sites[*j]] {
            writeln!(
                manifest,
                "{pair},{},{},{},{},{},{},{},{},{},{},{},{}",
                if site.essential {
                    "essential"
                } else {
                    "nonessential"
                },
                site.gene,
                site.id,
                site.feature.chrom,
                site.feature.island,
                site.beta[0],
                site.beta[1],
                site.expression[0],
                site.feature.sequence.unwrap(),
                site.feature.density,
                d,
                sites[*i].beta[4].is_finite() && sites[*j].beta[4].is_finite()
            )?;
        }
    }
    let mut contrasts = vec![];
    for (index, day) in days.iter().enumerate().skip(2) {
        let mut q = vec![];
        let mut c = vec![];
        for (pair, (i, j, _)) in pairs.iter().enumerate() {
            let a = recovery_fraction(&sites[*i].beta, index);
            let b = recovery_fraction(&sites[*j].beta, index);
            if !a.is_finite() || !b.is_finite() {
                continue;
            }
            q.push(a);
            c.push(b);
            writeln!(
                observations,
                "{pair},{},{},{day},{a},{b},{}",
                sites[*i].gene,
                sites[*j].gene,
                a - b
            )?;
        }
        let differences: Vec<_> = q.iter().zip(&c).map(|(a, b)| a - b).collect();
        let (lo, hi) = if differences.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            interval(&differences, 20261003)
        };
        writeln!(
            summary,
            "DNA_day{day},{},{},{},{},{lo},{hi}",
            q.len(),
            mean(&q),
            mean(&c),
            mean(&differences)
        )?;
        if *day == 42 {
            contrasts = differences;
        }
    }
    // Expression is an independently measured secondary readout; it never selects loci.
    let mut q = vec![];
    let mut c = vec![];
    for (i, j, _) in &pairs {
        let a = &sites[*i].expression;
        let b = &sites[*j].expression;
        if a.iter().chain(b).any(|v| !v.is_finite())
            || (a[0] - a[1]).abs() < 0.25
            || (b[0] - b[1]).abs() < 0.25
        {
            continue;
        }
        q.push(((a[0] - a[1]).abs() - (a[0] - a[4]).abs()) / (a[0] - a[1]).abs());
        c.push(((b[0] - b[1]).abs() - (b[0] - b[4]).abs()) / (b[0] - b[1]).abs());
    }
    let differences: Vec<_> = q.iter().zip(&c).map(|(a, b)| a - b).collect();
    let (lo, hi) = if differences.is_empty() {
        (f64::NAN, f64::NAN)
    } else {
        interval(&differences, 20261003)
    };
    writeln!(
        summary,
        "RNA_day42,{},{},{},{},{lo},{hi}",
        q.len(),
        mean(&q),
        mean(&c),
        mean(&differences)
    )?;
    // Mutation-specific baseline/loss normalization: descriptive positive control, no causal identification.
    for (index, day) in days.iter().enumerate().skip(2) {
        let mut wt = vec![];
        let mut mutant = vec![];
        for (i, j, _) in &pairs {
            for site in [&sites[*i], &sites[*j]] {
                let values = &matrix[&site.id];
                let beta: Vec<_> = mutant_indices.iter().map(|ix| values[ix[0]]).collect();
                if beta[0] < 0.7 || beta[0] - beta[1] < 0.15 {
                    continue;
                }
                let (a, b) = (
                    recovery_fraction(&site.beta, index),
                    recovery_fraction(&beta, index),
                );
                if a.is_finite() && b.is_finite() {
                    wt.push(a);
                    mutant.push(b);
                }
            }
        }
        let d: Vec<_> = wt.iter().zip(&mutant).map(|(a, b)| a - b).collect();
        let (lo, hi) = if d.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            interval(&d, 20261003)
        };
        writeln!(
            summary,
            "WT_minus_3BKO_day{day},{},{},{},{},{lo},{hi}",
            wt.len(),
            mean(&wt),
            mean(&mutant),
            mean(&d)
        )?;
    }
    let q_sites = sites.iter().filter(|s| s.essential).count();
    let c_genes: BTreeSet<_> = controls.iter().map(|i| &sites[*i].gene).collect();
    fs::write(format!("{out}/{prefix}_data_audit.md"),format!("# Recovery pilot data audit\n\nGSE51810: {} methylation arrays, one WT and one DNMT3B-KO array at each requested time. GSE51811: {} expression arrays; two WT arrays per time, averaged on log2(value+1) scale. Sample titles define genotype/time because GEO characteristics contain inconsistent labels. BODY-only/single-gene mapping uses the frozen hg19 GPL13534 table; expression mapping uses GPL10558 Symbol with no inferred aliases. No methylation detection-p filtering is possible from these series matrices. Baseline RNA detection filtering in this run: {detected_only}; detected-only sensitivity requires p<=0.01 in both WT baseline arrays and never uses later detection scores. Negative processed RNA intensities produce unavailable log2 measurements rather than being clipped. Gene essentiality does not establish the importance of its BODY methylation.\n\nExternal reference: {} core-essential and {} nonessential symbols. Eligible after baseline/loss/context/expression filters: {q_sites} essential probes/{} genes and {} nonessential probes/{} genes. Greedy deterministic nearest matching without reuse of either gene produces {} pairs; {} have finite primary outcomes. One probe pair per distinct gene pair balances genes. Matching fixes chromosome/CGI and calipers baseline beta 0.1, induced loss 0.1, baseline expression 2 log2 units, sequence rank 0.25 and CpG density 0.02. Alphabetical essential-gene order can affect support. Matching never uses later DNA/RNA outcomes.\n\nPrimary mean essential-minus-control day42 recovery: {}. Intervals in the CSV resample matched loci and describe heterogeneity; they do not estimate biological reproducibility. The primary interpretation is feasibility only.\n",titles.len(),rna_titles.len(),essential.len(),nonessential.len(),q_genes.len(),controls.len(),c_genes.len(),pairs.len(),contrasts.len(),mean(&contrasts)))?;
    println!(
        "Recovery pilot detected_only={detected_only}: {} gene pairs; day42 essential-minus-control recovery {}",
        pairs.len(),
        mean(&contrasts)
    );
    Ok(())
}

#[cfg(test)]
mod validation_tests {
    use super::*;
    #[test]
    fn recovery_titles_separate_genotypes_and_day_numbers() {
        let titles: Vec<_> = [
            "HCT116 Ctrl (2)",
            "HCT116 D5_1",
            "HCT116 D14_2",
            "HCT116 1KO Ctrl",
            "HCT116 3BKO Ctrl_1",
            "HCT116 3BKO  5-Aza-CdR D10_1",
            "HCT116 3BKO  5-Aza-CdR D14_2",
        ]
        .into_iter()
        .map(String::from)
        .collect();
        assert_eq!(
            recovery_indices(&titles, "WT", &[0, 5, 14]).unwrap(),
            vec![vec![0], vec![1], vec![2]]
        );
        assert_eq!(
            recovery_indices(&titles, "3BKO", &[0, 10, 14]).unwrap(),
            vec![vec![4], vec![5], vec![6]]
        );
        assert!(recovery_indices(&titles, "WT", &[10]).is_err());
        // Overshoot and continued loss are observable outcomes, not clipped recovery.
        assert!((recovery_fraction(&[0.8, 0.4, 1.0], 2) - 1.5).abs() < 1e-12);
        assert!((recovery_fraction(&[0.8, 0.4, 0.2], 2) + 0.5).abs() < 1e-12);
    }
    #[test]
    fn exact_multisite_killing_reduces_to_independent_chains() {
        let time = 3.0;
        let parameters = [0.15_f64.ln(), 0.3_f64.ln(), 0.0, 0.0];
        let (states, survival) = killed_states(&parameters, time);
        let beta = probability(0.0, 0.15, 0.3, time);
        assert!((survival - (-0.05 * time).exp()).abs() < 1e-10);
        for (state, p) in states.iter().enumerate() {
            let n = state.count_ones() as i32;
            assert!((p - beta.powi(n) * (1.0 - beta).powi(3 - n)).abs() < 1e-10);
        }
    }
    #[test]
    fn pulse_average_matches_numerical_integration_and_call_patterns_normalize() {
        for rate in [0.001, 0.2, 4.0, 50.0] {
            let integrated = (0..10000)
                .map(|i| 0.8 * (1.0 - (-rate * (2.0 + (i as f64 + 0.5) / 10000.0)).exp()))
                .sum::<f64>()
                / 10000.0;
            assert!((pulse_probability(rate, 0.8, 2.0, true) - integrated).abs() < 1e-7);
        }
        let (states, _) = killed_states(&[0.15_f64.ln(), 0.3_f64.ln(), 1.0, 2.0], 4.0);
        for mask in 1..8 {
            let sum: f64 = (0..8)
                .filter(|p| p & !mask == 0)
                .map(|p| observed_pattern(&states, p, mask))
                .sum();
            assert!((sum - 1.0).abs() < 1e-12);
        }
    }
}
