use epidrift::*;
use flate2::read::MultiGzDecoder;
use std::{
    collections::{BTreeMap, BTreeSet},
    error::Error,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};
type Result<T> = std::result::Result<T, Box<dyn Error>>;
fn reader(p: &str) -> Result<Box<dyn BufRead>> {
    let f = File::open(p)?;
    Ok(if p.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    })
}
// Input matrices use quoted fields but never embedded newlines. Handle escaped quotes.
fn fields(s: &str, sep: char) -> Vec<String> {
    let mut out = vec![];
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = s.trim_end_matches(['\r', '\n']).chars().peekable();
    while let Some(c) = chars.next() {
        if c == '"' {
            if quoted && chars.peek() == Some(&'"') {
                field.push('"');
                chars.next();
            } else {
                quoted = !quoted;
            }
        } else if c == sep && !quoted {
            out.push(std::mem::take(&mut field));
        } else {
            field.push(c);
        }
    }
    out.push(field);
    out
}
fn numeric(s: &str) -> Option<f64> {
    s.parse::<f64>().ok().filter(|x| x.is_finite())
}
fn writer(p: &str) -> Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(p)?))
}
#[derive(Clone)]
struct Cell {
    id: String,
    mouse: String,
    kind: String,
    age: f64,
    beta: f64,
}
fn cells(root: &str, out: &str) -> Result<Vec<Cell>> {
    let mut cells = vec![];
    for l in reader(&format!("{root}/polycomb_cells.txt"))?.lines() {
        let l = l?;
        let a = fields(&l, '\t');
        if a.len() != 4 {
            return Err(format!("bad cell row: {l}").into());
        }
        let mouse = a[0]
            .split('_')
            .find(|s| {
                s.starts_with("Blood") || (s.starts_with('B') && s[1..].parse::<u32>().is_ok())
            })
            .ok_or("cell without mouse ID")?
            .trim_start_matches("Blood")
            .trim_start_matches('B')
            .parse::<u32>()?;
        cells.push(Cell {
            id: a[0].clone(),
            mouse: format!("mouse{mouse}"),
            kind: a[1].clone(),
            age: a[3].parse()?,
            beta: a[2].parse()?,
        });
    }
    let mut w = writer(&format!("{out}/polycomb_age_summary.csv"))?;
    writeln!(
        w,
        "group,age_weeks,cells,mice,mean_beta,sd_beta,skewness,tail_above_0_089"
    )?;
    for group in ["all", "slow", "fast"] {
        for age in [10.0, 36.0, 77.0, 100.0] {
            let rows: Vec<_> = cells
                .iter()
                .filter(|c| c.age == age && matches_group(c, group))
                .collect();
            let y: Vec<f64> = rows.iter().map(|c| c.beta).collect();
            let m = mean(&y);
            let v = variance(&y);
            let skew =
                y.iter().map(|x| (x - m).powi(3)).sum::<f64>() / y.len() as f64 / v.powf(1.5);
            let mice: BTreeSet<_> = rows.iter().map(|c| &c.mouse).collect();
            writeln!(
                w,
                "{group},{age},{},{},{m},{},{skew},{}",
                y.len(),
                mice.len(),
                v.sqrt(),
                y.iter().filter(|v| **v > 0.089).count() as f64 / y.len() as f64
            )?;
        }
    }
    let mut w = writer(&format!("{out}/cell_manifest.csv"))?;
    writeln!(w, "cell,mouse,cell_type,age_weeks,beta")?;
    for c in &cells {
        writeln!(w, "{},{},{},{},{}", c.id, c.mouse, c.kind, c.age, c.beta)?;
    }
    Ok(cells)
}
fn matches_group(c: &Cell, g: &str) -> bool {
    match g {
        "all" => true,
        "slow" => ["BCell", "NveCD4T", "NveCd8T", "Mono"].contains(&c.kind.as_str()),
        _ => ["EffCD4T", "MemCD8T", "RegT", "NK"].contains(&c.kind.as_str()),
    }
}
// Five-point Gauss-Hermite integration of a shared mouse intercept.
// Scoring is joint within mouse; uncertainty is reported across mice.
fn mouse_logpdf(rows: &[&Cell], initial: f64, iv: f64, tau: f64, p: &[f64]) -> f64 {
    let nodes = [
        -2.8569700138728056,
        -1.355626179974266,
        0.0,
        1.355626179974266,
        2.8569700138728056,
    ];
    let weights: [f64; 5] = [
        0.0112574113277207,
        0.2220759220056126,
        0.5333333333333333,
        0.2220759220056126,
        0.0112574113277207,
    ];
    let terms: Vec<_> = nodes
        .iter()
        .zip(weights)
        .map(|(z, w)| {
            w.ln()
                + rows
                    .iter()
                    .map(|c| {
                        polycomb_logpdf(c.beta, (c.age - 10.0).max(0.0), initial + tau * z, iv, p)
                    })
                    .sum::<f64>()
        })
        .collect();
    logsumexp(&terms)
}
// Independent dense quadrature check of held-out mouse scoring. Fits still use
// the small quadrature; disagreement is an explicit numerical limitation.
fn mouse_logpdf_dense(rows: &[&Cell], initial: f64, iv: f64, tau: f64, p: &[f64]) -> f64 {
    if tau == 0.0 {
        return rows
            .iter()
            .map(|c| polycomb_logpdf(c.beta, (c.age - 10.0).max(0.0), initial, iv, p))
            .sum();
    }
    let weights: Vec<_> = (0..=80)
        .map(|j| {
            let z = j as f64 / 8.0 - 5.0;
            (-0.5 * z * z).exp() * if j == 0 || j == 80 { 0.5 } else { 1.0 }
        })
        .collect();
    let norm = weights.iter().sum::<f64>();
    let terms: Vec<_> = (0..=80)
        .map(|j| {
            let z = j as f64 / 8.0 - 5.0;
            (weights[j] / norm).ln()
                + rows
                    .iter()
                    .map(|c| {
                        polycomb_logpdf(c.beta, (c.age - 10.0).max(0.0), initial + tau * z, iv, p)
                    })
                    .sum::<f64>()
        })
        .collect();
    logsumexp(&terms)
}
fn baseline(cells: &[Cell], out: &str) -> Result<()> {
    let mut w = writer(&format!("{out}/model_comparison.csv"))?;
    writeln!(
        w,
        "dataset,group,held_out_replicate,model,n,log_predictive_density,parameters,status"
    )?;
    let mut ppc = writer(&format!("{out}/posterior_predictive_checks/polycomb.csv"))?;
    writeln!(ppc,"group,mouse,age,observed_mean,predicted_mean,observed_sd,predicted_sd,observed_tail,predicted_tail,wasserstein")?;
    let mut sensitivity = writer(&format!("{out}/random_effect_sensitivity.csv"))?;
    writeln!(
        sensitivity,
        "group,mouse,model,n,tau,fit_quadrature_lpd,dense_evaluation_lpd"
    )?;
    for group in ["slow", "fast"] {
        let mice: BTreeSet<_> = cells
            .iter()
            .filter(|c| matches_group(c, group))
            .map(|c| c.mouse.clone())
            .collect();
        for mouse in mice {
            let train: Vec<_> = cells
                .iter()
                .filter(|c| c.mouse != mouse && matches_group(c, group))
                .collect();
            let test: Vec<_> = cells
                .iter()
                .filter(|c| c.mouse == mouse && matches_group(c, group))
                .collect();
            if test.is_empty() {
                continue;
            }
            let young: Vec<_> = train
                .iter()
                .filter(|c| c.age == 10.0)
                .map(|c| c.beta)
                .collect();
            if young.len() < 3 {
                return Err("not enough independent young training cells".into());
            }
            let mut by_mouse: BTreeMap<String, Vec<&Cell>> = BTreeMap::new();
            for c in &train {
                by_mouse.entry(c.mouse.clone()).or_default().push(c);
            }
            let young_groups: Vec<Vec<f64>> = by_mouse
                .values()
                .filter(|cs| cs[0].age == 10.0)
                .map(|cs| cs.iter().map(|c| c.beta).collect())
                .collect();
            let initial = mean(&young_groups.iter().map(|ys| mean(ys)).collect::<Vec<_>>());
            let iv = young_groups
                .iter()
                .map(|ys| variance(ys) * (ys.len() - 1) as f64)
                .sum::<f64>()
                / young_groups.iter().map(|ys| ys.len() - 1).sum::<usize>() as f64;
            let mut age_means: BTreeMap<u32, Vec<f64>> = BTreeMap::new();
            for cs in by_mouse.values() {
                age_means
                    .entry(cs[0].age as u32)
                    .or_default()
                    .push(mean(&cs.iter().map(|c| c.beta).collect::<Vec<_>>()));
            }
            let residual_sum = by_mouse
                .values()
                .map(|cs| {
                    let ys: Vec<_> = cs.iter().map(|c| c.beta).collect();
                    (mean(&ys) - mean(&age_means[&(cs[0].age as u32)])).powi(2)
                })
                .sum::<f64>();
            let sampling_variance = by_mouse
                .values()
                .map(|cs| {
                    let ys: Vec<_> = cs.iter().map(|c| c.beta).collect();
                    variance(&ys) / ys.len() as f64
                })
                .sum::<f64>()
                / by_mouse.len() as f64;
            let tau = (residual_sum / (by_mouse.len() - age_means.len()).max(1) as f64
                - sampling_variance)
                .max(0.0)
                .sqrt();
            for model in ["one_state", "two_state"] {
                let bounds = if model == "one_state" {
                    vec![(0.005, 0.25), (0.0, 0.03)]
                } else {
                    vec![(0.005, 0.05), (0.05, 0.2), (0.0, 0.02), (0.0, 0.03)]
                };
                let starts = if model == "one_state" {
                    vec![vec![0.02, 0.001], vec![0.05, 0.01]]
                } else {
                    vec![
                        vec![0.02, 0.1, 0.002, 0.001],
                        vec![0.015, 0.15, 0.0065, 0.005],
                        vec![0.04, 0.07, 0.012, 0.01],
                    ]
                };
                // Joint likelihood integrates the shared intercept, then weights animals equally.
                let (_, p) = optimize(
                    |p| {
                        -by_mouse
                            .values()
                            .map(|cs| mouse_logpdf(cs, initial, iv, tau, p) / cs.len() as f64)
                            .sum::<f64>()
                    },
                    &starts,
                    &bounds,
                );
                let ll = mouse_logpdf(&test, initial, iv, tau, &p);
                let dense = mouse_logpdf_dense(&test, initial, iv, tau, &p);
                writeln!(
                    sensitivity,
                    "{group},{mouse},{model},{},{tau},{ll},{dense}",
                    test.len()
                )?;
                writeln!(w,"GSE225171,{group},{mouse},{model},{},{ll},\"{:?}\",continuous_switch_gaussian_random_intercept_approximation",test.len(),p)?;
                if model == "two_state" {
                    let age = test[0].age;
                    let ys: Vec<_> = test.iter().map(|c| c.beta).collect();
                    let grid: Vec<_> = (0..=10000).map(|j| j as f64 / 10000.0).collect();
                    let weights: Vec<_> = grid
                        .iter()
                        .map(|y| polycomb_logpdf(*y, age - 10.0, initial, iv + tau * tau, &p).exp())
                        .collect();
                    let norm = weights.iter().sum::<f64>();
                    let pm = grid.iter().zip(&weights).map(|(y, w)| y * w).sum::<f64>() / norm;
                    let pv = grid
                        .iter()
                        .zip(&weights)
                        .map(|(y, w)| (y - pm).powi(2) * w)
                        .sum::<f64>()
                        / norm;
                    let tail = grid
                        .iter()
                        .zip(&weights)
                        .filter(|(y, _)| **y > 0.089)
                        .map(|(_, w)| *w)
                        .sum::<f64>()
                        / norm;
                    let mut sorted = ys.clone();
                    sorted.sort_by(f64::total_cmp);
                    let (mut cdf, mut wi, mut distance) = (0.0, 0, 0.0);
                    for (j, y) in grid.iter().enumerate() {
                        cdf += weights[j] / norm;
                        while wi < sorted.len() && sorted[wi] <= *y {
                            wi += 1;
                        }
                        distance += (cdf - wi as f64 / sorted.len() as f64).abs() / 10000.0;
                    }
                    writeln!(
                        ppc,
                        "{group},{mouse},{age},{},{pm},{},{},{},{tail},{distance}",
                        mean(&ys),
                        variance(&ys).sqrt(),
                        pv.sqrt(),
                        ys.iter().filter(|y| **y > 0.089).count() as f64 / ys.len() as f64
                    )?;
                }
            }
        }
    }
    Ok(())
}
#[derive(Default, Clone)]
struct Sample {
    accession: String,
    key: String,
    donor: String,
    lineage: String,
    time: f64,
    treatment: String,
    sex: String,
    batch: String,
}
fn fibro_metadata(root: &str) -> Result<BTreeMap<String, Sample>> {
    let mut m = BTreeMap::new();
    let mut current = Sample::default();
    let mut attrs = BTreeMap::new();
    fn finish(c: &mut Sample, a: &BTreeMap<String, String>, m: &mut BTreeMap<String, Sample>) {
        if let (Some(slide), Some(array)) = (a.get("slide"), a.get("array")) {
            c.key = format!("{slide}_{array}");
            c.donor = a.get("cell_line").cloned().unwrap_or_default();
            c.lineage = a.get("cell_line_group").cloned().unwrap_or_default();
            c.time = a
                .get("days_grown_udays")
                .and_then(|s| numeric(s))
                .unwrap_or(f64::NAN);
            c.treatment = a.get("treatments").cloned().unwrap_or_default();
            c.sex = a.get("sex").cloned().unwrap_or_default();
            c.batch = a.get("dneasy_batch").cloned().unwrap_or_default();
            m.insert(c.key.clone(), c.clone());
        }
    }
    for l in reader(&format!("{root}/GSE179847_family.soft.gz"))?.lines() {
        let l = l?;
        if let Some(id) = l.strip_prefix("^SAMPLE = ") {
            finish(&mut current, &attrs, &mut m);
            attrs.clear();
            current = Sample {
                accession: id.into(),
                ..Sample::default()
            };
        } else if let Some(a) = l.strip_prefix("!Sample_characteristics_ch1 = ") {
            if let Some((k, v)) = a.split_once(": ") {
                attrs.insert(k.to_lowercase(), v.into());
            }
        }
    }
    finish(&mut current, &attrs, &mut m);
    Ok(m)
}
fn arrays(root: &str, out: &str, fibro: bool, limit: Option<usize>, matching: &str) -> Result<()> {
    let acc = if fibro { "GSE179847" } else { "GSE73115" };
    let filename = if fibro {
        "GSE179847_matrix.csv.gz"
    } else {
        "GSE73115_processed.txt.gz"
    };
    let sep = if fibro { ',' } else { '\t' };
    let mut r = reader(&format!("{root}/{filename}"))?;
    let mut l = String::new();
    r.read_line(&mut l)?;
    let header = fields(&l, sep);
    let md = if fibro {
        fibro_metadata(root)?
    } else {
        BTreeMap::new()
    };
    let mut columns = vec![];
    let mut samples = vec![];
    for (j, s) in header.iter().enumerate().skip(1) {
        if (fibro && s.ends_with(" beta")) || (!fibro && s.starts_with("Individual_")) {
            columns.push(j);
            let sample = if fibro {
                md.get(s.trim_end_matches(" beta"))
                    .cloned()
                    .ok_or_else(|| format!("matrix sample not in metadata: {s}"))?
            } else {
                let p: Vec<_> = s.split('_').collect();
                Sample {
                    key: s.clone(),
                    donor: p[1].into(),
                    lineage: p[1].into(),
                    time: p[2]
                        .split('-')
                        .next()
                        .ok_or("missing collection year")?
                        .parse::<f64>()?
                        - 1997.0,
                    treatment: "unadjusted_whole_blood".into(),
                    ..Sample::default()
                }
            };
            if !sample.time.is_finite() || sample.donor.is_empty() {
                return Err(format!("invalid time/donor: {s}").into());
            }
            samples.push(sample);
        }
    }
    let mut mw = writer(&format!("{out}/{acc}_samples.csv"))?;
    writeln!(mw, "accession,key,donor,lineage,time,treatment,sex,batch")?;
    for s in &samples {
        writeln!(
            mw,
            "{},{},{},{},{},{},{},{}",
            s.accession, s.key, s.donor, s.lineage, s.time, s.treatment, s.sex, s.batch
        )?;
    }
    let mut by_line: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (j, s) in samples.iter().enumerate() {
        by_line.entry(s.lineage.clone()).or_default().push(j);
    }
    for js in by_line.values_mut() {
        js.sort_by(|a, b| samples[*a].time.total_cmp(&samples[*b].time));
    }
    let pairs: Vec<_> = by_line
        .values()
        .filter(|js| js.len() > 1)
        .map(|js| (js[0], *js.last().unwrap()))
        .filter(|(a, b)| samples[*b].time > samples[*a].time)
        .collect();
    let feature_map = if Path::new("data/derived/human_features.tsv").exists() {
        load_features()?
    } else {
        BTreeMap::new()
    };
    let mut matched: BTreeMap<String, MatchedStratum> = BTreeMap::new();
    let mut sw = writer(&format!("{out}/{acc}_site_drift.csv"))?;
    writeln!(sw,"cpg,valid_pairs,mean_change,sd_change,mean_slope,early_mean,late_mean,early_variance,late_variance")?;
    let mut pw = writer(&format!("{out}/{acc}_replicate_drift.csv"))?;
    let mut sums = vec![(0usize, 0.0, 0.0, 0.0); pairs.len()];
    let (mut sites, mut valid_sites) = (0, 0);
    let mut subset = writer(&format!("{out}/{acc}_kinetic_subset.tsv"))?;
    writeln!(subset, "cpg\tsample\tdonor\tlineage\ttime\tbeta")?;
    let mut forecast = writer(&format!("{out}/{acc}_forecast_pairs.tsv"))?;
    writeln!(
        forecast,
        "cpg\tdonor\tlineage\tstart\tduration\tearly\tlate"
    )?;
    loop {
        l.clear();
        if r.read_line(&mut l)? == 0 {
            break;
        }
        let a = fields(&l, sep);
        if a.len() != header.len() {
            return Err(format!("{acc}: wrong matrix width at row {sites}").into());
        }
        let mut values: Vec<_> = columns
            .iter()
            .map(|j| {
                let b = numeric(&a[*j]);
                let p = a.get(*j + 1).and_then(|s| numeric(s));
                b.filter(|x| (0.0..=1.0).contains(x) && p.is_some_and(|p| p <= 0.01))
            })
            .collect();
        if !fibro {
            for js in by_line.values() {
                for time in [0.0, 10.0] {
                    let group: Vec<_> = js
                        .iter()
                        .copied()
                        .filter(|j| samples[*j].time == time)
                        .collect();
                    let observed: Vec<_> = group.iter().filter_map(|j| values[*j]).collect();
                    let average = if observed.is_empty() {
                        None
                    } else {
                        Some(mean(&observed))
                    };
                    for j in group {
                        values[j] = average;
                    }
                }
            }
        }
        let mut changes = vec![];
        let mut slopes = vec![];
        let mut early = vec![];
        let mut late = vec![];
        for (j, (ia, ib)) in pairs.iter().enumerate() {
            if let (Some(x), Some(y)) = (values[*ia], values[*ib]) {
                let d = y - x;
                changes.push(d);
                slopes.push(d / (samples[*ib].time - samples[*ia].time));
                early.push(x);
                late.push(y);
                let s = &mut sums[j];
                s.0 += 1;
                s.1 += d;
                s.2 += d * d;
                s.3 += d.abs();
            }
        }
        if changes.len() >= 3 {
            if let Some(f) = feature_map.get(&a[0]) {
                if let Some(score) = f.sequence {
                    if !f.island.is_empty()
                        && f.island != "Unknown"
                        && f.chrom.parse::<u32>().is_ok_and(|c| (1..=22).contains(&c))
                    {
                        let regulatory = if f.regulatory.contains("TSS") {
                            "promoter"
                        } else if f.regulatory.contains("Body") {
                            "body"
                        } else {
                            "other"
                        };
                        let contrast_pairs: Vec<_> = pairs
                            .iter()
                            .enumerate()
                            .filter_map(|(j, (ia, ib))| {
                                if fibro && samples[*ia].treatment != "Control" {
                                    return None;
                                }
                                Some((j, values[*ia]?, values[*ib]?))
                            })
                            .collect();
                        if contrast_pairs.len() >= 3 {
                            let contrast_early: Vec<_> =
                                contrast_pairs.iter().map(|(_, x, _)| *x).collect();
                            let (beta_bins, seq_bins, density_bins) = match matching {
                                "coarse" => (10.0, 5.0, 50.0),
                                "fine" => (40.0, 20.0, 200.0),
                                _ => (20.0, 10.0, 100.0),
                            };
                            let density_bin = if matching == "no_density" {
                                0
                            } else {
                                (f.density * density_bins).floor() as u32
                            };
                            let key = format!(
                                "{}:{}:{}:{}:{}:{}",
                                f.chrom,
                                f.island,
                                regulatory,
                                (mean(&contrast_early) * beta_bins).floor() as u32,
                                density_bin,
                                (score * seq_bins).floor() as u32
                            );
                            let g = matched.entry(key).or_insert_with(|| MatchedStratum {
                                values: vec![],
                                q_sum: vec![0.0; pairs.len()],
                                q_n: vec![0; pairs.len()],
                                c_sum: vec![0.0; pairs.len()],
                                c_n: vec![0; pairs.len()],
                            });
                            let absolute: Vec<_> = contrast_pairs
                                .iter()
                                .map(|(_, x, y)| (y - x).abs())
                                .collect();
                            g.values.push((f.esl, mean(&absolute)));
                            for (j, x, y) in contrast_pairs {
                                {
                                    if f.esl {
                                        g.q_sum[j] += (y - x).abs();
                                        g.q_n[j] += 1;
                                    } else {
                                        g.c_sum[j] += (y - x).abs();
                                        g.c_n[j] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if changes.len() >= 3 {
            valid_sites += 1;
            writeln!(
                sw,
                "{},{},{},{},{},{},{},{},{}",
                a[0],
                changes.len(),
                mean(&changes),
                variance(&changes).sqrt(),
                mean(&slopes),
                mean(&early),
                mean(&late),
                variance(&early),
                variance(&late)
            )?;
        }
        // Denser, outcome-independent paired subset for held-out stability prediction.
        if sites % 256 == 0 {
            for (ia, ib) in &pairs {
                if fibro && samples[*ia].treatment != "Control" {
                    continue;
                }
                if let (Some(x), Some(y)) = (values[*ia], values[*ib]) {
                    let s = &samples[*ia];
                    writeln!(
                        forecast,
                        "{}\t{}\t{}\t{}\t{}\t{x}\t{y}",
                        a[0],
                        s.donor,
                        s.lineage,
                        s.time,
                        samples[*ib].time - s.time
                    )?;
                }
            }
        }
        // Fixed, outcome-independent systematic subset for compact inference: every 4096th row.
        if sites % 4096 == 0 {
            for (j, v) in values.iter().enumerate() {
                if let Some(v) = v {
                    let s = &samples[j];
                    writeln!(
                        subset,
                        "{}\t{}\t{}\t{}\t{}\t{v}",
                        a[0], s.key, s.donor, s.lineage, s.time
                    )?;
                }
            }
        }
        sites += 1;
        if limit.is_some_and(|n| sites >= n) {
            break;
        }
    }
    writeln!(pw,"donor,lineage,treatment,time_start,time_end,cpg_pairs,mean_change,rms_change,mean_absolute_change")?;
    for ((a, b), (n, d, d2, abs)) in pairs.iter().zip(&sums) {
        let s = &samples[*a];
        writeln!(
            pw,
            "{},{},{},{},{},{n},{},{},{}",
            s.donor,
            s.lineage,
            s.treatment,
            s.time,
            samples[*b].time,
            d / (*n as f64),
            (d2 / (*n as f64)).sqrt(),
            abs / (*n as f64)
        )?;
    }
    if !matched.is_empty() {
        matched_results(acc, out, matched, &pairs, &samples)?;
    }
    let donors: BTreeSet<_> = samples.iter().map(|s| &s.donor).collect();
    println!("{acc}: {sites} sites, {valid_sites} with >=3 valid pairs, {} samples, {} donors, {} longitudinal groups",samples.len(),donors.len(),pairs.len());
    Ok(())
}
fn simulation_checks(out: &str) -> Result<()> {
    let site = Site {
        target: false,
        initial: 0.0,
        reference: 0.0,
        q: 1.0,
        context: vec![],
        away: 0.2,
        restoration: 0.1,
        old_multiplier: 4.0,
        loadings: vec![],
    };
    let mut w = writer(&format!("{out}/simulation_checks.csv"))?;
    writeln!(w,"mechanism,cells,survived,methylation_among_survivors,first_passage_fraction,analytic_survivor_probability,protection_fit_probability,selection_fit_probability,snapshot_identifiable")?;
    let mut text=String::from("<!doctype html><meta charset=utf-8><title>Identifiability checks</title><h1>Protection versus selection</h1><p>Seed 20261003; 20,000 cells per mechanism; one externally weighted binary site, horizon 5. Rates away=0.2, return=0.1. These are numerical counterexamples, not a comprehensive identifiability proof. Observation of survivors at one time does not determine the generating mechanism.</p><table><tr><th>Generator</th><th>Survived</th><th>Survivor beta</th><th>Protection fit</th><th>Selection fit</th></tr>");
    for (j, (name, m)) in [
        (
            "protection",
            Mechanism {
                protection: 1.2,
                recovery: 0.7,
                ..Mechanism::default()
            },
        ),
        (
            "selection",
            Mechanism {
                selection_baseline: 0.2,
                selection_strength: 2.0,
                ..Mechanism::default()
            },
        ),
        (
            "both",
            Mechanism {
                protection: 1.2,
                recovery: 0.7,
                selection_baseline: 0.2,
                selection_strength: 2.0,
                ..Mechanism::default()
            },
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let paths: Vec<_> = (0..20000)
            .map(|i| {
                simulate(
                    std::slice::from_ref(&site),
                    &m,
                    5.0,
                    0.5,
                    20261003 + j as u64 * 100000 + i,
                )
                .map_err(std::io::Error::other)
            })
            .collect::<std::result::Result<_, _>>()?;
        let alive: Vec<_> = paths.iter().filter(|p| p.death.is_none()).collect();
        let empirical = alive.iter().filter(|p| p.states[0]).count() as f64 / alive.len() as f64;
        let (a, b) = (
            site.away * (-m.protection).exp(),
            site.restoration * m.recovery.exp(),
        );
        let exact = survivor_probability(
            0.0,
            a,
            b,
            m.selection_baseline,
            m.selection_baseline * m.selection_strength.exp(),
            5.0,
        )
        .map_err(std::io::Error::other)?
        .0;
        // Both candidates fit survivor methylation only. Killing counts are deliberately absent.
        let (_, prot) = optimize(
            |p| (probability(0.0, 0.2 * (-p[0]).exp(), 0.1, 5.0) - empirical).powi(2),
            &[vec![1.0]],
            &[(0.0, 8.0)],
        );
        let prot_p = probability(0.0, 0.2 * (-prot[0]).exp(), 0.1, 5.0);
        let (_, sel) = optimize(
            |p| {
                let prob = survivor_probability(0.0, 0.2, 0.1, 0.2, 0.2 * p[0].exp(), 5.0)
                    .unwrap()
                    .0;
                (prob - empirical).powi(2)
            },
            &[vec![2.0]],
            &[(0.0, 8.0)],
        );
        let sel_p = survivor_probability(0.0, 0.2, 0.1, 0.2, 0.2 * sel[0].exp(), 5.0)
            .unwrap()
            .0;
        writeln!(
            w,
            "{name},20000,{},{empirical},{},{exact},{prot_p},{sel_p},false",
            alive.len(),
            paths.iter().filter(|p| p.first_passage.is_some()).count() as f64 / paths.len() as f64
        )?;
        text.push_str(&format!("<tr><td>{name}</td><td>{}</td><td>{empirical:.6}</td><td>{prot_p:.6}</td><td>{sel_p:.6}</td></tr>",alive.len()));
    }
    text.push_str("</table><p>Binary observations, binomial/beta-binomial read counts with fixed dispersion, and noisy aggregate arrays all have the same likelihood when endpoint probability is matched. This equivalence does not assert that realistic multisite longitudinal data are always unidentifiable. Observed disappearance or clone representation and generation-time methylation can help, but those require an expanded observation model. Simulated death differs strongly despite matched survivor methylation.</p>");
    fs::write(format!("{out}/identifiability_report.html"), text)?;
    let mut w = writer(&format!("{out}/continuous_sensitivity.csv"))?;
    writeln!(
        w,
        "dt,cells,endpoint_methylation,grid_first_passage_fraction"
    )?;
    let m = Mechanism {
        maintenance: Maintenance::Continuous {
            reversion: 0.2,
            mean: 0.0,
            noise: 0.3,
        },
        ..Mechanism::default()
    };
    for dt in [0.1, 0.025, 0.00625] {
        let ps: Vec<_> = (0..5000)
            .map(|i| {
                simulate_grid(std::slice::from_ref(&site), &m, 5.0, 0.5, dt, 20261003 + i).unwrap()
            })
            .collect();
        writeln!(
            w,
            "{dt},5000,{},{}",
            ps.iter().filter(|p| p.states[0]).count() as f64 / ps.len() as f64,
            ps.iter().filter(|p| p.first_passage.is_some()).count() as f64 / ps.len() as f64
        )?;
    }
    Ok(())
}

// Compact real-data M1 comparison, prespecified every 16th site in the systematic
// kinetic subset (every 65,536th input probe). External q inference is gated separately.
fn fit_real_baseline(out: &str) -> Result<()> {
    let mut w = writer(&format!("{out}/kinetic_comparison.csv"))?;
    writeln!(
        w,
        "dataset,cpg,held_out_donor,n,model,mean_log_predictive_density,rmse,status"
    )?;
    for acc in ["GSE179847", "GSE73115"] {
        let mut sample_map = BTreeMap::new();
        let mut lines = reader(&format!("{out}/{acc}_samples.csv"))?.lines();
        lines.next();
        for l in lines {
            let a = fields(&l?, ',');
            sample_map.insert(a[1].clone(), a[5].clone());
        }
        let mut data: BTreeMap<String, Vec<(String, String, f64, f64)>> = BTreeMap::new();
        let mut lines = reader(&format!("{out}/{acc}_kinetic_subset.tsv"))?.lines();
        lines.next();
        let mut seen = BTreeSet::new();
        let mut site_index = 0;
        let mut take = false;
        let mut previous = String::new();
        for l in lines {
            let a = fields(&l?, '\t');
            if a[0] != previous {
                take = site_index % 16 == 0;
                site_index += 1;
                previous = a[0].clone();
            }
            if !take {
                continue;
            }
            if acc == "GSE179847" && sample_map[&a[1]] != "Control" {
                continue;
            }
            // Blood technical replicate rows already share pooled beta. Keep one per person/time.
            let key = format!("{}_{}_{}", a[0], a[2], a[4]);
            if acc == "GSE73115" && !seen.insert(key) {
                continue;
            }
            data.entry(a[0].clone()).or_default().push((
                a[2].clone(),
                a[3].clone(),
                a[4].parse()?,
                a[5].parse()?,
            ));
        }
        for (cpg, rows) in data {
            let donors: BTreeSet<_> = rows.iter().map(|r| r.0.clone()).collect();
            for donor in donors {
                let train: Vec<_> = rows.iter().filter(|r| r.0 != donor).collect();
                let test: Vec<_> = rows.iter().filter(|r| r.0 == donor).collect();
                if train.len() < 5 || test.is_empty() {
                    continue;
                }
                let mut first: BTreeMap<&str, (f64, f64)> = BTreeMap::new();
                for r in &train {
                    let e = first.entry(&r.1).or_insert((r.2, r.3));
                    if r.2 < e.0 {
                        *e = (r.2, r.3);
                    }
                }
                let young: Vec<_> = first.values().map(|v| v.1).collect();
                let init = mean(&young);
                let error = variance(&young).sqrt().max(0.02);
                let sites = vec![Site {
                    target: init >= 0.5,
                    initial: init,
                    reference: init,
                    q: 0.0,
                    context: vec![],
                    away: 0.01,
                    restoration: 0.01,
                    old_multiplier: 1.0,
                    loadings: vec![],
                }];
                let kinetic: Vec<_> = train
                    .iter()
                    .map(|r| KineticRow {
                        site: 0,
                        time: r.2,
                        replicate: r.0.clone(),
                        observation: Observation::Array {
                            beta: r.3,
                            sd: error,
                        },
                    })
                    .collect();
                let fit =
                    fit_kinetics(&sites, &kinetic, false, 0.3).map_err(std::io::Error::other)?;
                for model in ["constant_young", "M1"] {
                    let mut ll = 0.0;
                    let mut se = 0.0;
                    for r in &test {
                        let row = KineticRow {
                            site: 0,
                            time: r.2,
                            replicate: r.0.clone(),
                            observation: Observation::Array {
                                beta: r.3,
                                sd: error,
                            },
                        };
                        let p = if model == "M1" {
                            kinetic_prediction(&fit, &sites, &row).map_err(std::io::Error::other)?
                        } else {
                            init
                        };
                        ll += log_likelihood(&row.observation, p).map_err(std::io::Error::other)?;
                        se += (r.3 - p).powi(2);
                    }
                    writeln!(
                        w,
                        "{acc},{cpg},{donor},{},{model},{},{},exploratory_unadjusted_arrays",
                        test.len(),
                        ll / test.len() as f64,
                        (se / test.len() as f64).sqrt()
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn reports(out: &str) -> Result<()> {
    use rand::{Rng, SeedableRng};
    fn interval(values: &[f64]) -> (f64, f64) {
        let mut rng = rand::rngs::StdRng::seed_from_u64(20261003);
        let mut boots: Vec<f64> = (0..4000)
            .map(|_| {
                (0..values.len())
                    .map(|_| values[rng.gen_range(0..values.len())])
                    .sum::<f64>()
                    / values.len() as f64
            })
            .collect();
        boots.sort_by(f64::total_cmp);
        (boots[100], boots[3900])
    }
    let mut w = writer(&format!("{out}/replicate_bootstrap.csv"))?;
    writeln!(
        w,
        "dataset,group,metric,biological_replicates,estimate,bootstrap_lower_95,bootstrap_upper_95"
    )?;
    let mut html=String::from("<!doctype html><meta charset=utf-8><title>Empirical drift report</title><h1>Empirical drift and predictive checks</h1><p>Exploratory analysis, seed 20261003; bootstrap resamples biological replicates, not CpGs or cells. Bulk results are unadjusted for batch, sex and cell composition. Read the limitations in results.md.</p><table><tr><th>Dataset / group</th><th>Metric</th><th>Estimate</th><th>95% replicate bootstrap</th></tr>");
    let mut emit = |dataset: &str, group: &str, metric: &str, values: &[f64]| -> Result<()> {
        let (lo, hi) = interval(values);
        let estimate = mean(values);
        writeln!(
            w,
            "{dataset},{group},{metric},{},{estimate},{lo},{hi}",
            values.len()
        )?;
        html.push_str(&format!("<tr><td>{dataset} / {group}</td><td>{metric}</td><td>{estimate:.6}</td><td>[{lo:.6}, {hi:.6}]</td></tr>"));
        Ok(())
    };
    for acc in ["GSE179847", "GSE73115"] {
        let mut lines = reader(&format!("{out}/{acc}_replicate_drift.csv"))?.lines();
        lines.next();
        let mut by: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for l in lines {
            let a = fields(&l?, ',');
            if acc == "GSE179847" && a[2] != "Control" {
                continue;
            }
            by.entry(a[0].clone()).or_default().push(a[8].parse()?);
        }
        let values: Vec<_> = by.values().map(|v| mean(v)).collect();
        emit(
            acc,
            "controls_or_whole_blood",
            "mean_absolute_endpoint_change",
            &values,
        )?;
    }
    for group in ["slow", "fast"] {
        let mut lines = reader(&format!("{out}/model_comparison.csv"))?.lines();
        lines.next();
        let mut pairs: BTreeMap<String, BTreeMap<String, (f64, usize)>> = BTreeMap::new();
        for l in lines {
            let a = fields(&l?, ',');
            if a[1] != group {
                continue;
            }
            pairs
                .entry(a[2].clone())
                .or_default()
                .insert(a[3].clone(), (a[5].parse()?, a[4].parse()?));
        }
        let values: Vec<_> = pairs
            .values()
            .map(|p| (p["two_state"].0 - p["one_state"].0) / p["one_state"].1 as f64)
            .collect();
        emit(
            "GSE225171",
            group,
            "two_minus_one_state_lpd_per_cell",
            &values,
        )?;
    }
    for acc in ["GSE179847", "GSE73115"] {
        let mut lines = reader(&format!("{out}/kinetic_comparison.csv"))?.lines();
        lines.next();
        let mut pairs: BTreeMap<(String, String), BTreeMap<String, f64>> = BTreeMap::new();
        for l in lines {
            let a = fields(&l?, ',');
            if a[0] != acc {
                continue;
            }
            pairs
                .entry((a[1].clone(), a[2].clone()))
                .or_default()
                .insert(a[4].clone(), a[5].parse()?);
        }
        let mut by: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for ((_, donor), p) in pairs {
            by.entry(donor)
                .or_default()
                .push(p["M1"] - p["constant_young"]);
        }
        let values: Vec<_> = by.values().map(|v| mean(v)).collect();
        emit(
            acc,
            "systematic_probe_subset",
            "M1_minus_constant_lpd",
            &values,
        )?;
    }
    for acc in ["GSE179847", "GSE73115"] {
        let path = format!("{out}/{acc}_matched_replicates.csv");
        if Path::new(&path).exists() {
            let mut lines = reader(&path)?.lines();
            lines.next();
            let mut by: BTreeMap<String, Vec<f64>> = BTreeMap::new();
            for line in lines {
                let a = fields(&line?, ',');
                by.entry(a[0].clone()).or_default().push(a[2].parse()?);
            }
            let values: Vec<_> = by.values().map(|v| mean(v)).collect();
            emit(
                acc,
                "q_esl_matched_context",
                "q_minus_control_absolute_change",
                &values,
            )?;
        }
    }
    html.push_str("</table><h2>Scientific gates</h2><p>F2/F8: external ESL stability masks compared with within-stratum random masks; matching uses chromosome, CpG island class, regulatory class, early beta, probe-sequence density and symmetrized sequence rank. F1 held-out mechanistic effects of essentiality or identity are not established. ESL stability is not functional ground truth. F3: one-snapshot protection/selection counterexample reproduced; multisite recovery not established. F4: approximate slow-fast distribution baseline evaluated with mouse-held-out prediction; final slow-group estimate is positive but fast-group gain is absent; see numerical approximation and bootstrap limitations. F5: low-rank numerical support tested, no covariate-adjusted biological conclusion. F6/F7: no independent identity or clock outcome comparison; first-passage simulation does not establish identity loss.</p><p><a href=identifiability_report.html>Identifiability report</a></p>");
    fs::write(format!("{out}/empirical_drift_report.html"), html)?;
    Ok(())
}

#[derive(Clone)]
struct Feature {
    chrom: String,
    position: String,
    density: f64,
    sequence: Option<f64>,
    island: String,
    regulatory: String,
    esl: bool,
}
fn annotations(root: &str, out: &str) -> Result<BTreeMap<String, Feature>> {
    let raw = std::process::Command::new("unzip")
        .args([
            "-p",
            &format!("{root}/esl_annotation.xlsx"),
            "xl/sharedStrings.xml",
        ])
        .output()?;
    if !raw.status.success() {
        return Err("cannot read ESL workbook shared strings".into());
    }
    let xml = String::from_utf8(raw.stdout)?;
    // CpG IDs contain only ASCII letters/digits and have no XML escaping; exclude all other strings.
    let esl: BTreeSet<String> = xml
        .split('>')
        .filter_map(|part| part.split('<').next())
        .filter(|s| {
            s.len() == 10 && s.starts_with("cg") && s[2..].bytes().all(|c| c.is_ascii_digit())
        })
        .map(String::from)
        .collect();
    if esl.len() != 31744 {
        return Err(format!("unexpected ESL annotation count: {}", esl.len()).into());
    }
    let pdf = format!("{root}/sequence_fidelity_supplement.pdf");
    let converted = std::process::Command::new("pdftotext")
        .args(["-layout", &pdf, "-"])
        .output()?;
    if !converted.status.success() {
        return Err("sequence supplement is not a readable PDF".into());
    }
    let mut ranking = BTreeMap::new();
    for l in String::from_utf8(converted.stdout)?.lines() {
        let a: Vec<_> = l.split_whitespace().collect();
        if a.len() == 4
            && a[1].len() == 6
            && a[1].bytes().all(|c| b"ACGT".contains(&c))
            && a[2].len() == 2
            && a[3].len() == 2
        {
            if let Ok(rank) = a[0].parse::<usize>() {
                if (1..=225).contains(&rank) {
                    ranking.insert(a[1].to_string(), 1.0 - (rank - 1) as f64 / 224.0);
                }
            }
        }
    }
    if ranking.len() != 225 {
        return Err(format!("unexpected sequence ranking count: {}", ranking.len()).into());
    }
    let mut features: BTreeMap<String, Feature> = BTreeMap::new();
    // Read the 450K annotation first; EPIC adds sequences/coordinates but does not overwrite genomic context.
    for acc in ["GSE73115", "GSE179847"] {
        let mut lines = reader(&format!("{root}/{acc}_family.soft.gz"))?.lines();
        let mut header = None;
        for line in lines.by_ref() {
            if line? == "!platform_table_begin" {
                header = Some(fields(
                    &lines.next().ok_or("missing platform header")??,
                    '\t',
                ));
                break;
            }
        }
        let header = header.ok_or("missing platform annotation")?;
        let indices: BTreeMap<&str, usize> = header
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        for line in lines {
            let line = line?;
            if line == "!platform_table_end" {
                break;
            }
            let a = fields(&line, '\t');
            if a.len() != header.len() {
                return Err("invalid platform annotation width".into());
            }
            let id = &a[indices["ID"]];
            if !id.starts_with("cg") {
                continue;
            }
            let seq = &a[indices["Forward_Sequence"]];
            let motif = flanking_motif(seq);
            let score = motif.and_then(|m| {
                let reverse: String = m
                    .chars()
                    .rev()
                    .map(|c| match c {
                        'A' => 'T',
                        'T' => 'A',
                        'C' => 'G',
                        'G' => 'C',
                        _ => unreachable!(),
                    })
                    .collect();
                Some((ranking.get(&m)? + ranking.get(&reverse)?) / 2.0)
            });
            let chrom = a[indices["CHR"]].clone();
            let position = a[indices["MAPINFO"]].clone();
            let island = indices
                .get("Relation_to_UCSC_CpG_Island")
                .map(|j| {
                    if a[*j].is_empty() {
                        "OpenSea".into()
                    } else {
                        a[*j].clone()
                    }
                })
                .unwrap_or_else(|| "Unknown".into());
            let regulatory = indices
                .get("UCSC_RefGene_Group")
                .map(|j| a[*j].clone())
                .unwrap_or_default();
            let clean = seq.replace("[CG]", "CG");
            let density = clean.as_bytes().windows(2).filter(|w| *w == b"CG").count() as f64
                / clean.len().max(1) as f64;
            if let Some(f) = features.get_mut(id) {
                if f.sequence.is_none() {
                    f.sequence = score;
                }
            } else {
                features.insert(
                    id.clone(),
                    Feature {
                        chrom,
                        position,
                        density,
                        sequence: score,
                        island,
                        regulatory,
                        esl: esl.contains(id),
                    },
                );
            }
        }
    }
    fs::create_dir_all("data/derived")?;
    let mut w = writer("data/derived/human_features.tsv")?;
    writeln!(w,"cpg\tchrom\tposition_hg19\tcpg_density_probe_sequence\tsequence_fidelity_rank_score\tcgi_context\tregulatory_class\tq_esl")?;
    for (id, f) in &features {
        writeln!(
            w,
            "{id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f.chrom,
            f.position,
            f.density,
            f.sequence
                .map(|v| v.to_string())
                .unwrap_or_else(|| "NA".into()),
            f.island,
            f.regulatory,
            u8::from(f.esl)
        )?;
    }
    let mut w = writer(&format!("{out}/annotation_coverage.csv"))?;
    writeln!(w, "features,esls,sequence_rank_covered,complete_context")?;
    writeln!(
        w,
        "{},{},{},{}",
        features.len(),
        features.values().filter(|f| f.esl).count(),
        features.values().filter(|f| f.sequence.is_some()).count(),
        features
            .values()
            .filter(|f| f.sequence.is_some() && !f.island.is_empty() && f.island != "Unknown")
            .count()
    )?;
    Ok(features)
}
fn flanking_motif(seq: &str) -> Option<String> {
    let (left, right) = seq.split_once("[CG]")?;
    if left.len() < 2 || right.len() < 2 {
        return None;
    }
    let m = format!("{}CG{}", &left[left.len() - 2..], &right[..2]);
    m.bytes().all(|c| b"ACGT".contains(&c)).then_some(m)
}

#[derive(Default)]
struct MatchedStratum {
    values: Vec<(bool, f64)>,
    q_sum: Vec<f64>,
    q_n: Vec<usize>,
    c_sum: Vec<f64>,
    c_n: Vec<usize>,
}
fn load_features() -> Result<BTreeMap<String, Feature>> {
    let mut m = BTreeMap::new();
    let mut lines = reader("data/derived/human_features.tsv")?.lines();
    lines.next();
    for line in lines {
        let a = fields(&line?, '\t');
        m.insert(
            a[0].clone(),
            Feature {
                chrom: a[1].clone(),
                position: a[2].clone(),
                density: a[3].parse()?,
                sequence: numeric(&a[4]),
                island: a[5].clone(),
                regulatory: a[6].clone(),
                esl: a[7] == "1",
            },
        );
    }
    Ok(m)
}
fn matched_results(
    acc: &str,
    out: &str,
    groups: BTreeMap<String, MatchedStratum>,
    pairs: &[(usize, usize)],
    samples: &[Sample],
) -> Result<()> {
    use rand::{seq::SliceRandom, SeedableRng};
    let supported: Vec<_> = groups
        .values()
        .filter(|g| g.values.iter().any(|v| v.0) && g.values.iter().any(|v| !v.0))
        .collect();
    let supported_q = supported
        .iter()
        .map(|g| g.values.iter().filter(|v| v.0).count())
        .sum::<usize>();
    let total_q = groups
        .values()
        .map(|g| g.values.iter().filter(|v| v.0).count())
        .sum::<usize>();
    let contrast =
        |groups: &[&MatchedStratum], shuffle: bool, rng: &mut rand::rngs::StdRng| -> f64 {
            let mut sum = 0.0;
            let mut weight = 0.0;
            for g in groups {
                let qn = g.values.iter().filter(|v| v.0).count();
                let mut y: Vec<_> = g.values.iter().map(|v| v.1).collect();
                let (q, c) = if shuffle {
                    y.shuffle(rng);
                    (mean(&y[..qn]), mean(&y[qn..]))
                } else {
                    (
                        mean(
                            &g.values
                                .iter()
                                .filter(|v| v.0)
                                .map(|v| v.1)
                                .collect::<Vec<_>>(),
                        ),
                        mean(
                            &g.values
                                .iter()
                                .filter(|v| !v.0)
                                .map(|v| v.1)
                                .collect::<Vec<_>>(),
                        ),
                    )
                };
                sum += qn as f64 * (q - c);
                weight += qn as f64;
            }
            sum / weight
        };
    if supported_q == 0 {
        return Err(format!("{acc}: no matched support for the ESL mask").into());
    }
    let mut rng = rand::rngs::StdRng::seed_from_u64(20261003);
    let actual = contrast(&supported, false, &mut rng);
    let mut nulls: Vec<_> = (0..500)
        .map(|_| contrast(&supported, true, &mut rng))
        .collect();
    let permutation_p =
        (1 + nulls.iter().filter(|x| x.abs() >= actual.abs()).count()) as f64 / 501.0;
    nulls.sort_by(f64::total_cmp);
    let mut w = writer(&format!("{out}/{acc}_matched_mask.csv"))?;
    writeln!(w,"mask,eligible_protected_sites,supported_protected_sites,supported_strata,contrast_absolute_change_q_minus_control,matched_null_lower_95,matched_null_upper_95,two_sided_permutation_p,status")?;
    writeln!(w,"q_esl,{total_q},{supported_q},{},{actual},{},{},{permutation_p},exploratory_stability_annotation_not_functional_essentiality",supported.len(),nulls[12],nulls[487])?;
    let mut w = writer(&format!("{out}/{acc}_matched_replicates.csv"))?;
    writeln!(w, "donor,lineage,contrast_absolute_change,matched_q_weight")?;
    for (j, (a, _)) in pairs.iter().enumerate() {
        let mut sum = 0.0;
        let mut weight = 0;
        for g in &supported {
            if g.q_n[j] > 0 && g.c_n[j] > 0 {
                sum +=
                    g.q_n[j] as f64 * (g.q_sum[j] / g.q_n[j] as f64 - g.c_sum[j] / g.c_n[j] as f64);
                weight += g.q_n[j];
            }
        }
        if weight > 0 {
            writeln!(
                w,
                "{},{},{},{}",
                samples[*a].donor,
                samples[*a].lineage,
                sum / weight as f64,
                weight
            )?;
        }
    }
    Ok(())
}

// Empirical F1 gate: forecast log absolute paired change, not infer microscopic rates.
#[derive(Clone)]
struct ForecastPair {
    cpg: String,
    donor: String,
    early: f64,
    change: f64,
    x: Vec<f64>,
    q: bool,
    chrom: u32,
    island: String,
    regulatory: String,
    sequence: f64,
    density: f64,
}
#[derive(Clone)]
struct NormalEquations {
    gram: Vec<Vec<f64>>,
    rhs: Vec<f64>,
    yy: f64,
    weight: f64,
}
impl NormalEquations {
    fn new(p: usize) -> Self {
        Self {
            gram: vec![vec![0.0; p]; p],
            rhs: vec![0.0; p],
            yy: 0.0,
            weight: 0.0,
        }
    }
    fn add(&mut self, row: &ForecastPair, weight: f64) {
        self.add_xy(&row.x, (row.change.abs() + 0.001).ln(), weight);
    }
    fn add_xy(&mut self, x: &[f64], y: f64, weight: f64) {
        self.yy += weight * y * y;
        self.weight += weight;
        for (i, xi) in x.iter().enumerate().filter(|(_, x)| **x != 0.0) {
            self.rhs[i] += weight * xi * y;
            for (j, xj) in x.iter().enumerate().take(i + 1).filter(|(_, x)| **x != 0.0) {
                self.gram[i][j] += weight * xi * xj;
                if i != j {
                    self.gram[j][i] += weight * xi * xj;
                }
            }
        }
    }
    fn subtract(&self, other: &Self) -> Self {
        Self {
            gram: self
                .gram
                .iter()
                .zip(&other.gram)
                .map(|(a, b)| a.iter().zip(b).map(|(x, y)| x - y).collect())
                .collect(),
            rhs: self
                .rhs
                .iter()
                .zip(&other.rhs)
                .map(|(x, y)| x - y)
                .collect(),
            yy: self.yy - other.yy,
            weight: self.weight - other.weight,
        }
    }
    fn fit(&self, p: usize, penalty: f64) -> Result<(Vec<f64>, f64)> {
        let mut a: Vec<Vec<f64>> = (0..p)
            .map(|i| {
                let mut row = self.gram[i][..p].to_vec();
                row[i] += if i == 0 { 1e-10 } else { penalty * self.weight };
                row.push(self.rhs[i]);
                row
            })
            .collect();
        for k in 0..p {
            let pivot = (k..p)
                .max_by(|i, j| a[*i][k].abs().total_cmp(&a[*j][k].abs()))
                .unwrap();
            a.swap(k, pivot);
            if a[k][k].abs() < 1e-12 {
                return Err("singular forecast normal equations".into());
            }
            let d = a[k][k];
            for j in k..=p {
                a[k][j] /= d;
            }
            for i in 0..p {
                if i != k {
                    let d = a[i][k];
                    for j in k..=p {
                        a[i][j] -= d * a[k][j];
                    }
                }
            }
        }
        let b: Vec<_> = a.iter().map(|r| r[p]).collect();
        let linear: f64 = b.iter().enumerate().map(|(i, v)| v * self.rhs[i]).sum();
        let quadratic: f64 = b
            .iter()
            .enumerate()
            .map(|(i, v)| {
                v * b
                    .iter()
                    .enumerate()
                    .map(|(j, u)| self.gram[i][j] * u)
                    .sum::<f64>()
            })
            .sum();
        let variance = ((self.yy - 2.0 * linear + quadratic) / self.weight).max(1e-6);
        Ok((b, variance))
    }
}
fn forecast_rows(
    acc: &str,
    out: &str,
    features: &BTreeMap<String, Feature>,
) -> Result<Vec<ForecastPair>> {
    let mut rows = vec![];
    let mut lines = reader(&format!("{out}/{acc}_forecast_pairs.tsv"))?.lines();
    lines.next();
    for line in lines {
        let a = fields(&line?, '\t');
        let Some(f) = features.get(&a[0]) else {
            continue;
        };
        let Some(sequence) = f.sequence else {
            continue;
        };
        let Ok(chrom) = f.chrom.parse::<u32>() else {
            continue;
        };
        if !(1..=22).contains(&chrom) || f.island.is_empty() || f.island == "Unknown" {
            continue;
        }
        let early: f64 = a[5].parse()?;
        let late: f64 = a[6].parse()?;
        let regulatory = if f.regulatory.contains("TSS") {
            "promoter"
        } else if f.regulatory.contains("Body") {
            "body"
        } else {
            "other"
        };
        // Only information available at the earlier observation and fixed external annotations.
        let mut x = vec![1.0, early, early * early];
        x.extend((1..10).map(|i| (early - i as f64 / 10.0).max(0.0)));
        x.extend([
            sequence,
            sequence * sequence,
            f.density * 10.0,
            a[3].parse::<f64>()? / 100.0,
            a[4].parse::<f64>()? / 100.0,
            f64::from(a[1].starts_with("SURF1")),
            f64::from(a[2].ends_with("ox3")),
        ]);
        x.extend(
            ["Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf"]
                .map(|s| f64::from(f.island == s)),
        );
        x.extend([
            f64::from(regulatory == "promoter"),
            f64::from(regulatory == "body"),
        ]);
        x.extend((2..=22).map(|c| f64::from(chrom == c)));
        x.push(f64::from(f.esl));
        rows.push(ForecastPair {
            cpg: a[0].clone(),
            donor: a[1].clone(),
            early,
            change: late - early,
            x,
            q: f.esl,
            chrom,
            island: f.island.clone(),
            regulatory: regulatory.into(),
            sequence,
            density: f.density,
        });
    }
    if rows.is_empty() {
        return Err(format!("{acc}: no forecast rows").into());
    }
    Ok(rows)
}
fn forecasting(out: &str) -> Result<()> {
    let features = load_features()?;
    let mut w = writer(&format!("{out}/forecast_comparison.csv"))?;
    writeln!(w,"dataset,split,penalty,held_out_donor,n,n_esl,log_change_mse_gain,esl_log_change_mse_gain,log_predictive_gain,esl_coefficient")?;
    let mut mw = writer(&format!("{out}/matching_sensitivity.csv"))?;
    writeln!(
        mw,
        "dataset,matching,held_out_donor,matched_esl_sites,absolute_change_contrast"
    )?;
    for acc in ["GSE179847", "GSE73115"] {
        let rows = forecast_rows(acc, out, &features)?;
        let donors: Vec<_> = rows
            .iter()
            .map(|r| r.donor.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let counts: BTreeMap<_, _> = donors
            .iter()
            .map(|d| (d.clone(), rows.iter().filter(|r| &r.donor == d).count()))
            .collect();
        let p = rows[0].x.len();
        let sites = rows.iter().map(|r| &r.cpg).collect::<BTreeSet<_>>().len();
        let esl = rows
            .iter()
            .filter(|r| r.q)
            .map(|r| &r.cpg)
            .collect::<BTreeSet<_>>()
            .len();
        println!(
            "{acc}: forecast {sites} annotated sites ({esl} ESL), {} pairs, {} donors",
            rows.len(),
            donors.len()
        );
        let mut all = NormalEquations::new(p);
        let mut by_donor: BTreeMap<_, _> = donors
            .iter()
            .map(|d| (d.clone(), NormalEquations::new(p)))
            .collect();
        for r in &rows {
            let weight = 1.0 / counts[&r.donor] as f64;
            all.add(r, weight);
            by_donor.get_mut(&r.donor).unwrap().add(r, weight);
        }
        for split in ["donor", "donor_and_chromosome"] {
            let folds = if split == "donor" { donors.len() } else { 5 };
            for fold in 0..folds {
                for genome_fold in 0..if split == "donor" { 1 } else { 5 } {
                    let heldout = |r: &ForecastPair| {
                        if split == "donor" {
                            r.donor == donors[fold]
                        } else {
                            donors.binary_search(&r.donor).unwrap() % 5 == fold
                        }
                    };
                    let test: Vec<_> = rows
                        .iter()
                        .filter(|r| {
                            heldout(r) && (split == "donor" || r.chrom as usize % 5 == genome_fold)
                        })
                        .collect();
                    let train = if split == "donor" {
                        all.subtract(&by_donor[&donors[fold]])
                    } else {
                        let mut train = NormalEquations::new(p);
                        for r in rows
                            .iter()
                            .filter(|r| !heldout(r) && r.chrom as usize % 5 != genome_fold)
                        {
                            train.add(r, 1.0 / counts[&r.donor] as f64);
                        }
                        train
                    };
                    for penalty in [0.001, 0.01, 0.1] {
                        let (b0, v0) = train.fit(p - 1, penalty)?;
                        let (b1, v1) = train.fit(p, penalty)?;
                        let mut scores: BTreeMap<String, (usize, usize, f64, f64, f64)> =
                            BTreeMap::new();
                        for r in &test {
                            let y = (r.change.abs() + 0.001).ln();
                            let predict =
                                |b: &[f64]| b.iter().zip(&r.x).map(|(a, b)| a * b).sum::<f64>();
                            let e0 = y - predict(&b0);
                            let e1 = y - predict(&b1);
                            let s = scores.entry(r.donor.clone()).or_default();
                            s.0 += 1;
                            s.1 += usize::from(r.q);
                            s.2 += e0 * e0 - e1 * e1;
                            if r.q {
                                s.3 += e0 * e0 - e1 * e1;
                            }
                            s.4 += normal_logpdf(y, predict(&b1), v1)
                                - normal_logpdf(y, predict(&b0), v0);
                        }
                        for (donor, (n, nq, mse, qmse, lp)) in scores {
                            writeln!(
                                w,
                                "{acc},{split},{penalty},{donor},{n},{nq},{},{},{},{}",
                                mse / n as f64,
                                if nq > 0 { qmse / nq as f64 } else { f64::NAN },
                                lp / n as f64,
                                b1[p - 1]
                            )?;
                        }
                    }
                }
            }
        }
        // Matching reference is computed from training donors only; held-out outcomes never set strata.
        for donor in &donors {
            let mut early: BTreeMap<String, (f64, usize)> = BTreeMap::new();
            for other in donors.iter().filter(|d| *d != donor) {
                let mut within: BTreeMap<&str, (f64, usize)> = BTreeMap::new();
                for r in rows.iter().filter(|r| &r.donor == other) {
                    let s = within.entry(&r.cpg).or_default();
                    s.0 += r.early;
                    s.1 += 1;
                }
                for (cpg, (sum, n)) in within {
                    let s = early.entry(cpg.into()).or_default();
                    s.0 += sum / n as f64;
                    s.1 += 1;
                }
            }
            for matching in [
                "coarse",
                "standard",
                "fine",
                "no_sequence",
                "no_density",
                "no_context",
            ] {
                let mut groups: BTreeMap<String, BTreeMap<String, (bool, f64, usize)>> =
                    BTreeMap::new();
                for r in rows.iter().filter(|r| &r.donor == donor) {
                    let Some((sum, n)) = early.get(&r.cpg) else {
                        continue;
                    };
                    let (beta_bins, seq_bins, density_bins) = match matching {
                        "coarse" => (10.0, 5.0, 50.0),
                        "fine" => (40.0, 20.0, 200.0),
                        _ => (20.0, 10.0, 100.0),
                    };
                    let context = if matching == "no_context" {
                        String::new()
                    } else {
                        format!("{}:{}", r.island, r.regulatory)
                    };
                    let key = format!(
                        "{}:{}:{}:{}:{}",
                        r.chrom,
                        context,
                        (sum / *n as f64 * beta_bins).floor(),
                        if matching == "no_sequence" {
                            0.0
                        } else {
                            (r.sequence * seq_bins).floor()
                        },
                        if matching == "no_density" {
                            0.0
                        } else {
                            (r.density * density_bins).floor()
                        }
                    );
                    let s = groups
                        .entry(key)
                        .or_default()
                        .entry(r.cpg.clone())
                        .or_insert((r.q, 0.0, 0));
                    s.1 += r.change.abs();
                    s.2 += 1;
                }
                let (mut contrast, mut weight) = (0.0, 0usize);
                for g in groups.values() {
                    let q: Vec<_> = g
                        .values()
                        .filter(|r| r.0)
                        .map(|r| r.1 / r.2 as f64)
                        .collect();
                    let c: Vec<_> = g
                        .values()
                        .filter(|r| !r.0)
                        .map(|r| r.1 / r.2 as f64)
                        .collect();
                    if !q.is_empty() && c.len() >= 2 {
                        contrast += q.len() as f64 * (mean(&q) - mean(&c));
                        weight += q.len();
                    }
                }
                writeln!(
                    mw,
                    "{acc},{matching},{donor},{weight},{}",
                    if weight > 0 {
                        contrast / weight as f64
                    } else {
                        f64::NAN
                    }
                )?;
            }
        }
    }
    Ok(())
}

// Restricted bulk M1/M2 test at externally unmethylated ESLs and low-baseline controls.
// Common gain/return rates are a deliberately simple pooling assumption.
fn kinetic_followup(out: &str) -> Result<()> {
    let features = load_features()?;
    let mut w = writer(&format!("{out}/pooled_kinetic_comparison.csv"))?;
    writeln!(w,"dataset,fold,held_out_donor,n,n_esl,model,log_predictive_gain_vs_m1,mse_gain_vs_m1,mse_gain_vs_constant,log_gain_rate,log_return_rate,protection_effect,restoration_effect")?;
    for acc in ["GSE179847", "GSE73115"] {
        let rows: Vec<_> = forecast_rows(acc, out, &features)?
            .into_iter()
            .filter(|r| r.early <= 0.2)
            .collect();
        let donors: Vec<_> = rows
            .iter()
            .map(|r| r.donor.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let counts: BTreeMap<_, _> = donors
            .iter()
            .map(|d| (d.clone(), rows.iter().filter(|r| &r.donor == d).count()))
            .collect();
        let prediction = |r: &ForecastPair, p: &[f64]| {
            probability(
                r.early,
                (p[0] - p.get(2).copied().unwrap_or(0.0) * f64::from(r.q)).exp(),
                (p[1] + p.get(3).copied().unwrap_or(0.0) * f64::from(r.q)).exp(),
                r.x[16],
            )
        };
        for fold in 0..5 {
            let test_fold = |r: &ForecastPair| donors.binary_search(&r.donor).unwrap() % 5 == fold;
            let train: Vec<_> = rows.iter().filter(|r| !test_fold(r)).collect();
            let sse = |p: &[f64]| {
                let (mut sum, mut weight) = (0.0, 0.0);
                for r in &train {
                    let wt = 1.0 / counts[&r.donor] as f64;
                    sum += wt * (r.early + r.change - prediction(r, p)).powi(2);
                    weight += wt;
                }
                sum / weight
            };
            let (_, m1) = optimize(
                sse,
                &[vec![-4.0, -2.0], vec![-2.0, -1.0], vec![-6.0, 0.0]],
                &[(-10.0, 4.0); 2],
            );
            let v1 = sse(&m1).max(1e-6);
            for (model, bound) in [
                ("m1", (0.0, 0.0)),
                ("m2_protection", (0.0, 6.0)),
                ("m2_signed", (-6.0, 6.0)),
            ] {
                let p = if model == "m1" {
                    m1.clone()
                } else {
                    let mut start = m1.clone();
                    start.extend([0.0, 0.0]);
                    let mut alt = start.clone();
                    alt[2] = 1.0;
                    alt[3] = 1.0;
                    optimize(
                        sse,
                        &[start, alt],
                        &[(-10.0, 4.0), (-10.0, 4.0), bound, bound],
                    )
                    .1
                };
                let v = sse(&p).max(1e-6);
                let mut scores: BTreeMap<String, (usize, usize, f64, f64, f64)> = BTreeMap::new();
                for r in rows.iter().filter(|r| test_fold(r)) {
                    let y = r.early + r.change;
                    let mu = prediction(r, &p);
                    let mu1 = prediction(r, &m1);
                    let s = scores.entry(r.donor.clone()).or_default();
                    s.0 += 1;
                    s.1 += usize::from(r.q);
                    s.2 += normal_logpdf(y, mu, v) - normal_logpdf(y, mu1, v1);
                    s.3 += (y - mu1).powi(2) - (y - mu).powi(2);
                    s.4 += r.change.powi(2) - (y - mu).powi(2);
                }
                for (donor, (n, nq, lp, mse, constant)) in scores {
                    writeln!(
                        w,
                        "{acc},{fold},{donor},{n},{nq},{model},{},{},{},{},{},{},{}",
                        lp / n as f64,
                        mse / n as f64,
                        constant / n as f64,
                        p[0],
                        p[1],
                        p.get(2).copied().unwrap_or(0.0),
                        p.get(3).copied().unwrap_or(0.0)
                    )?;
                }
            }
        }
        println!(
            "{acc}: pooled M1/M2 scored {} low-baseline pairs",
            rows.len()
        );
    }
    Ok(())
}

// Independent measured proliferation outcome; ESL is a stability mask, not an identity mask.
fn growth_followup(root: &str, out: &str) -> Result<()> {
    let features = load_features()?;
    let rows = forecast_rows("GSE179847", out, &features)?;
    let mut attributes: Vec<BTreeMap<String, String>> = vec![];
    let mut current = BTreeMap::new();
    for line in reader(&format!("{root}/GSE179847_family.soft.gz"))?.lines() {
        let line = line?;
        if line.starts_with("^SAMPLE = ") {
            if !current.is_empty() {
                attributes.push(std::mem::take(&mut current));
            }
        } else if let Some(v) = line.strip_prefix("!Sample_characteristics_ch1 = ") {
            if let Some((k, v)) = v.split_once(": ") {
                current.insert(k.to_lowercase(), v.into());
            }
        }
    }
    attributes.push(current);
    type GrowthEndpoints = (String, f64, f64, Vec<f64>, Vec<f64>);
    let mut pairs: BTreeMap<String, GrowthEndpoints> = BTreeMap::new();
    let mut lines = reader(&format!("{out}/GSE179847_forecast_pairs.tsv"))?.lines();
    lines.next();
    for line in lines {
        let a = fields(&line?, '\t');
        pairs.entry(a[2].clone()).or_insert((
            a[1].clone(),
            a[3].parse()?,
            a[4].parse()?,
            vec![],
            vec![],
        ));
    }
    // Use exactly the same annotated site universe as the forecast, for both distances.
    let mut lines = reader(&format!("{out}/GSE179847_forecast_pairs.tsv"))?.lines();
    lines.next();
    let eligible: BTreeSet<_> = rows.iter().map(|r| r.cpg.clone()).collect();
    for line in lines {
        let a = fields(&line?, '\t');
        if !eligible.contains(&a[0]) {
            continue;
        }
        let d = (a[6].parse::<f64>()? - a[5].parse::<f64>()?).abs();
        let pair = pairs.get_mut(&a[2]).unwrap();
        pair.3.push(d);
        if features[&a[0]].esl {
            pair.4.push(d);
        }
    }
    let mut scored: Vec<(String, String, Vec<f64>, f64)> = vec![];
    let mut w = writer(&format!("{out}/growth_pairs.csv"))?;
    writeln!(w,"donor,lineage,duration,early_doubling_hours,late_doubling_hours,log_doubling_ratio,unweighted_distance,esl_distance,n_sites,n_esl")?;
    for (lineage, (donor, start, duration, all, q)) in pairs {
        let endpoint = |time: f64| {
            attributes.iter().find_map(|a| {
                if a.get("cell_line_group")? != &lineage
                    || (numeric(a.get("days_grown_udays")?)? - time).abs() > 1e-5
                {
                    return None;
                }
                numeric(a.get("population_doubling_time_uhours_per_division")?).filter(|v| *v > 0.0)
            })
        };
        let (Some(early), Some(late)) = (endpoint(start), endpoint(start + duration)) else {
            continue;
        };
        if q.is_empty() || all.is_empty() {
            continue;
        }
        let y = (late / early).ln();
        let d = mean(&all);
        let qd = mean(&q);
        writeln!(
            w,
            "{donor},{lineage},{duration},{early},{late},{y},{d},{qd},{},{}",
            all.len(),
            q.len()
        )?;
        scored.push((
            donor,
            lineage,
            vec![
                1.0,
                duration / 100.0,
                (early / 48.0).ln(),
                d * 10.0,
                qd * 10.0,
            ],
            y,
        ));
    }
    let donors: BTreeSet<_> = scored.iter().map(|r| r.0.clone()).collect();
    let mut w = writer(&format!("{out}/growth_prediction.csv"))?;
    writeln!(w,"held_out_donor,n,penalty,unweighted_mse_gain_vs_time,esl_mse_gain_vs_time,combined_mse_gain_vs_time")?;
    for donor in &donors {
        let mut train = [
            NormalEquations::new(3),
            NormalEquations::new(4),
            NormalEquations::new(4),
            NormalEquations::new(5),
        ];
        for (d, _, x, y) in scored.iter().filter(|r| &r.0 != donor) {
            let weight = 1.0 / scored.iter().filter(|r| &r.0 == d).count() as f64;
            let designs = [
                x[..3].to_vec(),
                x[..4].to_vec(),
                vec![x[0], x[1], x[2], x[4]],
                x.clone(),
            ];
            for (t, x) in train.iter_mut().zip(designs) {
                t.add_xy(&x, *y, weight);
            }
        }
        for penalty in [0.001, 0.01, 0.1] {
            let b: Vec<_> = train
                .iter()
                .map(|t| t.fit(t.rhs.len(), penalty))
                .collect::<Result<Vec<_>>>()?;
            let mut errors = [0.0; 4];
            let mut n = 0;
            for (_, _, x, y) in scored.iter().filter(|r| &r.0 == donor) {
                let designs = [
                    x[..3].to_vec(),
                    x[..4].to_vec(),
                    vec![x[0], x[1], x[2], x[4]],
                    x.clone(),
                ];
                for ((err, (coef, _)), x) in errors.iter_mut().zip(&b).zip(designs) {
                    *err += (y - coef.iter().zip(x).map(|(a, b)| a * b).sum::<f64>()).powi(2);
                }
                n += 1;
            }
            writeln!(
                w,
                "{donor},{n},{penalty},{},{},{}",
                (errors[0] - errors[1]) / n as f64,
                (errors[0] - errors[2]) / n as f64,
                (errors[0] - errors[3]) / n as f64
            )?;
        }
    }
    println!(
        "Growth follow-up: {} culture endpoints, {} held-out donors",
        scored.len(),
        donors.len()
    );
    Ok(())
}

fn followup_report(out: &str) -> Result<()> {
    use rand::{Rng, SeedableRng};
    let mut grouped: BTreeMap<String, BTreeMap<String, (f64, f64)>> = BTreeMap::new();
    let mut add = |key: String, donor: String, value: f64, weight: f64| {
        if value.is_finite() && weight > 0.0 {
            let s = grouped.entry(key).or_default().entry(donor).or_default();
            s.0 += value * weight;
            s.1 += weight;
        }
    };
    for line in reader(&format!("{out}/forecast_comparison.csv"))?
        .lines()
        .skip(1)
    {
        let a = fields(&line?, ',');
        add(
            format!("{} forecast {} penalty={} MSE gain", a[0], a[1], a[2]),
            a[3].clone(),
            a[6].parse()?,
            a[4].parse()?,
        );
        add(
            format!("{} forecast {} penalty={} ESL MSE gain", a[0], a[1], a[2]),
            a[3].clone(),
            a[7].parse()?,
            a[5].parse()?,
        );
    }
    for line in reader(&format!("{out}/pooled_kinetic_comparison.csv"))?
        .lines()
        .skip(1)
    {
        let a = fields(&line?, ',');
        add(
            format!("{} {} log density gain vs M1", a[0], a[5]),
            a[2].clone(),
            a[6].parse()?,
            1.0,
        );
        add(
            format!("{} {} MSE gain vs constant", a[0], a[5]),
            a[2].clone(),
            a[8].parse()?,
            1.0,
        );
    }
    for line in reader(&format!("{out}/growth_prediction.csv"))?
        .lines()
        .skip(1)
    {
        let a = fields(&line?, ',');
        for (i, name) in [(3, "unweighted"), (4, "ESL"), (5, "combined")] {
            add(
                format!("growth {name} penalty={} MSE gain", a[2]),
                a[0].clone(),
                a[i].parse()?,
                1.0,
            );
        }
    }
    for acc in ["GSE179847", "GSE73115"] {
        for matching in ["standard", "coarse", "fine", "no_density"] {
            let file = if matching == "standard" {
                format!("{out}/{acc}_matched_replicates.csv")
            } else {
                format!("{out}/{acc}_{matching}_matched_replicates.csv")
            };
            for line in reader(&file)?.lines().skip(1) {
                let a = fields(&line?, ',');
                add(
                    format!("{acc} full matching {matching} ESL-minus-control"),
                    a[0].clone(),
                    a[2].parse()?,
                    1.0,
                );
            }
        }
    }
    let mut w = writer(&format!("{out}/followup_summary.csv"))?;
    writeln!(w,"comparison,biological_units,mean,bootstrap_lower,bootstrap_upper,units_with_positive_value,status")?;
    let mut md = writer(&format!("{out}/validation_followup.md"))?;
    writeln!(md,"# Continued validation (2026-10-03)\n\nThese exploratory tests strengthen the external **stability** signal, but do not validate essentiality, allocated maintenance, selection or biological identity loss. All gains compare predictions in excluded donors; positive MSE gain means lower error. Blood participants are twins: family IDs are unavailable in GEO, so person-held-out splits and person-bootstrap intervals are provisional.\n\n| Test | Units | Mean | 95% bootstrap |\n|---|---:|---:|---:|")?;
    let mut rng = rand::rngs::StdRng::seed_from_u64(20261003);
    for (key, units) in &grouped {
        let values: Vec<_> = units.values().map(|(sum, n)| sum / n).collect();
        let mut draws: Vec<_> = (0..4000)
            .map(|_| {
                (0..values.len())
                    .map(|_| values[rng.gen_range(0..values.len())])
                    .sum::<f64>()
                    / values.len() as f64
            })
            .collect();
        draws.sort_by(f64::total_cmp);
        let status = if key.starts_with("GSE73115") {
            "person_bootstrap_does_not_account_for_twins"
        } else {
            "fixed_fold_bootstrap_does_not_refit_training"
        };
        writeln!(
            w,
            "{key},{},{},{},{},{},{status}",
            values.len(),
            mean(&values),
            draws[100],
            draws[3899],
            values.iter().filter(|v| **v > 0.0).count()
        )?;
        if key.contains("penalty=0.01") && !key.contains("ESL MSE")
            || key.contains("full matching")
            || key.contains("log density") && !key.contains(" m1 ")
        {
            writeln!(
                md,
                "| {key} | {} | {:.6} | [{:.6}, {:.6}] |",
                values.len(),
                mean(&values),
                draws[100],
                draws[3899]
            )?;
        }
    }
    writeln!(md,"\n## Interpretation and limits\n\n- Forecasts use every 256th input probe, detection filtering, and complete external sequence/context annotations: 1,530 fibroblast probes (85 ESL) and 1,682 blood probes (88 ESL). The response is log(abs(endpoint change)+0.001). The empirical ridge baseline includes starting beta splines, sequence, density, CGI/regulatory class, chromosome, duration/start time and fibroblast condition/oxygen proxies. Its augmented version adds one external ESL indicator. Fixed penalties 0.001/0.01/0.1 are sensitivity settings, not tuned on held-out outcomes. Five chromosome groups and five deterministic donor groups are jointly excluded for the stricter split. Shared CpGs remain in the donor-only test. This is an empirical predictive gate, not a mechanistic M2 likelihood. Sex, array batch, genome-wide density and blood composition remain unadjusted.\n- Pooled kinetic tests use the exact binary CTMC marginal, individual earlier beta as the initial probability, and a Gaussian bulk-array residual. They restrict both ESLs and controls to earlier beta <=0.2, and hold out five donor groups. M1 pools gain/return rates; M2 multiplies gain by exp(-alpha*q) and return by exp(beta*q). Both positive-only and signed effects are fit from two starts after three M1 starts. Rates use 100 days/years as their unit. Site/context heterogeneity is not fully modeled here; q coefficients reaching optimization bounds and protection/restoration tradeoffs forbid parameter-level biological interpretation. Small predictive gains are evidence for this restricted pooling model only. In blood, allowing signed effects fits alpha around -0.93 to -0.98 (higher gain at ESLs) and beta around 2.05 to 2.14 (higher return), and scores better than positive-only protection. Thus the observed stability does not specifically support suppressed away rates. Fibroblast beta repeatedly reaches plus/minus 6; restoration is not identified.\n- Full-matrix matching varies beta, density and sequence bin widths together (coarse/standard/fine), plus omission of density. These runs use cohort early means and the same external context as the initial analysis; they remain transductive. The sampled, donor-training-only matching exercise is recorded separately in matching_sensitivity.csv: sparse strict strata often have zero support. Unsupported comparisons are NaN, not zero effects.\n- Growth is a first independent biological endpoint: log(late/early measured population-doubling time), using the same DNA sample times. Only 11 culture endpoint pairs from eight donors have finite positive doubling measurements at both times. Forecasts compare time/starting-doubling-time covariates against additional unweighted and ESL absolute methylation distances, with donors excluded. This small pilot tests association with proliferation slowing, not identity loss; the ESL distance is not an independently defined functional weight. No threshold, causal mechanism or selection-versus-protection claim follows.\n- Bootstrap intervals resample donor-level fixed test scores and omit uncertainty from refitting and spatial probe dependence. Blood twins require family mapping for valid family-held-out inference. No multiple-comparison correction is claimed; all settings/results are retained.\n\nThe correct paired fibroblast RNA study is GSE179848 (downloaded and checksummed); the previously downloaded GSE225172 expression file belongs to the mouse study. RNA outcomes are not analyzed in this follow-up.\n\nRun `make followup` after `make reproduce` to regenerate these extended tests. Detailed comparisons are in forecast_comparison.csv, pooled_kinetic_comparison.csv, growth_pairs.csv, growth_prediction.csv, full matching CSVs and followup_summary.csv.")?;
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let cmd = args.get(1).map(String::as_str).unwrap_or("reproduce");
    let root = args.get(2).map(String::as_str).unwrap_or("data/raw");
    let out = args.get(3).map(String::as_str).unwrap_or("code/results");
    fs::create_dir_all(format!("{out}/posterior_predictive_checks"))?;
    match cmd {
        "baseline" => {
            let c = cells(root, out)?;
            baseline(&c, out)?;
        }
        "arrays" => {
            arrays(root, out, true, None, "standard")?;
            arrays(root, out, false, None, "standard")?;
        }
        "matching-followup" => {
            for matching in ["coarse","fine","no_density"] {
                let temporary=format!("data/derived/matching_runs/{matching}"); fs::create_dir_all(&temporary)?;
                for (acc,fibro) in [("GSE179847",true),("GSE73115",false)] {
                    arrays(root,&temporary,fibro,None,matching)?;
                    for suffix in ["matched_mask","matched_replicates"] {
                        fs::copy(format!("{temporary}/{acc}_{suffix}.csv"),format!("{out}/{acc}_{matching}_{suffix}.csv"))?;
                    }
                }
            }
        }
        "blood" => arrays(root, out, false, None, "standard")?,
        "annotate" => {annotations(root,out)?;},
        "simulate" => simulation_checks(out)?,
        "report" => reports(out)?,
        "fit" => fit_real_baseline(out)?,
        "forecast" => forecasting(out)?,
        "kinetic-followup" => kinetic_followup(out)?,
        "growth-followup" => growth_followup(root,out)?,
        "followup-report" => followup_report(out)?,
        "reproduce" => {
            annotations(root,out)?;
            let c = cells(root, out)?;
            baseline(&c, out)?;
            arrays(root, out, true, None, "standard")?;
            arrays(root, out, false, None, "standard")?;
            fit_real_baseline(out)?;
            forecasting(out)?;
            simulation_checks(out)?;
            reports(out)?;
        }
        _ => {
            return Err(
                "usage: epidrift [baseline|arrays|blood|fit|forecast|kinetic-followup|growth-followup|matching-followup|followup-report|simulate|report|reproduce] [raw-directory] [results-directory]"
                    .into(),
            )
        }
    }
    if !Path::new(out).exists() {
        return Err("output directory missing".into());
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn regression_recovers_coefficients_and_removes_heldout_rows() {
        let mut all = NormalEquations::new(2);
        let mut heldout = NormalEquations::new(2);
        let mut train = NormalEquations::new(2);
        for i in 0..20 {
            let x = [1.0, i as f64 / 10.0];
            let y = 2.0 - 3.0 * x[1];
            all.add_xy(&x, y, 1.0);
            if i % 2 == 0 {
                heldout.add_xy(&x, y, 1.0);
            } else {
                train.add_xy(&x, y, 1.0);
            }
        }
        let removed = all.subtract(&heldout);
        let (b, v) = removed.fit(2, 0.0).unwrap();
        let (expected, _) = train.fit(2, 0.0).unwrap();
        assert!((b[0] - 2.0).abs() < 1e-8 && (b[1] + 3.0).abs() < 1e-8);
        assert!((b[0] - expected[0]).abs() < 1e-10 && (b[1] - expected[1]).abs() < 1e-10);
        assert!(v <= 1.000001e-6);
    }
    #[test]
    fn csv_quoted() {
        assert_eq!(
            fields("\"a,b\",\"c\"\"d\",1", ','),
            vec!["a,b", "c\"d", "1"]
        );
    }
}
