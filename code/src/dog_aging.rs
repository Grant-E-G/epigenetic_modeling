//! Frozen dog-aging proxy test. Bulk methylation cannot identify intrinsic repair or selection.
use super::{fields, reader, writer, Result};
use epidrift::mean;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufRead, Read, Seek, SeekFrom, Write},
};

#[derive(Clone, Debug)]
struct Sample {
    id: String,
    dog: u64,
    age: f64,
    sex: String,
    size: String,
    batch: String,
}
#[derive(Clone, Debug)]
struct Dog {
    first: usize,
    next: usize,
    sample: Sample,
    dt: f64,
    next_batch: String,
}
#[derive(Clone, Debug)]
struct Interval {
    start: u32,
    end: u32,
    gene: u32,
}
#[derive(Default)]
struct Track {
    intervals: Vec<Interval>,
    max_end: Vec<u32>,
}
impl Track {
    fn build(mut intervals: Vec<Interval>) -> Self {
        intervals.sort_by_key(|v| (v.start, v.end, v.gene));
        let mut end = 0;
        let max_end = intervals
            .iter()
            .map(|v| {
                end = end.max(v.end);
                end
            })
            .collect();
        Self { intervals, max_end }
    }
    fn overlaps(&self, start: u32, end: u32) -> Vec<&Interval> {
        let mut i = self.intervals.partition_point(|v| v.start < end);
        let mut out = vec![];
        while i > 0 && self.max_end[i - 1] > start {
            i -= 1;
            if self.intervals[i].end > start {
                out.push(&self.intervals[i]);
            }
        }
        out
    }
}
type Tracks = BTreeMap<String, Track>;
struct Annotation {
    promoters: Tracks,
    bodies: Tracks,
    islands: Tracks,
    functional: BTreeMap<u32, (bool, String)>,
}

struct TwoBit {
    file: File,
    little: bool,
    index: BTreeMap<String, u32>,
    current: String,
    packed: Vec<u8>,
    n_blocks: Vec<(u32, u32)>,
    size: u32,
}
fn word(f: &mut File, little: bool) -> Result<u32> {
    let mut b = [0; 4];
    f.read_exact(&mut b)?;
    Ok(if little {
        u32::from_le_bytes(b)
    } else {
        u32::from_be_bytes(b)
    })
}
impl TwoBit {
    fn open(path: &str) -> Result<Self> {
        let mut file = File::open(path)?;
        let mut sig = [0; 4];
        file.read_exact(&mut sig)?;
        let little = match u32::from_be_bytes(sig) {
            0x1a41_2743 => false,
            0x4327_411a => true,
            _ => return Err("invalid twoBit signature".into()),
        };
        if word(&mut file, little)? != 0 {
            return Err("unsupported twoBit version".into());
        }
        let count = word(&mut file, little)?;
        let _ = word(&mut file, little)?;
        let mut index = BTreeMap::new();
        for _ in 0..count {
            let mut len = [0];
            file.read_exact(&mut len)?;
            let mut name = vec![0; len[0] as usize];
            file.read_exact(&mut name)?;
            index.insert(String::from_utf8(name)?, word(&mut file, little)?);
        }
        Ok(Self {
            file,
            little,
            index,
            current: String::new(),
            packed: vec![],
            n_blocks: vec![],
            size: 0,
        })
    }
    fn load(&mut self, chrom: &str) -> Result<()> {
        if self.current == chrom {
            return Ok(());
        }
        let offset = *self
            .index
            .get(chrom)
            .ok_or("chromosome absent from twoBit")?;
        self.file.seek(SeekFrom::Start(u64::from(offset)))?;
        self.size = word(&mut self.file, self.little)?;
        let n = word(&mut self.file, self.little)?;
        let starts: Vec<_> = (0..n)
            .map(|_| word(&mut self.file, self.little))
            .collect::<Result<_>>()?;
        self.n_blocks.clear();
        for start in starts {
            self.n_blocks
                .push((start, start + word(&mut self.file, self.little)?));
        }
        let masks = word(&mut self.file, self.little)?;
        self.file.seek(SeekFrom::Current(i64::from(masks) * 8))?;
        let _ = word(&mut self.file, self.little)?;
        self.packed.resize(self.size.div_ceil(4) as usize, 0);
        self.file.read_exact(&mut self.packed)?;
        self.current = chrom.to_string();
        Ok(())
    }
    fn base(&self, pos: u32) -> u8 {
        let i = self.n_blocks.partition_point(|v| v.0 <= pos);
        if i > 0 && self.n_blocks[i - 1].1 > pos {
            return b'N';
        }
        b"TCAG"[((self.packed[(pos / 4) as usize] >> (6 - 2 * (pos % 4))) & 3) as usize]
    }
    fn context(&mut self, chrom: &str, start: u32, end: u32) -> Result<Option<(f64, f64)>> {
        self.load(chrom)?;
        if end > self.size || start >= end {
            return Err("invalid reference interval".into());
        }
        let mut gc = 0;
        let mut cpg = 0;
        for p in start..end {
            let b = self.base(p);
            if b == b'N' {
                return Ok(None);
            }
            gc += usize::from(b == b'C' || b == b'G');
            cpg += usize::from(b == b'C' && p + 1 < end && self.base(p + 1) == b'G');
        }
        Ok(Some((
            gc as f64 / f64::from(end - start),
            cpg as f64 / f64::from(end - start),
        )))
    }
}

fn annotations(root: &str, out: &str) -> Result<Annotation> {
    let mut human = BTreeMap::new();
    let mut audit = writer(&format!("{out}/dog_aging_orthology.csv"))?;
    writeln!(audit, "human_gene,human_symbol,class,dog_gene,status")?;
    for (essential, name) in [(true, "BAGEL_CEGv2.txt"), (false, "BAGEL_NEGv1.txt")] {
        for l in reader(&format!("{root}/{name}"))?.lines().skip(1) {
            let l = l?;
            if l.trim().is_empty() {
                continue;
            }
            let a = fields(&l, '\t');
            if a[2].is_empty() {
                writeln!(
                    audit,
                    ",{},{},,missing_human_gene_id",
                    a[0],
                    if essential { "essential" } else { "control" }
                )?;
                continue;
            }
            human.insert(a[2].parse::<u32>()?, (essential, a[0].clone()));
        }
    }
    let mut forward: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    let mut reverse: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for l in reader(&format!("{root}/dog_aging/gene_orthologs.gz"))?.lines() {
        let l = l?;
        if l.starts_with('#') {
            continue;
        }
        let a: Vec<_> = l.split('\t').collect();
        let ids = if a[0] == "9606" && a[3] == "9615" {
            Some((a[1], a[4]))
        } else if a[0] == "9615" && a[3] == "9606" {
            Some((a[4], a[1]))
        } else {
            None
        };
        if let Some((h, d)) = ids {
            let h = h.parse()?;
            let d = d.parse()?;
            forward.entry(h).or_default().insert(d);
            reverse.entry(d).or_default().insert(h);
        }
    }
    let mut functional = BTreeMap::new();
    for (h, (essential, symbol)) in human {
        let ds = forward.get(&h);
        let d = ds.and_then(|v| {
            if v.len() == 1 {
                v.first().copied()
            } else {
                None
            }
        });
        let valid = d.is_some_and(|d| reverse[&d].len() == 1);
        let status = if valid {
            "unique"
        } else if ds.is_none() {
            "absent"
        } else {
            "ambiguous"
        };
        writeln!(
            audit,
            "{h},{symbol},{},{},{status}",
            if essential { "essential" } else { "control" },
            d.map_or(String::new(), |d| d.to_string())
        )?;
        if valid {
            functional.insert(d.unwrap(), (essential, symbol));
        }
    }
    let mut links = BTreeMap::new();
    for l in reader(&format!("{root}/dog_aging/ncbiRefSeqLink.txt.gz"))?.lines() {
        let l = l?;
        let a: Vec<_> = l.split('\t').collect();
        if let Ok(g) = a[6].parse::<u32>() {
            links.insert(a[0].to_string(), g);
        }
    }
    let mut ps: BTreeMap<String, Vec<Interval>> = BTreeMap::new();
    let mut bs: BTreeMap<String, BTreeMap<u32, (u32, u32)>> = BTreeMap::new();
    for l in reader(&format!("{root}/dog_aging/ncbiRefSeq.txt.gz"))?.lines() {
        let l = l?;
        let a: Vec<_> = l.split('\t').collect();
        let Some(&gene) = links.get(a[1]) else {
            return Err("RefSeq transcript lacks gene ID".into());
        };
        let start = a[4].parse::<u32>()?;
        let end = a[5].parse::<u32>()?;
        let tss = if a[3] == "+" { start } else { end };
        ps.entry(a[2].to_string()).or_default().push(Interval {
            start: tss.saturating_sub(1000),
            end: tss + 1000,
            gene,
        });
        bs.entry(a[2].to_string())
            .or_default()
            .entry(gene)
            .and_modify(|v| {
                v.0 = v.0.min(start);
                v.1 = v.1.max(end);
            })
            .or_insert((start, end));
    }
    let promoters = ps.into_iter().map(|(c, v)| (c, Track::build(v))).collect();
    let bodies = bs
        .into_iter()
        .map(|(c, v)| {
            (
                c,
                Track::build(
                    v.into_iter()
                        .map(|(gene, (start, end))| Interval { start, end, gene })
                        .collect(),
                ),
            )
        })
        .collect();
    let mut islands: BTreeMap<String, Vec<Interval>> = BTreeMap::new();
    for l in reader(&format!("{root}/dog_aging/cpgIslandExt.txt.gz"))?.lines() {
        let l = l?;
        let a: Vec<_> = l.split('\t').collect();
        islands.entry(a[1].to_string()).or_default().push(Interval {
            start: a[2].parse()?,
            end: a[3].parse()?,
            gene: 0,
        });
    }
    Ok(Annotation {
        promoters,
        bodies,
        islands: islands
            .into_iter()
            .map(|(c, v)| (c, Track::build(v)))
            .collect(),
        functional,
    })
}

fn sample_info(root: &str, columns: &[String]) -> Result<Vec<Sample>> {
    let mut batches = BTreeMap::new();
    for l in reader("data/derived/dog_aging/author_metadata.csv")?
        .lines()
        .skip(1)
    {
        let a = fields(&l?, ',');
        batches.insert(a[0].clone(), (a[1].parse::<u64>()?, a[5].clone()));
    }
    let mut samples = BTreeMap::new();
    for l in reader(&format!("{root}/dog_aging/GSE306793_sample_info.csv.gz"))?
        .lines()
        .skip(1)
    {
        let a = fields(&l?, ',');
        let dog = a[1]
            .strip_prefix("dog_id: ")
            .ok_or("bad sample title")?
            .split(':')
            .next()
            .unwrap()
            .parse()?;
        let (author_dog, batch) = batches.get(&a[0]).ok_or("missing author sample ID")?;
        if *author_dog != dog {
            return Err("author/GEO dog ID mismatch".into());
        }
        let s = Sample {
            id: a[0].clone(),
            dog,
            age: a[6].parse()?,
            sex: a[7].clone(),
            size: a[8].clone(),
            batch: batch.clone(),
        };
        if samples.insert(s.id.clone(), s).is_some() {
            return Err("duplicate sample ID".into());
        }
    }
    if samples.len() != columns.len() - 1 {
        return Err("metadata/matrix sample count mismatch".into());
    }
    columns[1..]
        .iter()
        .map(|id| {
            samples
                .remove(id)
                .ok_or_else(|| format!("unmatched matrix ID {id}").into())
        })
        .collect()
}
fn dog_pairs(samples: &[Sample], min_age: f64) -> Vec<Dog> {
    let mut by_dog: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, s) in samples.iter().enumerate() {
        by_dog.entry(s.dog).or_default().push(i);
    }
    let mut out = vec![];
    for mut ids in by_dog.into_values() {
        ids.sort_by(|&i, &j| samples[i].age.total_cmp(&samples[j].age));
        if ids.len() < 2 || samples[ids[0]].age < min_age {
            continue;
        }
        let first = ids[0];
        let next = ids[1];
        let dt = samples[next].age - samples[first].age;
        if (0.5..=2.5).contains(&dt) {
            out.push(Dog {
                first,
                next,
                sample: samples[first].clone(),
                dt,
                next_batch: samples[next].batch.clone(),
            });
        }
    }
    out
}

#[derive(Clone)]
struct Region {
    id: String,
    chrom: String,
    start: u32,
    end: u32,
    gene: u32,
    essential: bool,
    arm: String,
    island: bool,
    beta: f64,
    log_cov: f64,
    gc: f64,
    cpg: f64,
    train_n: usize,
    methyl: Vec<u32>,
    coverage: Vec<u32>,
}
fn assignment(a: &Annotation, chrom: &str, start: u32, end: u32) -> Option<(u32, &'static str)> {
    let autosome = chrom.strip_prefix("chr")?.parse::<u32>().ok()?;
    if !(1..=38).contains(&autosome) {
        return None;
    }
    let ps = a.promoters.get(chrom)?.overlaps(start, end);
    let bs = a.bodies.get(chrom)?.overlaps(start, end);
    let all: BTreeSet<_> = ps.iter().chain(bs.iter()).map(|v| v.gene).collect();
    if all.len() != 1 {
        return None;
    }
    let gene = *all.first()?;
    if !a.functional.contains_key(&gene) {
        return None;
    }
    let arm = if ps.iter().any(|v| v.start <= start && v.end >= end) {
        "promoter"
    } else if ps.is_empty() && bs.iter().any(|v| v.start <= start && v.end >= end) {
        "body"
    } else {
        return None;
    };
    Some((gene, arm))
}

fn read_matrices(
    root: &str,
    out: &str,
    a: &Annotation,
    coverage_gate: u32,
    min_age: f64,
    flank_context: bool,
) -> Result<(Vec<Region>, Vec<Sample>)> {
    let mut m = reader(&format!("{root}/dog_aging/GSE306792_methylation.csv.gz"))?;
    let mut n = reader(&format!("{root}/dog_aging/GSE306792_coverage.csv.gz"))?;
    let mut ml = String::new();
    let mut nl = String::new();
    m.read_line(&mut ml)?;
    n.read_line(&mut nl)?;
    let header = fields(&ml, ',');
    if header != fields(&nl, ',') {
        return Err("matrix headers differ".into());
    }
    let samples = sample_info(root, &header)?;
    let training: Vec<_> = dog_pairs(&samples, min_age)
        .into_iter()
        .filter(|d| d.sample.dog % 5 != 0)
        .map(|d| d.first)
        .collect();
    let mut genome = TwoBit::open(&format!("{root}/dog_aging/canFam4.2bit"))?;
    let mut regions = vec![];
    let mut rows = 0;
    let mut candidate = 0;
    let mut cpg0 = 0;
    let mut cpg1 = 0;
    let mut coord_n = 0;
    let mut zero = 0u64;
    let mut calls = 0u64;
    let mut aggregate_rows = 0;
    let mut seen = BTreeSet::new();
    loop {
        ml.clear();
        nl.clear();
        let mr = m.read_line(&mut ml)?;
        let nr = n.read_line(&mut nl)?;
        if mr == 0 && nr == 0 {
            break;
        }
        if mr == 0 || nr == 0 {
            return Err("matrix row counts differ".into());
        }
        rows += 1;
        let mut mi = ml.trim_end().split(',');
        let mut ni = nl.trim_end().split(',');
        let id = mi.next().unwrap().trim_matches('"');
        if !seen.insert(id.to_string()) {
            return Err("duplicate regional identifier".into());
        }
        if ni.next().unwrap().trim_matches('"') != id {
            return Err("matrix row IDs differ".into());
        }
        let v: Vec<_> = id.split('_').collect();
        let (start, end) = if v.len() == 3 {
            (v[1].parse::<u32>()?, v[2].parse::<u32>()?)
        } else {
            // Supplied promoter/TE aggregates reuse CpGs from the simple region universe.
            // Validate their counts, but never count overlapping aggregate rows as new loci.
            aggregate_rows += 1;
            (0, 0)
        };
        // RegionFinder is an R genomic-range workflow. Verify the convention against CpG sequence.
        if v.len() == 3 && coord_n < 2000 {
            genome.load(v[0])?;
            if start == 0 || start + 1 >= genome.size {
                return Err("bad coordinate audit interval".into());
            }
            cpg0 += usize::from(genome.base(start) == b'C' && genome.base(start + 1) == b'G');
            cpg1 += usize::from(genome.base(start - 1) == b'C' && genome.base(start) == b'G');
            coord_n += 1;
            if coord_n == 2000 && (cpg1 < 1900 || cpg1 <= cpg0 * 5) {
                return Err("one-based CpG coordinates not verified".into());
            }
        }
        let ass = if v.len() == 3 {
            assignment(a, v[0], start.saturating_sub(1), end)
        } else {
            None
        };
        let mut methyl = vec![];
        let mut coverage = vec![];
        for j in 0..samples.len() {
            let m = mi.next().ok_or("short methylation row")?.parse::<u32>()?;
            let n = ni.next().ok_or("short coverage row")?.parse::<u32>()?;
            if m > n {
                return Err(format!("methylated > total at {id},column {j}").into());
            }
            zero += u64::from(n == 0);
            calls += 1;
            if ass.is_some() {
                methyl.push(m);
                coverage.push(n);
            }
        }
        if mi.next().is_some() || ni.next().is_some() {
            return Err("extra matrix values".into());
        }
        let Some((gene, arm)) = ass else {
            continue;
        };
        candidate += 1;
        let baseline: Vec<_> = training
            .iter()
            .copied()
            .filter(|&i| coverage[i] >= coverage_gate)
            .collect();
        if baseline.len() < 20 {
            continue;
        }
        let beta = baseline
            .iter()
            .map(|&i| f64::from(methyl[i]) / f64::from(coverage[i]))
            .sum::<f64>()
            / baseline.len() as f64;
        if (arm == "promoter" && beta > 0.3) || (arm == "body" && beta < 0.7) {
            continue;
        }
        let context_start = if flank_context {
            (start - 1).saturating_sub(250)
        } else {
            start - 1
        };
        genome.load(v[0])?;
        let context_end = if flank_context {
            (end + 250).min(genome.size)
        } else {
            end
        };
        let Some((gc, cpg)) = genome.context(v[0], context_start, context_end)? else {
            continue;
        };
        let log_cov = baseline
            .iter()
            .map(|&i| f64::from(coverage[i]).ln())
            .sum::<f64>()
            / baseline.len() as f64;
        let island = a
            .islands
            .get(v[0])
            .is_some_and(|t| !t.overlaps(start - 1, end).is_empty());
        regions.push(Region {
            id: id.to_string(),
            chrom: v[0].to_string(),
            start: start - 1,
            end,
            gene,
            essential: a.functional[&gene].0,
            arm: arm.to_string(),
            island,
            beta,
            log_cov,
            gc,
            cpg,
            train_n: baseline.len(),
            methyl,
            coverage,
        });
        if rows % 25000 == 0 {
            eprintln!(
                "dog audit: {rows} rows, {} eligible proxy regions",
                regions.len()
            );
        }
    }
    let mut w = writer(&format!(
        "{out}/dog_aging_matrix_audit_{coverage_gate}_{min_age}.csv"
    ))?;
    writeln!(w, "metric,value")?;
    for (key, value) in [
        ("samples", samples.len() as u64),
        ("matrix_regions", rows),
        ("excluded_labeled_aggregate_rows", aggregate_rows),
        ("annotated_proxy_candidates", candidate),
        ("baseline_eligible_regions", regions.len() as u64),
        ("zero_coverage_calls", zero),
        ("matrix_calls", calls),
        ("coordinate_audit_regions", coord_n as u64),
        ("zero_based_start_CpG", cpg0 as u64),
        ("one_based_start_CpG", cpg1 as u64),
    ] {
        writeln!(w, "{key},{value}")?;
    }
    Ok((regions, samples))
}

fn match_distance(e: &Region, c: &Region) -> Option<f64> {
    if e.chrom != c.chrom || e.arm != c.arm || e.island != c.island {
        return None;
    }
    let db = (e.beta - c.beta).abs();
    let dn = (e.log_cov - c.log_cov).abs();
    let dg = (e.gc - c.gc).abs();
    let dp = (e.cpg - c.cpg).abs();
    let dl = (f64::from(e.end - e.start) / f64::from(c.end - c.start))
        .ln()
        .abs();
    if db > 0.05 || dn > 0.5 || dg > 0.05 || dp > 0.02 || dl > 2f64.ln() {
        return None;
    }
    Some(
        (db / 0.05).powi(2)
            + (dn / 0.5).powi(2)
            + (dg / 0.05).powi(2)
            + (dp / 0.02).powi(2)
            + (dl / 2f64.ln()).powi(2),
    )
}
fn matched_regions(regions: &[Region], arm: &str) -> Vec<(usize, usize, f64)> {
    let es: Vec<_> = regions
        .iter()
        .enumerate()
        .filter(|(_, r)| r.essential && r.arm == arm)
        .collect();
    let cs: Vec<_> = regions
        .iter()
        .enumerate()
        .filter(|(_, r)| !r.essential && r.arm == arm)
        .collect();
    let mut candidates = vec![];
    for (i, e) in es {
        for &(j, c) in &cs {
            if let Some(d) = match_distance(e, c) {
                candidates.push((i, j, d));
            }
        }
    }
    candidates.sort_by(|a, b| {
        a.2.total_cmp(&b.2)
            .then_with(|| regions[a.0].id.cmp(&regions[b.0].id))
            .then_with(|| regions[a.1].id.cmp(&regions[b.1].id))
    });
    let mut used = BTreeSet::new();
    let mut out = vec![];
    for (i, j, d) in candidates {
        if used.contains(&regions[i].gene) || used.contains(&regions[j].gene) {
            continue;
        }
        used.insert(regions[i].gene);
        used.insert(regions[j].gene);
        out.push((i, j, d));
    }
    out
}

#[derive(Clone, Copy, Default)]
struct Change {
    adjusted: f64,
    raw: f64,
    signed: f64,
    noise: f64,
    cov: f64,
}
fn change(m0: u32, n0: u32, m1: u32, n1: u32, dt: f64) -> Change {
    let p0 = f64::from(m0) / f64::from(n0);
    let p1 = f64::from(m1) / f64::from(n1);
    let raw = ((p1 - p0) / dt).powi(2);
    let noise =
        (p0 * (1. - p0) / f64::from(n0 - 1) + p1 * (1. - p1) / f64::from(n1 - 1)) / dt.powi(2);
    Change {
        adjusted: raw - noise,
        raw,
        signed: (p1 - p0) / dt,
        noise,
        cov: (f64::from(n0).ln() + f64::from(n1).ln()) / 2.,
    }
}
#[derive(Clone)]
struct Observation {
    dog: usize,
    pair: usize,
    e: Change,
    c: Change,
}
fn observations(
    regions: &[Region],
    matches: &[(usize, usize, f64)],
    dogs: &[Dog],
    gate: u32,
) -> Vec<Observation> {
    let mut out = vec![];
    for (dog, d) in dogs.iter().enumerate() {
        for (pair, &(i, j, _)) in matches.iter().enumerate() {
            let e = &regions[i];
            let c = &regions[j];
            if [
                e.coverage[d.first],
                e.coverage[d.next],
                c.coverage[d.first],
                c.coverage[d.next],
            ]
            .iter()
            .any(|&n| n < gate)
            {
                continue;
            }
            out.push(Observation {
                dog,
                pair,
                e: change(
                    e.methyl[d.first],
                    e.coverage[d.first],
                    e.methyl[d.next],
                    e.coverage[d.next],
                    d.dt,
                ),
                c: change(
                    c.methyl[d.first],
                    c.coverage[d.first],
                    c.methyl[d.next],
                    c.coverage[d.next],
                    d.dt,
                ),
            });
        }
    }
    out
}
fn dog_means(obs: &[Observation], n: usize, raw: bool) -> Vec<(f64, f64, usize)> {
    let mut sums = vec![(0., 0., 0); n];
    for o in obs {
        sums[o.dog].0 += if raw { o.e.raw } else { o.e.adjusted };
        sums[o.dog].1 += if raw { o.c.raw } else { o.c.adjusted };
        sums[o.dog].2 += 1;
    }
    sums.into_iter()
        .filter(|s| s.2 > 0)
        .map(|(e, c, n)| (e / n as f64, c / n as f64, n))
        .collect()
}
fn ratio(e: f64, c: f64) -> f64 {
    if c > 0. {
        1. - e / c
    } else {
        f64::NAN
    }
}
fn percentile(mut x: Vec<f64>, p: f64) -> f64 {
    x.retain(|v| v.is_finite());
    x.sort_by(f64::total_cmp);
    if x.is_empty() {
        return f64::NAN;
    }
    x[((x.len() - 1) as f64 * p).round() as usize]
}
fn bootstrap(
    obs: &[Observation],
    ndogs: usize,
    npairs: usize,
    mode: &str,
    raw: bool,
    rng: &mut StdRng,
) -> (f64, f64, f64, f64, usize) {
    let summaries = dog_means(obs, ndogs, raw);
    let mut draws = Vec::with_capacity(10000);
    for _ in 0..10000 {
        let value = if mode == "dog" {
            let mut e = 0.;
            let mut c = 0.;
            for _ in 0..summaries.len() {
                let row = summaries[rng.gen_range(0..summaries.len())];
                e += row.0;
                c += row.1;
            }
            ratio(e, c)
        } else {
            let mut pw = vec![0.; npairs];
            for _ in 0..npairs {
                pw[rng.gen_range(0..npairs)] += 1.;
            }
            let mut dw = vec![0.; ndogs];
            if mode == "crossed" {
                for _ in 0..ndogs {
                    dw[rng.gen_range(0..ndogs)] += 1.;
                }
            } else {
                dw.fill(1.);
            }
            let mut rows = vec![(0., 0., 0.); ndogs];
            for o in obs {
                rows[o.dog].0 += pw[o.pair] * if raw { o.e.raw } else { o.e.adjusted };
                rows[o.dog].1 += pw[o.pair] * if raw { o.c.raw } else { o.c.adjusted };
                rows[o.dog].2 += pw[o.pair];
            }
            let mut e = 0.;
            let mut c = 0.;
            for (i, (es, cs, w)) in rows.into_iter().enumerate() {
                if w > 0. {
                    e += dw[i] * es / w;
                    c += dw[i] * cs / w;
                }
            }
            ratio(e, c)
        };
        draws.push(value);
    }
    let valid = draws.iter().filter(|v| v.is_finite()).count();
    (
        percentile(draws.clone(), 0.025),
        percentile(draws.clone(), 0.975),
        percentile(draws.clone(), 0.0125),
        percentile(draws, 0.9875),
        valid,
    )
}

pub fn run(root: &str, out: &str, part: &str) -> Result<()> {
    fs::create_dir_all("data/derived/dog_aging")?;
    let a = annotations(root, out)?;
    let settings = match part {
        "all" => vec![(20, 2.), (10, 2.), (30, 2.), (20, 1.), (20, 4.)],
        "primary" | "interval" => vec![(20, 2.)],
        _ => return Err("dog-aging expects primary, interval or all".into()),
    };
    let mut summary = writer(&format!("{out}/dog_aging_summary.csv"))?;
    writeln!(summary,"coverage,min_age,arm,endpoint,resampling,matched_genes,evaluation_dogs,observations,essential_mean,control_mean,protection_fraction,ci95_low,ci95_high,ci97_5_low,ci97_5_high,valid_draws,gate")?;
    let mut rng = StdRng::seed_from_u64(20261009);
    for (coverage_gate, min_age) in settings {
        eprintln!("dog audit: coverage {coverage_gate}, minimum age {min_age}");
        let (regions, samples) =
            read_matrices(root, out, &a, coverage_gate, min_age, part != "interval")?;
        let dogs: Vec<_> = dog_pairs(&samples, min_age)
            .into_iter()
            .filter(|d| d.sample.dog % 5 == 0)
            .collect();
        let mut dw = writer(&format!(
            "{out}/dog_aging_dogs_{coverage_gate}_{min_age}.csv"
        ))?;
        writeln!(
            dw,
            "dog,first_id,next_id,initial_age,elapsed_years,sex,size,first_batch,next_batch"
        )?;
        for d in &dogs {
            writeln!(
                dw,
                "{},{},{},{},{},{},{},{},{}",
                d.sample.dog,
                d.sample.id,
                samples[d.next].id,
                d.sample.age,
                d.dt,
                d.sample.sex,
                d.sample.size,
                d.sample.batch,
                d.next_batch
            )?;
        }
        for arm in ["promoter", "body"] {
            let matches = matched_regions(&regions, arm);
            let mut mw = writer(&format!(
                "{out}/dog_aging_matches_{coverage_gate}_{min_age}_{arm}.csv"
            ))?;
            writeln!(mw,"essential_gene,control_gene,essential_symbol,control_symbol,essential_region,control_region,chromosome,island,essential_beta,control_beta,essential_log_coverage,control_log_coverage,essential_length,control_length,essential_gc,control_gc,essential_cpg,control_cpg,essential_training_dogs,control_training_dogs,distance")?;
            for &(i, j, d) in &matches {
                let e = &regions[i];
                let c = &regions[j];
                writeln!(
                    mw,
                    "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    e.gene,
                    c.gene,
                    a.functional[&e.gene].1,
                    a.functional[&c.gene].1,
                    e.id,
                    c.id,
                    e.chrom,
                    e.island,
                    e.beta,
                    c.beta,
                    e.log_cov,
                    c.log_cov,
                    e.end - e.start,
                    c.end - c.start,
                    e.gc,
                    c.gc,
                    e.cpg,
                    c.cpg,
                    e.train_n,
                    c.train_n,
                    d
                )?;
            }
            let obs = observations(&regions, &matches, &dogs, coverage_gate);
            let cache_prefix = if part == "interval" {
                "interval_observations"
            } else {
                "observations"
            };
            let mut ow = writer(&format!(
                "data/derived/dog_aging/{cache_prefix}_{coverage_gate}_{min_age}_{arm}.csv"
            ))?;
            writeln!(ow,"dog,pair,essential_adjusted,control_adjusted,essential_raw,control_raw,essential_signed,control_signed,essential_noise,control_noise,essential_log_coverage,control_log_coverage")?;
            for o in &obs {
                writeln!(
                    ow,
                    "{},{},{},{},{},{},{},{},{},{},{},{}",
                    dogs[o.dog].sample.dog,
                    o.pair,
                    o.e.adjusted,
                    o.c.adjusted,
                    o.e.raw,
                    o.c.raw,
                    o.e.signed,
                    o.c.signed,
                    o.e.noise,
                    o.c.noise,
                    o.e.cov,
                    o.c.cov
                )?;
            }
            for raw in [false, true] {
                let sums = dog_means(&obs, dogs.len(), raw);
                let em = mean(&sums.iter().map(|s| s.0).collect::<Vec<_>>());
                let cm = mean(&sums.iter().map(|s| s.1).collect::<Vec<_>>());
                let gate = if matches.len() < 20 || sums.len() < 30 {
                    "insufficient_support"
                } else if cm <= 0. {
                    "nonpositive_control"
                } else {
                    "proxy_only"
                };
                for mode in ["dog", "gene", "crossed"] {
                    let (lo, hi, slo, shi, valid) = if sums.is_empty() || matches.is_empty() {
                        (f64::NAN, f64::NAN, f64::NAN, f64::NAN, 0)
                    } else {
                        bootstrap(&obs, dogs.len(), matches.len(), mode, raw, &mut rng)
                    };
                    writeln!(summary,"{coverage_gate},{min_age},{arm},{},{mode},{},{},{},{em},{cm},{},{lo},{hi},{slo},{shi},{valid},{gate}",if raw {"raw_squared"} else {"binomial_adjusted"},matches.len(),sums.len(),obs.len(),ratio(em,cm))?;
                }
            }
            eprintln!(
                "dog audit: {arm}, {} matched gene pairs, {} observed dog/pair records",
                matches.len(),
                obs.len()
            );
        }
    }
    Ok(())
}

/// Conditional power is for the observed proxy endpoint, not a latent kinetic coefficient.
/// The wild generator keeps dog-level paired residual magnitudes and assumes independent dogs.
fn conditional_power(sums: &[(f64, f64, usize)], rng: &mut StdRng) -> Vec<(f64, f64, f64)> {
    if sums.len() < 2 {
        return vec![];
    }
    let control = mean(&sums.iter().map(|r| r.1).collect::<Vec<_>>());
    if control <= 0. {
        return vec![];
    }
    let delta: Vec<_> = sums.iter().map(|r| r.1 - r.0).collect();
    let center = mean(&delta);
    let residual: Vec<_> = delta.iter().map(|v| v - center).collect();
    let statistic = |effect: f64, rng: &mut StdRng| {
        let y: Vec<_> = residual
            .iter()
            .map(|&r| if rng.gen_bool(0.5) { r } else { -r })
            .map(|r| r + effect * control)
            .collect();
        mean(&y) / (epidrift::variance(&y) / y.len() as f64).sqrt()
    };
    let critical = percentile((0..10000).map(|_| statistic(0., rng)).collect(), 0.9875);
    [0., 0.05, 0.10, 0.20]
        .into_iter()
        .map(|effect| {
            let detected = (0..10000)
                .filter(|_| statistic(effect, rng) > critical)
                .count();
            let p = detected as f64 / 10000.;
            (effect, p, (p * (1. - p) / 10000.).sqrt())
        })
        .collect()
}

pub fn audit(out: &str) -> Result<()> {
    let mut power = writer(&format!("{out}/dog_aging_power.csv"))?;
    writeln!(power,"coverage,min_age,arm,endpoint,dogs,injected_protection,detection_probability,monte_carlo_se,draws,scope")?;
    let mut missing = writer(&format!("{out}/dog_aging_measurement_audit.csv"))?;
    writeln!(missing,"coverage,min_age,arm,evaluation_dogs,observed_dogs,matched_gene_pairs,observed_gene_pairs,paired_records,paired_record_fraction,median_pairs_per_observed_dog,median_dogs_per_observed_pair,same_batch_dogs,missing_batch_dogs,mean_essential_noise,mean_control_noise,mean_essential_log_coverage,mean_control_log_coverage")?;
    let mut balance = writer(&format!("{out}/dog_aging_balance.csv"))?;
    writeln!(balance,"coverage,min_age,arm,covariate,essential_mean,control_mean,mean_difference,pooled_standardized_difference")?;
    let mut subgroup = writer(&format!("{out}/dog_aging_subgroups.csv"))?;
    writeln!(subgroup,"coverage,min_age,arm,group,endpoint,dogs,essential_mean,control_mean,protection_fraction,scope")?;
    let mut rng = StdRng::seed_from_u64(20261010);
    for (gate, min_age) in [(20, 2.), (10, 2.), (30, 2.), (20, 1.), (20, 4.)] {
        let mut dogs = vec![];
        for l in reader(&format!("{out}/dog_aging_dogs_{gate}_{min_age}.csv"))?
            .lines()
            .skip(1)
        {
            dogs.push(fields(&l?, ','));
        }
        let lookup: BTreeMap<_, _> = dogs
            .iter()
            .enumerate()
            .map(|(i, d)| (d[0].parse::<u64>().unwrap(), i))
            .collect();
        for arm in ["promoter", "body"] {
            let mut matches = vec![];
            for l in reader(&format!(
                "{out}/dog_aging_matches_{gate}_{min_age}_{arm}.csv"
            ))?
            .lines()
            .skip(1)
            {
                matches.push(fields(&l?, ','));
            }
            for (label, ei, ci) in [
                ("baseline_beta", 8, 9),
                ("log_coverage", 10, 11),
                ("region_length", 12, 13),
                ("flank_gc", 14, 15),
                ("flank_cpg", 16, 17),
                ("training_coverage_dogs", 18, 19),
            ] {
                let e: Vec<_> = matches
                    .iter()
                    .map(|r| r[ei].parse::<f64>().unwrap())
                    .collect();
                let c: Vec<_> = matches
                    .iter()
                    .map(|r| r[ci].parse::<f64>().unwrap())
                    .collect();
                let em = mean(&e);
                let cm = mean(&c);
                let sd = ((epidrift::variance(&e) + epidrift::variance(&c)) / 2.).sqrt();
                writeln!(
                    balance,
                    "{gate},{min_age},{arm},{label},{em},{cm},{},{}",
                    em - cm,
                    (em - cm) / sd
                )?;
            }
            let mut obs = vec![];
            for l in reader(&format!(
                "data/derived/dog_aging/observations_{gate}_{min_age}_{arm}.csv"
            ))?
            .lines()
            .skip(1)
            {
                let a = fields(&l?, ',');
                obs.push(Observation {
                    dog: lookup[&a[0].parse()?],
                    pair: a[1].parse()?,
                    e: Change {
                        adjusted: a[2].parse()?,
                        raw: a[4].parse()?,
                        signed: a[6].parse()?,
                        noise: a[8].parse()?,
                        cov: a[10].parse()?,
                    },
                    c: Change {
                        adjusted: a[3].parse()?,
                        raw: a[5].parse()?,
                        signed: a[7].parse()?,
                        noise: a[9].parse()?,
                        cov: a[11].parse()?,
                    },
                });
            }
            let mut dc = vec![0; dogs.len()];
            let mut pc = vec![0; matches.len()];
            let mut per_dog = vec![(0., 0., 0., 0., 0); dogs.len()];
            for o in &obs {
                dc[o.dog] += 1;
                pc[o.pair] += 1;
                let v = &mut per_dog[o.dog];
                v.0 += o.e.noise;
                v.1 += o.c.noise;
                v.2 += o.e.cov;
                v.3 += o.c.cov;
                v.4 += 1;
            }
            let observed_dogs = dc.iter().filter(|&&n| n > 0).count();
            let observed_pairs = pc.iter().filter(|&&n| n > 0).count();
            let med_dc = percentile(
                dc.iter().filter(|&&n| n > 0).map(|&n| n as f64).collect(),
                0.5,
            );
            let med_pc = percentile(
                pc.iter().filter(|&&n| n > 0).map(|&n| n as f64).collect(),
                0.5,
            );
            let same_batch = dogs
                .iter()
                .filter(|d| !d[7].is_empty() && d[7] == d[8])
                .count();
            let missing_batch = dogs
                .iter()
                .filter(|d| d[7].is_empty() || d[8].is_empty())
                .count();
            let weights: Vec<_> = per_dog.iter().filter(|r| r.4 > 0).collect();
            let mn = |column: usize| {
                mean(&weights.iter().map(|r|match column {0=>r.0,1=>r.1,2=>r.2,_=>r.3}/r.4 as f64).collect::<Vec<_>>())
            };
            writeln!(missing,"{gate},{min_age},{arm},{},{observed_dogs},{},{observed_pairs},{},{},{med_dc},{med_pc},{same_batch},{missing_batch},{},{},{},{}",dogs.len(),matches.len(),obs.len(),obs.len() as f64/(dogs.len()*matches.len()) as f64,mn(0),mn(1),mn(2),mn(3))?;
            for raw in [false, true] {
                let endpoint = if raw {
                    "raw_squared"
                } else {
                    "binomial_adjusted"
                };
                let sums = dog_means(&obs, dogs.len(), raw);
                for (effect, p, se) in conditional_power(&sums, &mut rng) {
                    writeln!(power,"{gate},{min_age},{arm},{endpoint},{},{effect},{p},{se},10000,conditional_independent_dog_proxy",sums.len())?;
                }
                for group in [
                    "all", "Female", "Male", "Small", "Medium", "Standard", "Large", "Giant",
                ] {
                    let subset: Vec<_> = obs
                        .iter()
                        .filter(|o| {
                            group == "all" || dogs[o.dog][5] == group || dogs[o.dog][6] == group
                        })
                        .cloned()
                        .collect();
                    let s = dog_means(&subset, dogs.len(), raw);
                    let e = mean(&s.iter().map(|r| r.0).collect::<Vec<_>>());
                    let c = mean(&s.iter().map(|r| r.1).collect::<Vec<_>>());
                    writeln!(
                        subgroup,
                        "{gate},{min_age},{arm},{group},{endpoint},{},{e},{c},{},descriptive_only",
                        s.len(),
                        ratio(e, c)
                    )?;
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sampling_correction_is_unbiased_under_independent_binomial_calls() {
        let p0: f64 = 0.25;
        let p1: f64 = 0.5;
        let n = 8u32;
        let prob = |k: u32, p: f64| {
            let choose = (0..k).fold(1., |x, j| x * f64::from(n - j) / f64::from(j + 1));
            choose * p.powi(k as i32) * (1. - p).powi((n - k) as i32)
        };
        let mut expected = 0.;
        for m0 in 0..=n {
            for m1 in 0..=n {
                expected += prob(m0, p0) * prob(m1, p1) * change(m0, n, m1, n, 2.).adjusted;
            }
        }
        assert!((expected - ((p1 - p0) / 2.).powi(2)).abs() < 1e-12);
    }
    #[test]
    fn interval_lookup_keeps_long_overlapping_transcripts() {
        let t = Track::build(vec![
            Interval {
                start: 0,
                end: 1000,
                gene: 1,
            },
            Interval {
                start: 50,
                end: 60,
                gene: 2,
            },
            Interval {
                start: 100,
                end: 110,
                gene: 3,
            },
        ]);
        assert_eq!(
            t.overlaps(200, 300)
                .iter()
                .map(|v| v.gene)
                .collect::<Vec<_>>(),
            vec![1]
        );
        assert!(t.overlaps(1000, 1001).is_empty());
    }
    #[test]
    fn matching_and_measurement_exclude_ineligible_context_and_counts() {
        assert!(change(4, 20, 4, 20, 1.).adjusted < 0.);
        let samples = vec![
            Sample {
                id: "a".into(),
                dog: 5,
                age: 3.,
                sex: "Female".into(),
                size: "Small".into(),
                batch: "b".into(),
            },
            Sample {
                id: "b".into(),
                dog: 5,
                age: 4.,
                sex: "Female".into(),
                size: "Small".into(),
                batch: "b".into(),
            },
        ];
        assert_eq!(dog_pairs(&samples, 2.).len(), 1);
        assert!(dog_pairs(&samples, 4.).is_empty());
    }

    #[test]
    fn future_counts_cannot_choose_matches_and_zero_coverage_is_missing() {
        let region = Region {
            id: "chr1_100_199".into(),
            chrom: "chr1".into(),
            start: 99,
            end: 199,
            gene: 1,
            essential: true,
            arm: "body".into(),
            island: false,
            beta: 0.8,
            log_cov: 4.,
            gc: 0.5,
            cpg: 0.02,
            train_n: 30,
            methyl: vec![18, 18],
            coverage: vec![20, 20],
        };
        let mut control = region.clone();
        control.gene = 2;
        control.essential = false;
        control.id = "chr1_300_399".into();
        let mut regions = vec![region, control];
        let original = matched_regions(&regions, "body");
        regions[0].methyl[1] = 0;
        regions[1].coverage[1] = 0;
        regions[1].methyl[1] = 0;
        assert_eq!(matched_regions(&regions, "body"), original);
        let dog = Dog {
            first: 0,
            next: 1,
            sample: Sample {
                id: "sample".into(),
                dog: 5,
                age: 3.,
                sex: "Female".into(),
                size: "Small".into(),
                batch: "a".into(),
            },
            dt: 1.,
            next_batch: "b".into(),
        };
        assert!(observations(&regions, &original, &[dog], 20).is_empty());
        regions[1].beta = 0.9;
        assert!(matched_regions(&regions, "body").is_empty());
    }
}
