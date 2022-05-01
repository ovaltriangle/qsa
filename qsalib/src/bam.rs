use std::path::{PathBuf, Path};

use ndarray::{Array, ArrayView, Ix1, Ix2, ShapeBuilder, Axis};
use bam::BamReader;

use crate::matrices::Matrices;
use crate::utils::expand_dir;
use crate::error::{QSAError, Result};

pub struct BamDataBuilder {
    bams: Vec<PathBuf>,
    dirs: Vec<PathBuf>,
    range: (i32, i32),
    threshold: f64,
    checks: bool,
}

impl Default for BamDataBuilder {
    fn default() -> Self {
        BamDataBuilder {
            bams: Vec::default(),
            dirs: Vec::default(),
            range: (i32::default(), i32::default()),
            threshold: f64::default(),
            checks: true,
        }
    }
}

impl BamDataBuilder {
    pub fn add_bam<P>(&mut self, bam: P) -> Result<&mut Self>
        where P: AsRef<Path>
    {
        if bam.as_ref().exists() {
            self.bams.push(bam.as_ref().to_path_buf());
        } else {
            return Err(QSAError::BAMNotFound);
        }

        Ok(self)
    }

    pub fn add_bams<P>(&mut self, bams: Vec<P>) -> Result<&mut Self>
        where P: AsRef<Path>
    {
        for bam in bams {
            self.add_bam(bam)?;
        }

        Ok(self)
    }

    pub fn add_dir<P>(&mut self, dir: P) -> Result<&mut Self>
        where P: AsRef<Path>
    {
        if !dir.as_ref().is_dir() {
            return Err(QSAError::DirNotFound);
        }

        self.dirs.push(dir.as_ref().to_path_buf());

        Ok(self)
    }

    pub fn add_dirs<P>(&mut self, dirs: Vec<P>) -> Result<&mut Self>
        where P: AsRef<Path>
    {
        for dir in dirs {
            self.add_dir(dir)?;
        }

        Ok(self)
    }

    pub fn in_range(&mut self, range: (i32, i32)) -> &mut Self {
        self.range = range;

        self
    }

    pub fn with_threshold(&mut self, threshold: f64) -> &mut Self {
        self.threshold = threshold;

        self
    }

    pub fn with_checks(&mut self, checks: bool) -> &mut Self {
        self.checks = checks;

        self
    }

    pub fn build(&mut self) -> Result<BamData> {
        for dir in &self.dirs {
            self.bams.append(&mut expand_dir(dir.to_str().unwrap(), "bam"));
        }

        let mut bams: Vec<Bam> = Vec::new();
        for bamp in &self.bams {
            let bam = Bam::new(bamp, self.range, self.threshold)?;

            bams.push(bam);
        }

        BamData::from_bams(bams, self.checks)
    }
}

#[derive(Default)]
pub struct BamData {
    bams: Vec<Bam>,
    checks: bool,
    alpha: Array<f64, Ix1>,
    beta: Array<f64, Ix2>,
}

impl BamData {
    fn alpha(bams: &Vec<Bam>) -> Array<f64, Ix1> {
        Array::from_vec(bams.iter().map(|x| x.alpha_diversity()).collect::<Vec<f64>>())
    }

    fn alpha_add(&mut self, bam: &Bam) -> &mut Self {
        self.alpha.append(Axis(0), ArrayView::from(&[bam.alpha_diversity()])).unwrap();

        self
    }

    fn beta(alpha: ArrayView<f64, Ix1>) -> Array<f64, Ix2> {
        let cols = *alpha.shape().get(0).unwrap();
        let mut beta = Array::<f64, Ix2>::zeros((cols, cols).f());

        for i in 0..cols {
            for j in 0..cols {
                unsafe { *beta.uget_mut([i, j]) = (alpha[[i]] - alpha[[j]]).abs(); }
            }
        }

        beta
    }

    fn beta_upd(&mut self) -> &mut Self {
        let (rest, last) = self.alpha.view().split_at(Axis(0), self.alpha.len() - 1);
        let last = last.get(0).unwrap();

        let values = rest.iter()
            .map(|x| (x - last).abs())
            .collect::<Vec<f64>>();

        // CHECK: Does this work?
        self.beta.push_column(ArrayView::from(&values.as_slice()[..values.len() - 1])).unwrap();
        self.beta.push_row(ArrayView::from(values.as_slice())).unwrap();

        self
    }

    pub fn from_bams(bams: Vec<Bam>, checks: bool) -> Result<Self> {
        if checks {
            let checked = bams.iter()
                .map(|x| &x.sqsn)
                .collect::<Vec<_>>()
                .windows(2)
                .all(|w| w[0] == w[1] || w[0].is_empty() || w[1].is_empty());

            if !checked {
                return Err(QSAError::BamChecksFailed);
            }
        }

        let alpha = BamData::alpha(&bams);
        let beta = BamData::beta(alpha.view());

        Ok(
            BamData {
                bams,
                checks,
                alpha,
                beta
            }
        )
    }

    pub fn push(&mut self, bam: Bam) -> Result<()> {
        if self.checks {
            let sqsn = &bam.sqsn;

            if !sqsn.is_empty() {
                let checked = self.bams.iter()
                    .map(|x| &x.sqsn)
                    .any(|x| x == sqsn);

                if !checked {
                    return Err(QSAError::BamChecksFailed);
                }
            }
        }

        self.alpha_add(&bam).beta_upd().bams.push(bam);

        Ok(())
    }

    pub fn alpha_diversity(&self) -> ArrayView<f64, Ix1> {
        self.alpha.view()
    }

    pub fn beta_diversity(&self) -> ArrayView<f64, Ix2> {
        self.beta.view()
    }

    pub fn get_names(&self) -> Vec<String> {
        let mut rv: Vec<String> = Vec::new();
        for bam in &self.bams {
            rv.push(bam.name.clone());
        }

        rv
    }
}

pub struct BamDataIntoIterator {
    bam: std::vec::IntoIter<Bam>,
}

impl IntoIterator for BamData {
    type Item = Bam;
    type IntoIter = BamDataIntoIterator;

    fn into_iter(self) -> Self::IntoIter {
        BamDataIntoIterator {
            bam: self.bams.into_iter(),
        }
    }
}

impl Iterator for BamDataIntoIterator {
    type Item = Bam;

    fn next(&mut self) -> Option<Self::Item> {
        self.bam.next()
    }
}

pub struct BamDataIterator<'a> {
    bam: std::slice::Iter<'a, Bam>,
}

impl<'a> IntoIterator for &'a BamData {
    type Item = &'a Bam;
    type IntoIter = BamDataIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        BamDataIterator {
            bam: self.bams.iter()
        }
    }
}

impl<'a> Iterator for BamDataIterator<'a> {
    type Item = &'a Bam;

    fn next(&mut self) -> Option<Self::Item> {
        self.bam.next()
    }
}

pub struct Bam {
    pub name: String,
    pub matrices: Matrices,
    pub(crate) sqsn: String,
}

impl Bam {
    pub fn new<P>(bam: P, range: (i32, i32), threshold: f64) -> Result<Self>
        where P: AsRef<Path>
    {
        let name = bam.as_ref().iter().nth(1).unwrap().to_str().unwrap();
        let name = name[..name.len() - 4].to_string();

        let bam = BamReader::from_path(bam, 0).unwrap();

        let sqsn =
            if bam.header().n_references() > 0 {
                bam.header().reference_name(0).unwrap().to_string()
            } else {
                "".to_string()
            };

        let matrices = Matrices::new(bam, range, threshold)?;

        Ok(
            Bam {
                name,
                matrices,
                sqsn,
            }
        )
    }

    pub fn alpha_diversity(&self) -> f64 {
        self.matrices.get_efficiency().sum() / self.matrices.get_coverage().len() as f64
    }

    pub fn set_name(&mut self, name: String) -> &mut Self {
        self.name = name;

        self
    }

    pub fn pfm_to_csv<P>(&self, path: P, filename: &str)
        where P: AsRef<Path>
    {
        self.matrices.pfm_to_csv(path, filename);
    }
}