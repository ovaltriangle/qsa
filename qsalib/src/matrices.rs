use std::fs::File;
use std::path::Path;

use ndarray::{Array, ArrayView, Ix1, Ix2, ShapeBuilder, s, Axis};
use bam::{BamReader, Record, RecordReader};
use csv::{Writer, WriterBuilder};

use crate::error::{Result, QSAError};

pub struct Matrices {
    pfm: Array<u64, Ix2>,
    coverage: Array<f64, Ix1>,
    ppm: Array<f64, Ix2>,
    efficiency: Array<f64, Ix1>,
}

impl Matrices {
    fn pfm_coverage(mut bam: BamReader<File>, range: (i32, i32)) -> Result<(Array<u64, Ix2>, Array<f64, Ix1>)> {
        let (start, end) = range;

        let mut pfm = Array::<u64, Ix2>::zeros((4, (end - start) as usize).f());
        let mut coverage = Array::<f64, Ix1>::zeros(((end - start) as usize).f());

        let mut record = Record::new();

        // calculate PFM
        loop {
            match bam.read_into(&mut record) {
                Ok(true) => {
                    let sequence = record.sequence()
                        .to_vec_acgtn_only()
                        .iter()
                        .map(|v| {
                            match v % 32 {
                                1 => 0,     3 => 1,     // a|A  c|C
                                7 => 2,     20 => 3,    // g|G  t|T
                                21 => 3,    _ => 4,     // u|U  n|N
                            }
                        })
                        .collect::<Vec<u8>>();

                    let (seq_start, seq_end) = (record.start(), record.start() + sequence.len() as i32);

                    if seq_start < start || seq_end > end {
                        continue
                    }

                    let fcol = seq_start - start;
                    for (i, row) in sequence.iter().enumerate() {
                        if *row == 4 {
                            continue
                        }

                        let col = fcol as usize + i;

                        let cell = pfm.get_mut((*row as usize, col))
                            .expect(format!("could not access ({}, {}) record", row, col).as_str());
                        *cell += 1;
                    }
                },
                Ok(false) => break,
                Err(why) => panic!("{}", why),
            }
        }

        // calculate coverage
        for col in 0..pfm.ncols() {
            let nt = pfm.column(col).sum() as f64;

            /*
            if nt == 0. {
                return Err(QSAError::CoverageHole);
            }

             */

            *coverage.get_mut(col).unwrap() = nt;
        }

        let max_val = coverage.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        // coverage.map_inplace(|x| *x /= max_val);
        coverage /= max_val;  // broadcast

        // return arrays
        Ok((pfm, coverage))
    }

    fn ppm(pfm: ArrayView<u64, Ix2>) -> Array<f64, Ix2> {
        let mut ppm = pfm.map(|x| *x as f64);

        for col in 0..ppm.ncols() {
            let nt = ppm.column(col).sum();

            ppm.column_mut(col).map_inplace(|x| *x /= nt);
        }

        ppm
    }

    fn efficiency(ppm: ArrayView<f64, Ix2>) -> Array<f64, Ix1> {
        let size = ppm.len_of(Axis(1));
        let mut efficiency = Array::<f64, Ix1>::zeros(size);

        for i in 0..size {
            let col = ppm.column(i);

            let norm_shann = - (col.map(|x| (x * x.log2()) / (4_f64.log2())).sum());

            *efficiency.get_mut(i).unwrap() = norm_shann;
        }

        efficiency
    }

    pub(crate) fn new(bam: BamReader<File>, range: (i32, i32), threshold: f64) -> Result<Matrices> {
        let (pfm, coverage) = Matrices::pfm_coverage(bam, range)?;

        let left_t = coverage.iter().position(|&x| x > threshold).unwrap();
        let right_t = coverage.len() - coverage.iter().rev().position(|&x| x > threshold).unwrap();

        let (pfm, coverage) =
            (
                pfm.slice(s![.., left_t..=right_t]).to_owned(),
                coverage.slice(s![left_t..=right_t]).to_owned(),
            );

        let ppm = Matrices::ppm(pfm.view());
        let efficiency = Matrices::efficiency(ppm.view());

        Ok (
            Matrices {
                pfm,
                coverage,
                ppm,
                efficiency,
            }
        )
    }

    pub fn get_pfm(&self) -> ArrayView<u64, Ix2> {
        self.pfm.view()
    }

    pub fn get_coverage(&self) -> ArrayView<f64, Ix1> {
        self.coverage.view()
    }

    pub fn get_ppm(&self) -> ArrayView<f64, Ix2> {
        self.ppm.view()
    }

    pub fn get_efficiency(&self) -> ArrayView<f64, Ix1> {
        self.efficiency.view()
    }

    pub fn pfm_to_csv<P>(&self, path: P, filename: &str)
        where P: AsRef<Path>
    {
        let file = File::create(path.as_ref().join(Path::new(filename))).expect("could not open file");
        let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);

        writer.write_record(&["A", "C", "G", "T"]).unwrap();

        for col in self.pfm.columns() {
            writer.write_record(col.to_slice().unwrap().iter().map(|x| x.to_string()).collect::<Vec<_>>()).unwrap();
        }
    }
}