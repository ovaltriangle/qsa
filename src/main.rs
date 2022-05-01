use std::path::PathBuf;
use std::process::exit;

use structopt::StructOpt;
use inline_python::python;

use qsalib::prelude::*;

#[derive(Debug, StructOpt)]  // TODO: What **does** qsa exactly?
/// QuasiSpecies Analyser (QSA) is a bioinformatics tool which enables the analysis
/// of quasispecies viruses with ease.
struct QSAArgs {
    /// BAM files to be analysed or the directory containing them.
    ///
    /// The tool may give better and more complete results scaling with the amount
    /// of files utilised in the analysis.
    /// You might select multiple files and multiple directories by simply typing
    /// the name of those you want analysed.
    bams: Vec<PathBuf>,
    /// Starting range to be considered when selecting the reads to analyse.
    ///
    /// The starting range makes possible to discard all sequences starting before
    /// its value. The default value of 0 removes sequence selection by starting
    /// range altogether.
    #[structopt(short, long, default_value)]
    start: i32,
    /// Ending range to be considered when selecting the reads to analyse.
    ///
    /// The ending range makes possible to discard all sequences ending after
    /// its value. The default value of 0 removes sequence selection by ending
    /// range altogether.
    #[structopt(short, long, default_value)]
    end: i32,
    /// Minimum coverage value for the position to be considered valid.
    ///
    /// When reading a BAM file, the relative coverage of the reads over the
    /// reference sequence is calculated. When the threshold is set, the tool
    /// will discard the region of all read sequences where this threshold is
    /// not met. Use a value of 0 to disable this function.
    #[structopt(short, long, default_value = "0.65")]
    threshold: f64,
    /// Disables checks.
    ///
    /// As of now, the only check the program performs is that all BAMs have
    /// the same reference sequence. This check is performed to ensure there
    /// is no mistake in BAM files when starting the analysis process.
    /// If the BAM file has no header, this check is skipped.
    #[structopt(short, long)]
    no_checks: bool,
    /// Selects the output directory for the generated files.
    ///
    /// All output of the program will be put in a specific folder you can
    /// later check to see the result of all operations. If the specified
    /// folder does not exist, it will be created.
    #[structopt(short, long, default_value = "qsaout")]
    out_dir: PathBuf,
}

impl QSAArgs {
    fn validate(&self) -> () {
        if self.bams.is_empty() {
            eprintln!("You need to specify at least one or two BAM files.");
            exit(3)
        }
    }

    fn into_bamdata(self) -> Result<BamData> {
        let mut bams: Vec<PathBuf> = Vec::default();
        let mut dirs: Vec<PathBuf> = Vec::default();
        for path in self.bams {
            if path.is_dir() {
                dirs.push(path);
            } else if path.is_file() {
                bams.push(path);
            }
        }

        BamDataBuilder::default()
            .add_bams(bams)?
            .add_dirs(dirs)?
            .in_range((self.start, self.end))
            .with_threshold(self.threshold)
            .with_checks(!self.no_checks)
            .build()
    }
}

fn efficiency2graph(path: String, bam: &Bam) {
    let efficiency = bam.matrices.get_efficiency().to_vec();
    let filename = path + "/" + bam.name.as_str() + "-efficiency.png";

    python! {
        import matplotlib.pyplot as plt

        plt.figure(figsize=[10, 5])

        ax = plt.subplot(111)

        ax.plot('efficiency, "r", linewidth=1)

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        ax.set_axisbelow(True)
        ax.yaxis.grid(color="gray", linestyle="dashed")
        //ax.xaxis.grid(color="gray")

        plt.xlabel("position")
        plt.ylabel("efficiency")

        //plt.show()
        plt.margins(x=0.0075)
        plt.tight_layout(pad=0.15)
        plt.savefig('filename, transparent=True, bbox_inches="tight")
        plt.clf()
    }
}

fn alphadiv2graph(path: String, bamdata: &BamData) {
    let alpha = bamdata.alpha_diversity().to_vec();
    let labels = bamdata.get_names();
    let filename = path + "/alpha-diversity.png";

    python! {
        import matplotlib.pyplot as plt

        COLORS = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080", "#000000"]

        plt.figure(figsize=[8, 5])

        ax = plt.subplot(111)

        ax.bar('labels, 'alpha, color=COLORS)

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        ax.set_axisbelow(True)
        ax.yaxis.grid(color="gray", linestyle="dashed")

        //plt.xlabel("sample")
        plt.ylabel("α-diversity")

        //plt.show()
        plt.margins(x=0.015)
        plt.tight_layout(pad=0.15)
        plt.savefig('filename, transparent=True, bbox_inches="tight")
        plt.clf()
    }
}

fn betadiv2graph(path: String, bamdata: &BamData) {
    let beta = bamdata.beta_diversity();
    let mut betav: Vec<f64> = Vec::new();
    for i in 0..beta.nrows() {
        for j in (i+1)..beta.nrows() {
            betav.push(beta[[i, j]]);
        }
    }

    let labels = bamdata.get_names();
    let filename = path + "/beta-diversity.png";

    python! {
        import networkx as nx
        from networkx.drawing.nx_agraph import graphviz_layout
        import matplotlib.pyplot as plt
        from itertools import combinations as cmb

        COLORS = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080", "#000000"]

        plt.figure(figsize=[6.4, 4.8])

        labels_comb = list(cmb('labels, 2))
        edges_len = {labels_comb[i]: {"len": 'betav[i]} for i in range(len('betav))}

        G = nx.Graph()

        G.add_nodes_from('labels)
        G.add_edges_from(labels_comb)

        nx.set_edge_attributes(G, edges_len)

        ax = plt.subplot(111)

        ax.set_title("β-diversity")

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

        nx.draw_networkx_nodes(G, pos=graphviz_layout(G, prog="neato"), node_color=COLORS[:len('labels)])
        nx.draw_networkx_labels(G, pos=graphviz_layout(G, prog="neato"))
        //nx.draw_networkx(G, pos=graphviz_layout(G, prog="neato"), node_color=COLORS[:len('labels)])

        plt.tight_layout(pad=0.15)
        plt.savefig('filename, transparent=True, bbox_inches="tight")
        plt.clf()
    }
}

fn main() {
    let args: QSAArgs = QSAArgs::from_args();

    args.validate();

    let out_dir = args.out_dir.to_str().unwrap().to_owned();
    std::fs::create_dir_all(out_dir.clone()).expect("could not create output directory");

    let bam_data = args.into_bamdata();

    match bam_data {
        Ok(data) => {
            println!("All is OK, data built successfully");

            for bam in &data {
                efficiency2graph(out_dir.clone(), bam);
                bam.pfm_to_csv(out_dir.clone(), (bam.name.clone() + ".csv").as_str());
            }

            alphadiv2graph(out_dir.clone(), &data);
            
            betadiv2graph(out_dir.clone(), &data);
        },
        Err(why) => {
            eprintln!("{}", why);
            exit(4)
        }
    }
}