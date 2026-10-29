# 🍄 Fungalcaller: End-to-End Nextflow Pipeline for Fungal Comparative Genomics & Phylogenomics

An automated, containerized **Nextflow (DSL2)** pipeline designed for comprehensive processing, structural/functional annotation, secretome/effector profiling, and high-resolution phylogenomic reconstruction of fungal genomes (optimized and tested on *Trichoderma, Fusarium and Clonostachys*).

---

## 🚀 Key Features

* **Flexible Input:** Seamlessly handles both raw short-read sequencing data (`FASTQ`) via automated assembly or pre-assembled genome sequences (`FASTA`).
* **Quality Control & Assembly:** FastQC, FastP trimming, SPAdes assembly, and **Tiara**-based contig cleaning/filtering.
* **Genome Completeness & Statistics:** Automated **BUSCO** evaluation and **Jellyfish** k-mer profiling with summary aggregation.
* **Structural & Functional Annotation:** RepeatMasker (repeat masking), Augustus (ab initio gene prediction), and AntiSMASH (biosynthetic gene clusters).
* **Functional Profiling:** High-throughput functional annotation via **EggNOG-mapper** and **InterProScan**.
* **Secretome & Effector Discovery:** Dedicated prediction pipeline using **DeepSig** (signal peptides) and **EffectorP** (fungal effectors).
* **Comparative & Phylogenomic Pipeline:** Orthologue group inference via **OrthoFinder**, multiple sequence alignment with **MAFFT**, gap trimming with **trimAl**, concatenated supermatrix generation, and Maximum Likelihood tree building with **RAxML-NG** (1000 bootstrap replicates).
* **Full Reproducibility:** Containerized execution via **Docker** (`mordziarz/fungi_pipeline:latest`) managed through Nextflow process resources.

---

## ⚙️ Input Data (samples_table.csv)

The pipeline expects a CSV file specifying the sample names, input type (`fasta` or `fastq`), and paths to reads or assembly files.

### Structure & Column Descriptions:
- **`sample`**: Unique sample identifier / strain name.
- **`type`**: Input data type (`fasta` or `fastq`).
- **`reads_1`**: Path to forward FASTQ file (leave empty if type is `fasta`).
- **`reads_2`**: Path to reverse FASTQ file (leave empty if type is `fasta`).
- **`fasta`**: Path to pre-assembled FASTA genome (leave empty if type is `fastq`).

### Example `samples_table.csv`:

```csv
sample,type,reads_1,reads_2,fasta
sample1,fasta,,,/path/to/genomes/strain_A.fna
sample2,fasta,,,/path/to/genomes/strain_B.fna
Trichoderma_paratroviride,fastq,/path/to/reads/strain_C_R1.fastq.gz,/path/to/reads/strain_C_R2.fastq.gz,
```

## 🛠️ Requirements & Dependencies

```text
- Nextflow (version 22.0 or higher)
- Docker (with access to container: mordziarz/fungi_pipeline:latest)
- Recommended Hardware:
  - CPU: 8+ cores
  - RAM: 16+ GB (recommended for annotation tools and OrthoFinder)
  - Storage: ~150 GB free space (required for local databases: EggNOG, antiSMASH, InterProScan)
```

## 📥 Installation & Usage

### Step 1: Download required databases
Before running the pipeline, ensure your biological databases are downloaded and configured locally (e.g., via the provided `download_databases.sh` script):
```bash
./download_databases.sh
```

### Step 2: Run the pipeline
You can run **Fungalcaller** directly from GitHub using Nextflow:

```bash
nextflow run Mordziarz/Fungalcaller \
    -r main \
    --samples_table samples_table.csv \
    --outgroup "sample4" \
    -profile docker \
    -resume
```

## 📂 Output Structure

Results are automatically saved to your specified output directory (`results/` by default):

```text
results/
├── antismash/                  # Biosynthetic gene cluster predictions (GenBank & HTML reports)
├── busco/                      # BUSCO completeness evaluations & aggregated summary
├── concatenated_alignment/     # Trimmed supermatrix alignments & partition files
├── deepsig/                    # Predicted signal peptides & secretory proteins
├── effectorp/                  # Candidate fungal effector predictions
├── eggnog/                     # High-throughput EggNOG functional annotations
├── fastp/                      # Adapter trimming & sequence quality filtering
├── fastqc/                     # Raw & trimmed read quality control reports
├── interproscan/               # Protein domain & functional family annotations
├── orthofinder/                # Orthogroup inferences & single-copy orthologue sequences
├── raxml/                      # Maximum Likelihood species tree (.support & .bestTree)
├── repeatmasker/               # Repeat masking reports & masked FASTA assemblies
├── spades/                     # De novo genome assemblies (for FASTQ inputs)
├── summary/                    # Aggregated quality & genome assembly statistics
└── tiara/                      # Cleaned & non-organellar/contaminant filtered contigs
```

## 🔮 Future Perspectives

* Integration of Oxford Nanopore long-read sequencing data for hybrid and automated de novo nanopore genome assemblies.
* Expansion of downstream comparative genomics modules.

---

## 📖 Citation & Paper in Preparation

If you use **Fungalcaller** in your research, a manuscript describing the pipeline is currently in preparation. 

---

## 📬 Contact & Support

For any questions, bug reports, or feature requests, please feel free to reach out:

* **Mateusz Maździarz**
* Email: `mateusz.mazdziarz@uwm.edu.pl` or `mazdziarzm@gmail.com`
* GitHub: [@Mordziarz](https://github.com/Mordziarz)