```markdown
# DynPocket

Dynamic pocket analysis pipeline software – an integrated toolkit for flexible protein dynamic active-site feature extraction and analysis, specifically designed for enzymes with a single major active pocket that exhibits significant conformational changes during catalysis.

Traditional static crystal-structure–based approaches often fail to capture the dynamic features of enzymes in their catalytic cycles. DynPocket bridges this gap by combining deep‑learning conformational generation models (BioEmu) or long‑timescale MD trajectories with classical pocket detection algorithms (fpocket). It batch‑generates protein conformational ensembles, tracks the dynamic trajectory of the dominant pocket, and extracts a core pocket ensemble (Top 30%). The pipeline outputs high‑frequency lining residues and a comprehensive set of physicochemical properties, making it an ideal auxiliary tool for enzymatic reaction mechanism studies, dynamic feature site identification, targeted molecular design, and classification model construction.

---

## Key Features

- **Multi‑mode intelligent pocket tracking** – Supports either "unsupervised global clustering (Mode A)" based on centroid distances from a reference conformation, or "targeted residue‑guided (Mode B)" by specifying catalytic residues, to precisely pinpoint the Top 30% core active pocket.

- **Conformational reduction and redundancy removal** – Employs a greedy screening algorithm with pairwise Cα RMSD > 5.0 Å to clean large‑scale BioEmu‑generated trajectories, extracting representative non‑redundant conformational frames.

- **Mandatory consistency safety check** – Rigorously maps FASTA 1‑based indices to PDB genuine residue numbers; terminates immediately upon any mismatch, preventing catastrophic coordinate misalignment downstream.

- **Multi‑process accelerated detection** – Leverages Python `multiprocessing` process pools to run fpocket calculations in parallel across the reduced conformational frames.

- **Automated feature extraction & visualisation** – Automatically aggregates dynamic physicochemical properties (e.g., pocket volume, polar surface area), generates line plots showing property variations across conformations, and produces one‑click PyMOL rendering scripts to highlight high‑frequency lining residues.

- **Enterprise‑grade engineering design** – Natively supports multi‑FASTA batch parallelisation with isolated per‑protein directories, built‑in checkpoint/resume (`--resume`) for interrupted jobs, and full persistent logging for complete traceability.

---

## Core Workflow

1. **Input & validation**  
   Parses sequences and structures; by default uses BioEmu frame 0 (or a user‑provided PDB) as the reference for coordinate alignment. A hard limit of 1500 AA is enforced to prevent GPU memory overflow (can be overridden with `--force`).

2. **Conformation generation & redundancy reduction**  
   Generates an ensemble (default 1000 frames) via BioEmu (time scales as O(N²) with sequence length), followed by initial RMSD filtering and iterative greedy selection. (Alternatively, MD `.xtc` trajectories plus a topology PDB can be used.)

3. **Pocket calculation & filtering**  
   Runs fpocket in parallel; maps predicted pocket‑lining residues back to the reference conformation, then selects the Top 30% core pockets by Euclidean distance to the global centroid or to the target residues.

4. **Feature extraction & output**  
   Produces structured CSV tables summarising mean/variance of pocket properties across frames, and outputs high‑frequency core lining residues sorted by occurrence frequency.

---

## Command‑Line Interface (CLI)

Invoke the main pipeline via the entry script `main.py`:

### Required Arguments

| Argument | Description |
| --- | --- |
| `-i`, `--fasta` | Input protein sequence file path (supports single or multi‑FASTA). |

### Optional & Environment Control Arguments

| Argument | Description |
| --- | --- |
| `-p`, `--pdb` | User‑provided reference PDB file for coordinate mapping (single‑sequence mode). |
| `--pdb_dir` | Directory containing PDB files for batch mode; filenames must be `{FASTA_header}.pdb`. |
| `-o`, `--out_dir` | Main output directory; defaults to `./dynpocket_out`. |
| `--frame_num` | Number of conformations to generate with BioEmu; defaults to `1000`. |
| `--force` | Force execution, bypassing the 1500 AA length safety limit. |
| `--check-env` | Boolean flag; performs a self‑check of runtime dependencies (BioEmu, fpocket) only. |
| `--resume` | Enables global checkpoint/resume, reading `progress.json` to restore an interrupted job. |

### Pocket Filtering & Clustering Arguments

| Argument | Description |
| --- | --- |
| `--target_residues` | **[Advanced Mode B]** Specify residues near or within the target pocket (e.g., `"15,48,102"`) to guide focused tracking. |
| `--centroid_top_ratio` | Ratio of top frames to retain based on centroid distance; defaults to `0.3` (i.e., keep the top 30% core pockets). |
| `--preset` | fpocket search parameters; defaults to `-i 80 -m 5.0 -M 10.0 -D 3.0`. |
| `--inner_shell_residues_num` | Number of top high‑frequency lining residues to output; defaults to `20` (set to `-1` to output all). |

### Standalone Visualisation Sub‑command

To generate highlighting scripts for a specific conformation frame and a specific pocket, use the standalone `vis` sub‑command:

```bash
dynpocket-analyzer vis <protein_dir> --pocket <number>
```

This command parses the results from `pockets.csv` (e.g., `--pocket 10_1`) and generates a dedicated PyMOL visualisation script.

---

## Output Directory Structure

Upon successful completion, the pipeline creates an isolated directory per protein (based on the FASTA ID) under the specified `--out_dir`. The standard structure is as follows:

```text
--out_dir/        
 └── protein_ID_1/
     ├── bioemu_out/                  # Raw BioEmu output files
     ├── conformations_raw/           # Parsed BioEmu initial valid frames
     ├── conformations_representative/# Representative frames after RMSD reduction
     ├── pocket_data/                 # fpocket raw outputs for each representative frame
     ├── inner_residues.csv           # Lining residues for all detected pockets
     ├── pockets.csv                  # Top 30% pocket averages and per‑frame physicochemical property scores
     ├── visualization_assets/        # Visualisation auxiliary files
     │   ├── properties_variation_plot.png  # Dynamic variation plot of pocket volume, polarity, etc.
     │   ├── inner_residues_list.txt        # Plain‑text list of top high‑frequency lining residues (default first 20)
     │   └── highlight_pocket.pml           # One‑click PyMOL script for automated display
     ├── progress.json                # Checkpoint state for resume functionality
     └── run.log                      # Persistent log with validation status, estimated runtime, and processing details
```
```markdown
# DynPocket

Dynamic pocket analysis pipeline software – an integrated toolkit for flexible protein dynamic active-site feature extraction and analysis, specifically designed for enzymes with a single major active pocket that exhibits significant conformational changes during catalysis.

Traditional static crystal-structure–based approaches often fail to capture the dynamic features of enzymes in their catalytic cycles. DynPocket bridges this gap by combining deep‑learning conformational generation models (BioEmu) or long‑timescale MD trajectories with classical pocket detection algorithms (fpocket). It batch‑generates protein conformational ensembles, tracks the dynamic trajectory of the dominant pocket, and extracts a core pocket ensemble (Top 30%). The pipeline outputs high‑frequency lining residues and a comprehensive set of physicochemical properties, making it an ideal auxiliary tool for enzymatic reaction mechanism studies, dynamic feature site identification, targeted molecular design, and classification model construction.

---

## Key Features

- **Multi‑mode intelligent pocket tracking** – Supports either "unsupervised global clustering (Mode A)" based on centroid distances from a reference conformation, or "targeted residue‑guided (Mode B)" by specifying catalytic residues, to precisely pinpoint the Top 30% core active pocket.

- **Conformational reduction and redundancy removal** – Employs a greedy screening algorithm with pairwise Cα RMSD > 5.0 Å to clean large‑scale BioEmu‑generated trajectories, extracting representative non‑redundant conformational frames.

- **Mandatory consistency safety check** – Rigorously maps FASTA 1‑based indices to PDB genuine residue numbers; terminates immediately upon any mismatch, preventing catastrophic coordinate misalignment downstream.

- **Multi‑process accelerated detection** – Leverages Python `multiprocessing` process pools to run fpocket calculations in parallel across the reduced conformational frames.

- **Automated feature extraction & visualisation** – Automatically aggregates dynamic physicochemical properties (e.g., pocket volume, polar surface area), generates line plots showing property variations across conformations, and produces one‑click PyMOL rendering scripts to highlight high‑frequency lining residues.

- **Enterprise‑grade engineering design** – Natively supports multi‑FASTA batch parallelisation with isolated per‑protein directories, built‑in checkpoint/resume (`--resume`) for interrupted jobs, and full persistent logging for complete traceability.

---

## Core Workflow

1. **Input & validation**  
   Parses sequences and structures; by default uses BioEmu frame 0 (or a user‑provided PDB) as the reference for coordinate alignment. A hard limit of 1500 AA is enforced to prevent GPU memory overflow (can be overridden with `--force`).

2. **Conformation generation & redundancy reduction**  
   Generates an ensemble (default 1000 frames) via BioEmu (time scales as O(N²) with sequence length), followed by initial RMSD filtering and iterative greedy selection. (Alternatively, MD `.xtc` trajectories plus a topology PDB can be used.)

3. **Pocket calculation & filtering**  
   Runs fpocket in parallel; maps predicted pocket‑lining residues back to the reference conformation, then selects the Top 30% core pockets by Euclidean distance to the global centroid or to the target residues.

4. **Feature extraction & output**  
   Produces structured CSV tables summarising mean/variance of pocket properties across frames, and outputs high‑frequency core lining residues sorted by occurrence frequency.

---

## Command‑Line Interface (CLI)

Invoke the main pipeline via the entry script `main.py`:

### Required Arguments

| Argument | Description |
| --- | --- |
| `-i`, `--fasta` | Input protein sequence file path (supports single or multi‑FASTA). |

### Optional & Environment Control Arguments

| Argument | Description |
| --- | --- |
| `-p`, `--pdb` | User‑provided reference PDB file for coordinate mapping (single‑sequence mode). |
| `--pdb_dir` | Directory containing PDB files for batch mode; filenames must be `{FASTA_header}.pdb`. |
| `-o`, `--out_dir` | Main output directory; defaults to `./dynpocket_out`. |
| `--frame_num` | Number of conformations to generate with BioEmu; defaults to `1000`. |
| `--force` | Force execution, bypassing the 1500 AA length safety limit. |
| `--check-env` | Boolean flag; performs a self‑check of runtime dependencies (BioEmu, fpocket) only. |
| `--resume` | Enables global checkpoint/resume, reading `progress.json` to restore an interrupted job. |

### Pocket Filtering & Clustering Arguments

| Argument | Description |
| --- | --- |
| `--target_residues` | **[Advanced Mode B]** Specify residues near or within the target pocket (e.g., `"15,48,102"`) to guide focused tracking. |
| `--centroid_top_ratio` | Ratio of top frames to retain based on centroid distance; defaults to `0.3` (i.e., keep the top 30% core pockets). |
| `--preset` | fpocket search parameters; defaults to `-i 80 -m 5.0 -M 10.0 -D 3.0`. |
| `--inner_shell_residues_num` | Number of top high‑frequency lining residues to output; defaults to `20` (set to `-1` to output all). |

### Standalone Visualisation Sub‑command

To generate highlighting scripts for a specific conformation frame and a specific pocket, use the standalone `vis` sub‑command:

```bash
dynpocket-analyzer vis <protein_dir> --pocket <number>
```

This command parses the results from `pockets.csv` (e.g., `--pocket 10_1`) and generates a dedicated PyMOL visualisation script.

---

## Output Directory Structure

Upon successful completion, the pipeline creates an isolated directory per protein (based on the FASTA ID) under the specified `--out_dir`. The standard structure is as follows:

```text
--out_dir/        
 └── protein_ID_1/
     ├── bioemu_out/                  # Raw BioEmu output files
     ├── conformations_raw/           # Parsed BioEmu initial valid frames
     ├── conformations_representative/# Representative frames after RMSD reduction
     ├── pocket_data/                 # fpocket raw outputs for each representative frame
     ├── inner_residues.csv           # Lining residues for all detected pockets
     ├── pockets.csv                  # Top 30% pocket averages and per‑frame physicochemical property scores
     ├── visualization_assets/        # Visualisation auxiliary files
     │   ├── properties_variation_plot.png  # Dynamic variation plot of pocket volume, polarity, etc.
     │   ├── inner_residues_list.txt        # Plain‑text list of top high‑frequency lining residues (default first 20)
     │   └── highlight_pocket.pml           # One‑click PyMOL script for automated display
     ├── progress.json                # Checkpoint state for resume functionality
     └── run.log                      # Persistent log with validation status, estimated runtime, and processing details
```
