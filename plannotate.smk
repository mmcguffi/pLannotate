# pLannotate Snakemake Pipeline
# Modular plasmid annotation pipeline using core pLannotate functions

import os
import glob
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Global variables
INPUT_DIR = config.get("input_dir", "input")
OUTPUT_DIR = config.get("output_dir", "output")
TEMP_DIR = config.get("temp_dir", "temp")
LOG_DIR = config.get("log_dir", "logs")
YAML_FILE = config.get("yaml_file", "default")

# Get sample names from input files
def get_sample_names():
    """Get sample names from input files."""
    input_pattern = f"{INPUT_DIR}/*.{{fa,fasta,fas,fna,gbk,gb,gbf,gbff}}"
    input_files = glob.glob(input_pattern)
    return [Path(f).stem for f in input_files]

# Helper function to get file extension
def get_file_ext(sample):
    """Get file extension for a given sample."""
    for ext in ['fa', 'fasta', 'fas', 'fna', 'gbk', 'gb', 'gbf', 'gbff']:
        if os.path.exists(f"{INPUT_DIR}/{sample}.{ext}"):
            return ext
    return 'fa'  # default

SAMPLES = get_sample_names()

# Rule definitions
rule all:
    """Main rule that runs the complete pipeline."""
    input:
        expand(f"{OUTPUT_DIR}/{{sample}}_annotated.gbk", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}_annotations.csv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}_summary.txt", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}_plot.html", sample=SAMPLES) if config.get("generate_plots", True) else []
    message:
        "pLannotate modular pipeline completed successfully"

rule setup_directories:
    """Create necessary directories and check database setup."""
    output:
        directory(f"{OUTPUT_DIR}"),
        directory(f"{TEMP_DIR}"),
        directory(f"{LOG_DIR}"),
        touch(f"{TEMP_DIR}/databases_checked.txt")
    log:
        f"{LOG_DIR}/setup_directories.log"
    shell:
        """
        mkdir -p {output[0]} {output[1]} {output[2]} && \\
        python3 -c "
import sys
sys.path.insert(0, '.')
from plannotate import resources
if not resources.databases_exist():
    print('ERROR: Databases not found. Run plannotate setupdb first.', file=sys.stderr)
    sys.exit(1)
else:
    print('Databases validated successfully')
        " \\
        >& {log}
        """

rule validate_input:
    """Validate input sequence files and prepare for annotation."""
    input:
        seq_file = lambda wildcards: f"{INPUT_DIR}/{wildcards.sample}.{get_file_ext(wildcards.sample)}"
    output:
        validated = f"{TEMP_DIR}/{{sample}}_validated.fasta"
    params:
        max_length = config.get("max_sequence_length", 50000)
    log:
        f"{LOG_DIR}/validate_input_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
from plannotate.validation import validate_and_write_fasta
try:
    validate_and_write_fasta('{input.seq_file}', '{output.validated}', {params.max_length})
    print('Validation successful: {wildcards.sample}')
except Exception as e:
    print(f'Error validating {input.seq_file}: {{e}}', file=sys.stderr)
    sys.exit(1)
        " \\
        >& {log}
        """

rule infernal:
    """Run Infernal search against Rfam database for ncRNA annotations."""
    input:
        seq = f"{TEMP_DIR}/{{sample}}_validated.fasta",
        db_check = f"{TEMP_DIR}/databases_checked.txt"
    output:
        hits = f"{TEMP_DIR}/{{sample}}_infernal_hits.csv"
    params:
        linear = config.get("linear", False),
        yaml_file = lambda wildcards: YAML_FILE if YAML_FILE != "default" else ""
    log:
        f"{LOG_DIR}/infernal_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
import pandas as pd
from Bio import SeqIO
from plannotate import annotate, resources

# Read sequence
with open('{input.seq}', 'r') as f:
    records = list(SeqIO.parse(f, 'fasta'))
    seq = str(records[0].seq)

# Get database config
yaml_path = resources.get_yaml_path() if '{params.yaml_file}' == '' else '{params.yaml_file}'
databases = resources.get_yaml(yaml_path)

# Run Infernal search specifically
if 'Rfam' in databases:
    db_config = databases['Rfam']
    if not {params.linear}:
        query_seq = seq + seq  # double for circular
    else:
        query_seq = seq
        
    hits = annotate.infernal(query_seq, db_config)
    hits['db'] = 'Rfam'
    
    if not hits.empty:
        hits = annotate.calculate(hits, is_linear={params.linear})
        
    hits.to_csv('{output.hits}', index=False)
    print(f'Infernal search complete: {{len(hits)}} hits found')
else:
    # Create empty CSV if no Rfam database
    pd.DataFrame().to_csv('{output.hits}', index=False)
    print('No Rfam database configured')
        " \\
        >& {log}
        """

rule snapgene:
    """Run BLASTn search against SnapGene database."""
    input:
        seq = f"{TEMP_DIR}/{{sample}}_validated.fasta",
        db_check = f"{TEMP_DIR}/databases_checked.txt"
    output:
        hits = f"{TEMP_DIR}/{{sample}}_snapgene_hits.csv"
    params:
        linear = config.get("linear", False),
        yaml_file = lambda wildcards: YAML_FILE if YAML_FILE != "default" else ""
    log:
        f"{LOG_DIR}/snapgene_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
import pandas as pd
from Bio import SeqIO
from plannotate import annotate, resources

# Read sequence
with open('{input.seq}', 'r') as f:
    records = list(SeqIO.parse(f, 'fasta'))
    seq = str(records[0].seq)

# Get database config
yaml_path = resources.get_yaml_path() if '{params.yaml_file}' == '' else '{params.yaml_file}'
databases = resources.get_yaml(yaml_path)

# Run SnapGene search specifically
if 'snapgene' in databases:
    db_config = databases['snapgene']
    if not {params.linear}:
        query_seq = seq + seq  # double for circular
    else:
        query_seq = seq
        
    hits = annotate.blast(query_seq, db_config)
    hits['db'] = 'snapgene'
    
    if not hits.empty:
        hits = annotate.calculate(hits, is_linear={params.linear})
        
    hits.to_csv('{output.hits}', index=False)
    print(f'SnapGene search complete: {{len(hits)}} hits found')
else:
    # Create empty CSV if no snapgene database
    pd.DataFrame().to_csv('{output.hits}', index=False)
    print('No SnapGene database configured')
        " \\
        >& {log}
        """

rule swissprot:
    """Run DIAMOND search against SwissProt database."""
    input:
        seq = f"{TEMP_DIR}/{{sample}}_validated.fasta",
        db_check = f"{TEMP_DIR}/databases_checked.txt"
    output:
        hits = f"{TEMP_DIR}/{{sample}}_swissprot_hits.csv"
    params:
        linear = config.get("linear", False),
        yaml_file = lambda wildcards: YAML_FILE if YAML_FILE != "default" else ""
    log:
        f"{LOG_DIR}/swissprot_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
import pandas as pd
from Bio import SeqIO
from plannotate import annotate, resources

# Read sequence
with open('{input.seq}', 'r') as f:
    records = list(SeqIO.parse(f, 'fasta'))
    seq = str(records[0].seq)

# Get database config
yaml_path = resources.get_yaml_path() if '{params.yaml_file}' == '' else '{params.yaml_file}'
databases = resources.get_yaml(yaml_path)

# Run SwissProt search specifically
if 'swissprot' in databases:
    db_config = databases['swissprot']
    if not {params.linear}:
        query_seq = seq + seq  # double for circular
    else:
        query_seq = seq
        
    hits = annotate.diamond(query_seq, db_config)
    hits['db'] = 'swissprot'
    
    if not hits.empty:
        hits = annotate.calculate(hits, is_linear={params.linear})
        
    hits.to_csv('{output.hits}', index=False)
    print(f'SwissProt search complete: {{len(hits)}} hits found')
else:
    # Create empty CSV if no swissprot database
    pd.DataFrame().to_csv('{output.hits}', index=False)
    print('No SwissProt database configured')
        " \\
        >& {log}
        """

rule fpbase:
    """Run DIAMOND search against FPbase database."""
    input:
        seq = f"{TEMP_DIR}/{{sample}}_validated.fasta",
        db_check = f"{TEMP_DIR}/databases_checked.txt"
    output:
        hits = f"{TEMP_DIR}/{{sample}}_fpbase_hits.csv"
    params:
        linear = config.get("linear", False),
        yaml_file = lambda wildcards: YAML_FILE if YAML_FILE != "default" else ""
    log:
        f"{LOG_DIR}/fpbase_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
import pandas as pd
from Bio import SeqIO
from plannotate import annotate, resources

# Read sequence
with open('{input.seq}', 'r') as f:
    records = list(SeqIO.parse(f, 'fasta'))
    seq = str(records[0].seq)

# Get database config
yaml_path = resources.get_yaml_path() if '{params.yaml_file}' == '' else '{params.yaml_file}'
databases = resources.get_yaml(yaml_path)

# Run FPbase search specifically
if 'fpbase' in databases:
    db_config = databases['fpbase']
    if not {params.linear}:
        query_seq = seq + seq  # double for circular
    else:
        query_seq = seq
        
    hits = annotate.diamond(query_seq, db_config)
    hits['db'] = 'fpbase'
    
    if not hits.empty:
        hits = annotate.calculate(hits, is_linear={params.linear})
        
    hits.to_csv('{output.hits}', index=False)
    print(f'FPbase search complete: {{len(hits)}} hits found')
else:
    # Create empty CSV if no fpbase database
    pd.DataFrame().to_csv('{output.hits}', index=False)
    print('No FPbase database configured')
        " \\
        >& {log}
        """

rule clean_and_combine:
    """Combine hits from all databases and apply pLannotate filtering."""
    input:
        infernal = f"{TEMP_DIR}/{{sample}}_infernal_hits.csv",
        snapgene = f"{TEMP_DIR}/{{sample}}_snapgene_hits.csv", 
        swissprot = f"{TEMP_DIR}/{{sample}}_swissprot_hits.csv",
        fpbase = f"{TEMP_DIR}/{{sample}}_fpbase_hits.csv",
        seq = f"{TEMP_DIR}/{{sample}}_validated.fasta"
    output:
        combined = f"{TEMP_DIR}/{{sample}}_combined_hits.csv",
        filtered = f"{TEMP_DIR}/{{sample}}_filtered_hits.csv"
    params:
        linear = config.get("linear", False),
        detailed = config.get("detailed", False),
        yaml_file = lambda wildcards: YAML_FILE if YAML_FILE != "default" else ""
    log:
        f"{LOG_DIR}/clean_and_combine_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
import pandas as pd
from Bio import SeqIO
from plannotate import annotate, resources

# Read all hit files
hit_files = ['{input.infernal}', '{input.snapgene}', '{input.swissprot}', '{input.fpbase}']
all_hits = []

for hit_file in hit_files:
    try:
        hits = pd.read_csv(hit_file)
        if not hits.empty:
            all_hits.append(hits)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        print(f'No hits in {{hit_file}}')
        continue

# Combine all hits
if all_hits:
    combined_df = pd.concat(all_hits, ignore_index=True)
    
    # Add feature details for each database
    yaml_path = resources.get_yaml_path() if '{params.yaml_file}' == '' else '{params.yaml_file}'
    
    # Group by database and add details
    detailed_hits = []
    for db_name in combined_df['db'].unique():
        db_hits = combined_df[combined_df['db'] == db_name]
        if not db_hits.empty:
            feat_details = annotate.get_details(db_hits, yaml_path)
            db_hits = db_hits.merge(feat_details, on='sseqid', how='left', suffixes=('_x', None))
            db_hits = db_hits[db_hits.columns.drop(list(db_hits.filter(regex='_x')))]
            detailed_hits.append(db_hits)
    
    if detailed_hits:
        combined_df = pd.concat(detailed_hits, ignore_index=True)
    
    # Sort by score
    combined_df = combined_df.sort_values(
        by=['score', 'length', 'percmatch'], 
        ascending=[False, False, False]
    )
    
    combined_df.to_csv('{output.combined}', index=False)
    
    # Apply cleaning and filtering
    if {params.detailed}:
        combined_df['kind'] = combined_df['Type']
    else:
        combined_df['kind'] = 1
        
    filtered_df = annotate.clean(combined_df)
    filtered_df.to_csv('{output.filtered}', index=False)
    
    print(f'Combined {{len(combined_df)}} hits, filtered to {{len(filtered_df)}} hits')
else:
    # Create empty files if no hits
    empty_df = pd.DataFrame()
    empty_df.to_csv('{output.combined}', index=False)
    empty_df.to_csv('{output.filtered}', index=False)
    print('No hits found in any database')
        " \\
        >& {log}
        """

rule generate_outputs:
    """Generate final output files (GenBank, CSV, summary)."""
    input:
        filtered = f"{TEMP_DIR}/{{sample}}_filtered_hits.csv",
        seq = f"{TEMP_DIR}/{{sample}}_validated.fasta"
    output:
        gbk = f"{OUTPUT_DIR}/{{sample}}_annotated.gbk",
        csv = f"{OUTPUT_DIR}/{{sample}}_annotations.csv",
        summary = f"{OUTPUT_DIR}/{{sample}}_summary.txt"
    params:
        linear = config.get("linear", False),
        detailed = config.get("detailed", False),
        yaml_file = lambda wildcards: YAML_FILE if YAML_FILE != "default" else ""
    log:
        f"{LOG_DIR}/generate_outputs_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
import pandas as pd
from Bio import SeqIO
from plannotate.models import Construct
from plannotate import resources

# Read sequence
with open('{input.seq}', 'r') as f:
    records = list(SeqIO.parse(f, 'fasta'))
    seq_record = records[0]

# Read filtered hits
try:
    hits_df = pd.read_csv('{input.filtered}')
except (pd.errors.EmptyDataError, FileNotFoundError):
    hits_df = pd.DataFrame()

# Get yaml path
yaml_path = resources.get_yaml_path() if '{params.yaml_file}' == '' else '{params.yaml_file}'

# Create Construct object
construct = Construct(
    seq=seq_record.seq,
    linear={params.linear},
    detailed={params.detailed},
    db_options=yaml_path
)

# Generate outputs
gbk_content = construct.to_genbank()
with open('{output.gbk}', 'w') as f:
    f.write(gbk_content)

csv_df = construct.to_csv()
csv_df.to_csv('{output.csv}', index=False)

# Generate summary
with open('{output.summary}', 'w') as f:
    f.write(f'pLannotate Annotation Summary\\n')
    f.write(f'Sample: {wildcards.sample}\\n')
    f.write(f'Sequence length: {{len(seq_record.seq)}} bp\\n')
    f.write(f'Linear: {params.linear}\\n')
    f.write(f'Detailed: {params.detailed}\\n')
    f.write(f'Total annotations: {{len(csv_df)}}\\n')
    
    if not csv_df.empty:
        type_counts = csv_df['Type'].value_counts()
        f.write('\\nAnnotation types:\\n')
        for ann_type, count in type_counts.items():
            f.write(f'  {{ann_type}}: {{count}}\\n')

print(f'Generated outputs for {wildcards.sample}')
        " \\
        >& {log}
        """

rule plot:
    """Generate interactive HTML plot (optional)."""
    input:
        filtered = f"{TEMP_DIR}/{{sample}}_filtered_hits.csv",
        seq = f"{TEMP_DIR}/{{sample}}_validated.fasta"
    output:
        html = f"{OUTPUT_DIR}/{{sample}}_plot.html"
    params:
        linear = config.get("linear", False),
        detailed = config.get("detailed", False),
        yaml_file = lambda wildcards: YAML_FILE if YAML_FILE != "default" else ""
    log:
        f"{LOG_DIR}/plot_{{sample}}.log"
    shell:
        """
        python3 -c "
import sys
sys.path.insert(0, '.')
from Bio import SeqIO
from plannotate.models import Construct
from plannotate import resources

# Read sequence
with open('{input.seq}', 'r') as f:
    records = list(SeqIO.parse(f, 'fasta'))
    seq_record = records[0]

# Get yaml path
yaml_path = resources.get_yaml_path() if '{params.yaml_file}' == '' else '{params.yaml_file}'

# Create Construct object
construct = Construct(
    seq=seq_record.seq,
    linear={params.linear},
    detailed={params.detailed},
    db_options=yaml_path
)

# Generate HTML plot
html_content = construct.to_html(htmlfull=True)
with open('{output.html}', 'w') as f:
    f.write(html_content)

print(f'Generated plot for {wildcards.sample}')
        " \\
        >& {log}
        """

# Cleanup rule
rule cleanup:
    """Remove temporary files after pipeline completion."""
    input:
        expand(f"{OUTPUT_DIR}/{{sample}}_annotated.gbk", sample=SAMPLES)
    output:
        touch(f"{OUTPUT_DIR}/cleanup_complete.txt")
    params:
        keep_temp = config.get("keep_temp", False)
    log:
        f"{LOG_DIR}/cleanup.log"
    shell:
        """
        if [ "{params.keep_temp}" != "True" ]; then \\
            echo "Removing temporary files..." && \\
            rm -rf {TEMP_DIR} && \\
            echo "Cleanup complete" \\
        else \\
            echo "Keeping temporary files as requested" \\
        fi \\
        >& {log}
        """