# Streamlit UI
import streamlit as st
from streamlit_extras.colored_header import colored_header
from streamlit_extras.stylable_container import stylable_container
from streamlit_file_browser import st_file_browser
import datetime
import glob
import gzip
import io
import os
import re
import shutil
import sqlite3
import subprocess
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from pathlib import Path
import pandas as pd
from Bio import Entrez, SeqIO
from ete3 import NCBITaxa
from tqdm import tqdm
import os
import sys
import gzip
import shutil
import subprocess
from pathlib import Path
import tempfile
import webbrowser
import platform
import numpy as np
import pyarrow.parquet as pq
import pyarrow as pa
from stqdm import stqdm
import math
from functools import lru_cache

# suppress ncbi warning
import warnings
warnings.filterwarnings("ignore", message=r"taxid \d+ was translated into \d+")

### Set page config first
st.set_page_config(page_title="APSCALE-blast Database Creator", page_icon="🧬", layout="wide")
path_to_apscale_blast = Path(__file__).resolve().parent

## accession to taxid conversion (or update the existing file)

# Initialize ete3 NCBI database
ncbi = NCBITaxa()

# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
# db_path = '/Volumes/Coruscant/APSCALE_raw_databases/accession2taxid/db/nucl_gb.accession2taxid.db'
# gz_path = '/Volumes/Coruscant/APSCALE_raw_databases/accession2taxid/nucl_gb.accession2taxid.gz'

def create_acc2taxid_db(gz_path, db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS acc2taxid")
    c.execute("CREATE TABLE acc2taxid (accession TEXT PRIMARY KEY, taxid INTEGER)")

    with gzip.open(gz_path, "rt") as f:
        next(f)  # skip header
        batch = []
        for i, line in enumerate(f):
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                accession = fields[0]
                taxid = int(fields[2])
                batch.append((accession, taxid))
            if i % 100000 == 0:
                c.executemany("INSERT OR IGNORE INTO acc2taxid VALUES (?, ?)", batch)
                conn.commit()
                batch = []
        if batch:
            c.executemany("INSERT OR IGNORE INTO acc2taxid VALUES (?, ?)", batch)
            conn.commit()

    conn.close()

def get_taxonomy_from_accession_ete3(accession, db_path):

    def get_taxid_from_sqlite(accession, db_path):
        conn = sqlite3.connect(db_path)
        c = conn.cursor()
        c.execute("SELECT taxid FROM acc2taxid WHERE accession=?", (accession,))
        result = c.fetchone()
        conn.close()
        return result[0] if result else None

    try:
        # Step 1: Try local lookup
        taxid = get_taxid_from_sqlite(accession, db_path)

        # Step 2: If not found, fallback to Entrez
        if not taxid:
            handle = Entrez.esummary(db="nuccore", id=accession)
            summary = Entrez.read(handle)
            handle.close()
            taxid = int(summary[0]["TaxId"])

        # Step 3: Use ete3 to get lineage
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        # Map desired ranks
        wanted_ranks = {
            'superkingdom': 'unclassified',
            'phylum': 'unclassified',
            'class': 'unclassified',
            'order': 'unclassified',
            'family': 'unclassified',
            'genus': 'unclassified',
            'species': 'unclassified'
        }

        for tid in lineage:
            rank = ranks.get(tid)
            name = names.get(tid)
            if rank in wanted_ranks:
                wanted_ranks[rank] = name

        wanted_ranks['Accession'] = accession
        return wanted_ranks

    except Exception as e:
        return {'Accession': accession, 'Error': str(e)}

## NCBI ACCESSION CONVERTER

def accession_to_taxonomy(accession):
    # Fetch the record summary from NCBI
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Try to get TaxID from features
    tax_id = None
    for feature in record.features:
        if feature.type == "source":
            tax_id = feature.qualifiers.get("db_xref", [None])[0]
            if tax_id:
                if tax_id.startswith("taxon:"):
                    tax_id = tax_id.split(":")[1]
                    break

    if not tax_id:
        return {"Accession": accession, "superkingdom": None, "phylum": None,
                "class": None, "order": None, "family": None, "genus": None, "species": None}

    # Fetch taxonomy from Taxonomy database
    handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
    tax_record = Entrez.read(handle)[0]
    handle.close()

    # Build a dict of taxonomy
    lineage = {d["Rank"]: d["ScientificName"] for d in tax_record["LineageEx"]}

    result = {
        "Accession": accession,
        "superkingdom": lineage.get("superkingdom"),
        "phylum": lineage.get("phylum"),
        "class": lineage.get("class"),
        "order": lineage.get("order"),
        "family": lineage.get("family"),
        "genus": lineage.get("genus"),
        "species": tax_record.get("ScientificName")
    }

    return result

def make_blast_db(fasta_file, db_folder, title="db", dbtype="nucl"):

    fasta_file = Path(fasta_file)
    db_folder = Path(db_folder)
    db_folder.parent.mkdir(parents=True, exist_ok=True)

    is_gz = fasta_file.suffix == ".gz"
    is_windows = os.name == "nt"

    # -------------------------
    # WINDOWS: unzip first
    # -------------------------
    if is_windows and is_gz:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
            tmp_fasta = Path(tmp.name)

        with gzip.open(fasta_file, "rb") as f_in, open(tmp_fasta, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        cmd = [
            "makeblastdb",
            "-in", str(tmp_fasta),
            "-title", title,
            "-dbtype", dbtype,
            "-out", str(db_folder)
        ]

        subprocess.run(cmd, check=True)
        tmp_fasta.unlink()

    # -------------------------
    # UNIX: stream via zcat
    # -------------------------
    elif not is_windows and is_gz:
        cmd = (
            f"zcat < {fasta_file} | "
            f"makeblastdb -in - -title {title} -dbtype {dbtype} -out {db_folder}"
        )
        subprocess.run(cmd, shell=True, check=True)

    # -------------------------
    # Plain FASTA (all systems)
    # -------------------------
    else:
        cmd = [
            "makeblastdb",
            "-in", str(fasta_file),
            "-title", title,
            "-dbtype", dbtype,
            "-out", str(db_folder)
        ]
        subprocess.run(cmd, check=True)

def get_fasta_headers(fasta_file):
    fasta_file = Path(fasta_file)
    open_func = gzip.open if fasta_file.suffix == ".gz" else open

    headers = []
    with open_func(fasta_file, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip())

    return headers

def collect_accession_taxonomy(fasta_file, accession2taxonomy_parquet):
    TAX_COLS = [
        'Accession', 'superkingdom', 'phylum',
        'class', 'order', 'family', 'genus', 'species'
    ]

    if accession2taxonomy_parquet.exists():
        accession_df = pd.read_parquet(accession2taxonomy_parquet)
    else:
        accession_df = pd.DataFrame(columns=TAX_COLS)

    existing_accessions = set(accession_df['Accession'])

    fasta_accessions = set(get_fasta_headers(fasta_file))

    missing = fasta_accessions - existing_accessions
    cached = fasta_accessions & existing_accessions

    cached_df = accession_df[
        accession_df['Accession'].isin(cached)
    ]

    new_rows = []
    for acc in stqdm(missing, desc="Fetching taxonomy"):
        row = accession_to_taxonomy(acc)
        if row:
            new_rows.append(row)

    new_df = pd.DataFrame(new_rows, columns=TAX_COLS)

    out_df = pd.concat(
        [accession_df, new_df],
        ignore_index=True
    ).drop_duplicates(subset='Accession')

    return out_df

def check_ncbi():
    try:
        ncbi = NCBITaxa()
        ncbi.get_taxid_translator([65376])
        return True
    except Exception:
        return False

def check_blast_db(db_path):
    st.write('{}: Validating database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    try:
        subprocess.run(
            ["blastdbcmd", "-db", str(db_path), "-info"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )
        st.balloons()
        st.write('{}: Finished building database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    except:
        st.error('{}: Database could not be verified!'.format(datetime.datetime.now().strftime('%H:%M:%S')))

## DIATBARCODE

def run_diat_barcode(output_path, diat_barcode_xlsx):

    st.write('{} : Starting to collect accession numbers from .xlsx file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    diat_barcode_df = pd.read_excel(diat_barcode_xlsx, sheet_name='diatbarcode v12').fillna('')

    st.write('{} : Writing fasta file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    fasta_file = diat_barcode_xlsx.replace('.xlsx', '.fasta.gz').replace(' ', '_')
    with gzip.open(fasta_file, 'wt') as f:
        for line in diat_barcode_df[['Sequence ID', 'Sequence']].values.tolist():
            if line[0] != '':
                f.write(f'>{line[0]}\n')
                f.write(f'{line[1]}\n')

    st.write('{} : Finished to convert fasta format.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    st.write('{} : Starting to generate taxonomy file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    records = []
    for line in diat_barcode_df[["Species", "Genus", "Family (following Round, Crawford & Mann 1990)", "Order (following Round, Crawford & Mann 1990)", "Class (following Round, Crawford & Mann 1990)", "Phylum (following Algaebase 2018)", "Subkingdom (following Algaebase 2018)", "Sequence ID"]].values.tolist():
        if line[0] != '':
            records.append(line[::-1])
    records_df = pd.DataFrame(records, columns=['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    st.write('{} : Finished to convert accession numbers to taxonomy.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    # create tmp file
    temp_fasta = Path(Path(fasta_file).name.replace('.fasta.gz', '.fasta'))
    with gzip.open(fasta_file, 'rb') as f_in, open(temp_fasta, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # create db
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
    os.system(command)
    # delete tmp file
    if temp_fasta.exists():
        os.remove(temp_fasta)

    # move taxonomy file
    taxonomy_file = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    records_df.to_parquet(taxonomy_file)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    return db_folder

## MIDORI2

def midori2_taxonomy(fasta_file):
    all_accession_numbers = []
    with gzip.open(fasta_file, 'rt') as myfile:
        data = myfile.read()
        sequences = SeqIO.parse(StringIO(data), 'fasta')
        for record in sequences:
            accession = record.id
            taxonomy = []
            for i in record.description.split(';')[1:]:
                record_split = i.split('_')
                if len(record_split) == 2:
                    taxonomy.append(i.split('_')[0])
                else:
                    taxonomy.append(' '.join(i.split('_')[:2]))
            all_accession_numbers.append([accession] + taxonomy)

    accession_df = pd.DataFrame(all_accession_numbers, columns=['Accession', 'superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species'])
    return accession_df

def run_midori2(output_path, fasta_file):

    st.write('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = midori2_taxonomy(fasta_file)

    st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    # create tmp file
    temp_fasta = Path(Path(fasta_file).name.replace('.fasta.gz', '.fasta'))
    with gzip.open(fasta_file, 'rb') as f_in, open(temp_fasta, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # create db
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
    os.system(command)
    # delete tmp file
    if temp_fasta.exists():
        os.remove(temp_fasta)

    # write taxonomy file
    taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    accession_df.to_parquet(taxonomy_file_snappy)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    return db_folder

## PR2

def pr2_taxonomy(fasta_file):
    all_accession_numbers = []
    with gzip.open(fasta_file, 'rt') as myfile:
        data = myfile.read()
        sequences = SeqIO.parse(StringIO(data), 'fasta')
        for record in sequences:
            accession = record.id
            taxonomy_dict = {i.split(':')[0]:i.split(':')[1] for i in record.id.split(';')[1].split('=')[1].split(',')}
            taxonomy = []
            for taxon in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
                try:
                    taxonomy.append(taxonomy_dict[taxon].replace('_', ' '))
                except KeyError:
                    taxonomy.append('')
            all_accession_numbers.append([accession] + taxonomy)

    accession_df = pd.DataFrame(all_accession_numbers, columns=['Accession', 'superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species'])
    return accession_df

def run_pr2(output_path, fasta_file):

    st.write('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = pr2_taxonomy(fasta_file)

    st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    # create tmp file
    temp_fasta = Path(Path(fasta_file).name.replace('.fasta.gz', '.fasta'))
    with gzip.open(fasta_file, 'rb') as f_in, open(temp_fasta, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # create db
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
    os.system(command)
    # delete tmp file
    if temp_fasta.exists():
        os.remove(temp_fasta)

    # write taxonomy file
    taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    accession_df.to_parquet(taxonomy_file_snappy)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    return db_folder

## trl

def trnl_taxonomy(fasta_file, taxonomy_df):
    # Read sequences
    with gzip.open(fasta_file, 'rt') as handle:
        sequences = list(SeqIO.parse(handle, 'fasta'))
    accessions = [rec.id for rec in sequences]

    # Collect relevant accession numbers
    accession_df = taxonomy_df[taxonomy_df["Accession"].isin(accessions)]

    return accession_df

def run_trnl(output_path, fasta_file):

    st.write('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    taxonomy_df = pd.read_parquet(st.session_state['accession2taxonomy_parquet'])
    accession_df = trnl_taxonomy(fasta_file, taxonomy_df)

    st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    # create tmp file
    temp_fasta = Path(Path(fasta_file).name.replace('.fasta.gz', '.fasta'))
    with gzip.open(fasta_file, 'rb') as f_in, open(temp_fasta, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # create db
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
    os.system(command)
    # delete tmp file
    if temp_fasta.exists():
        os.remove(temp_fasta)

    # write taxonomy file
    taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    accession_df.to_parquet(taxonomy_file_snappy)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    return db_folder

## UNITE

def unite_taxonomy(fasta_file):
    all_accession_numbers = []
    with gzip.open(fasta_file, 'rt') as myfile:
        data = myfile.read()
        sequences = SeqIO.parse(StringIO(data), 'fasta')
        for record in sequences:
            accession = record.id
            taxonomy_dict = {i.split('__')[0]:i.split('__')[1] for i in record.id.split('|')[4].split(';')}
            taxonomy = []
            for taxon in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
                try:
                    taxonomy.append(taxonomy_dict[taxon].replace('_', ' '))
                except KeyError:
                    taxonomy.append('')
            all_accession_numbers.append([accession] + taxonomy)

    accession_df = pd.DataFrame(all_accession_numbers, columns=['Accession', 'superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species'])
    return accession_df

def run_unite(output_path, fasta_file):

    st.write('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    ## collect accession numbers
    accession_df = unite_taxonomy(fasta_file)

    st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    # create tmp file
    temp_fasta = Path(Path(fasta_file).name.replace('.fasta.gz', '.fasta'))
    with gzip.open(fasta_file, 'rb') as f_in, open(temp_fasta, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # create db
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
    os.system(command)
    # delete tmp file
    if temp_fasta.exists():
        os.remove(temp_fasta)

    # write taxonomy file
    taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    accession_df.to_parquet(taxonomy_file_snappy)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    return db_folder

## SILVA

def build_taxonomy_lookup(taxids, ncbi):
    TARGET_RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxids = list({int(t) for t in taxids if pd.notna(t)})
    lineages = {}
    for t in taxids:
        try:
            lineages[t] = ncbi.get_lineage(t)
        except Exception:
            lineages[t] = None
    all_lineage_taxids = {tid for lin in lineages.values() if lin for tid in lin}
    ranks = ncbi.get_rank(list(all_lineage_taxids)) if all_lineage_taxids else {}
    names = ncbi.get_taxid_translator(list(all_lineage_taxids)) if all_lineage_taxids else {}
    return lineages, ranks, names

def silva_taxonomy(taxid_df, acc_col='Accession', taxid_col='taxid', ncbi=ncbi):
    TARGET_RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    columns = ['Accession'] + TARGET_RANKS
    taxid_series = pd.to_numeric(taxid_df[taxid_col], errors='coerce')
    lineages, ranks, names = build_taxonomy_lookup(taxid_series.dropna().unique(), ncbi)
    @lru_cache(maxsize=None)
    def resolve(taxid):
        lineage = lineages.get(taxid)
        if not lineage:
            return tuple(f'Unknown {r}' for r in TARGET_RANKS)
        rank_map = {ranks[t]: names[t] for t in lineage if t in ranks and t in names}
        return tuple(rank_map.get(r, '') for r in TARGET_RANKS)
    records = [
        (acc, *(resolve(int(taxid)) if pd.notna(taxid) else tuple(f'Unknown {r}' for r in TARGET_RANKS)))
        for acc, taxid in zip(taxid_df[acc_col], taxid_series)
    ]
    return pd.DataFrame.from_records(records, columns=columns)

def run_silva(output_path, fasta_file, selected_tax_file):

    ## collect taxonomy
    st.write('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    taxid_df = pd.read_csv(selected_tax_file, sep="\t", header=None)
    taxid_df.columns = ['Accession', 'taxid']
    accession_df = silva_taxonomy(taxid_df)

    ## collect accession numbers
    st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # create database
    db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    # create tmp file
    temp_fasta = Path(Path(fasta_file).name.replace('.fasta.gz', '.fasta'))
    with gzip.open(fasta_file, 'rb') as f_in, open(temp_fasta, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # create db
    if not os.path.isdir(db_folder):
        os.mkdir(db_folder)
    db_folder = db_folder.joinpath('db')
    command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
    os.system(command)
    # delete tmp file
    if temp_fasta.exists():
        os.remove(temp_fasta)

    # write taxonomy file
    taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
    accession_df.to_parquet(taxonomy_file_snappy)

    ## zip the folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    return db_folder

# CUSTOM APSCALE DB

def custom_db_taxonomy(fasta_file, taxonomy_df):
    # Read sequences
    with open(fasta_file, 'r') as handle:
        sequences = list(SeqIO.parse(handle, 'fasta'))

    taxonomy_res = []
    for record in sequences:
        desc = record.description
        if desc.startswith('Taxonomy;;'):
            taxonomy = desc.split(';;')
            taxonomy_res.append(taxonomy)
        elif desc.startswith('Accession;;'):
            accession = desc.split(';;')[1]
            taxonomy = taxonomy_df[taxonomy_df["Accession"]==accession].head(1).values.tolist()[0]
            taxonomy_res.append(taxonomy)
    accession_df = pd.DataFrame(taxonomy_res, columns=taxonomy_df.columns)
    return accession_df

def run_custom_db(output_path, input_file, suffix):

    if suffix == '.xlsx':
        st.write('{} : Starting to collect accession numbers from .xlsx file.'.format(
            datetime.datetime.now().strftime('%H:%M:%S')))

        custom_barcodes_df = pd.read_excel(input_file).fillna('')
        records_df = custom_barcodes_df[['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]

        # check if accession numbers are unique
        if records_df['Accession'].duplicated().any() == True:
            st.error('{} : Please use unique accession numbers!'.format(
                datetime.datetime.now().strftime('%H:%M:%S')))
            return

        st.write('{} : Writing fasta file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

        fasta_file = str(input_file).replace('.xlsx', '.fasta.gz').replace(' ', '_')
        with gzip.open(fasta_file, 'wt') as f:
            for line in custom_barcodes_df[['Accession', 'sequence']].values.tolist():
                if line[0] != '':
                    f.write(f'>{line[0]}\n')
                    f.write(f'{line[1]}\n')

        st.write('{} : Finished to convert fasta format.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

        st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        # create database
        db_name = Path(fasta_file).name.replace('.fasta.gz', '').replace('.', '_')
        db_folder = Path(output_path).joinpath(f'db_{db_name}')
        # create tmp file
        temp_fasta = Path(Path(fasta_file).name.replace('.fasta.gz', '.fasta'))
        with gzip.open(fasta_file, 'rb') as f_in, open(temp_fasta, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        # create db
        if not os.path.isdir(db_folder):
            os.mkdir(db_folder)
        db_folder = db_folder.joinpath('db')
        command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
        os.system(command)
        # delete tmp file
        if temp_fasta.exists():
            os.remove(temp_fasta)

        # move taxonomy file
        taxonomy_file = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
        records_df.to_parquet(taxonomy_file)

        ## zip the folder
        output = Path(output_path).joinpath(f'db_{db_name}')
        shutil.make_archive(output, 'zip', output)

        st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

        return db_folder

    elif suffix == '.fasta':
        st.write('{} : Starting to collect taxonomy from fasta files.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

        ## collect accession numbers
        taxonomy_df = pd.read_parquet(st.session_state['accession2taxonomy_parquet'])
        accession_df = custom_db_taxonomy(input_file, taxonomy_df)

        st.write('{} : Starting to create database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

        # create database
        db_name = Path(input_file).name.replace('.fasta', '').replace('.', '_')
        db_folder = Path(output_path).joinpath(f'db_{db_name}')
        # create db
        if not os.path.isdir(db_folder):
            os.mkdir(db_folder)
        db_folder = db_folder.joinpath('db')
        command = f'makeblastdb -in {input_file} -title db -dbtype nucl -out {db_folder}'
        os.system(command)

        # write taxonomy file
        taxonomy_file_snappy = db_folder.parent.joinpath('db_taxonomy.parquet.snappy')
        accession_df.to_parquet(taxonomy_file_snappy)

        ## zip the folder
        output = Path(output_path).joinpath(f'db_{db_name}')
        shutil.make_archive(output, 'zip', output)

        st.write('{} : Finished database creation.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

        return db_folder

    else:
        st.warning('Please import an uncompressed .fasta or .xlsx file!')

# BOLDsystem DB
def run_BOLD(output_path, parquet_file):

    st.write(f'{datetime.datetime.now():%H:%M:%S} : Reading sequences from BOLD snapshot.')

    # only pull the columns we actually use out of the parquet file (column pruning) --
    # this cuts I/O and memory since BOLD snapshots typically carry many more columns than this
    columns = ['processid', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'nuc']
    taxonomy_columns = [c for c in columns if c != 'nuc']

    batch_size = 10_000
    pf = pq.ParquetFile(parquet_file)
    n_batches = math.ceil(pf.metadata.num_rows / batch_size)

    # set up output locations up front so we can stream results straight to disk
    db_name = Path(parquet_file).name.replace('.parquet', '').replace('.', '_')
    db_folder = Path(output_path).joinpath(f'db_{db_name}')
    os.makedirs(db_folder, exist_ok=True)
    taxonomy_file_snappy = db_folder.joinpath('db_taxonomy.parquet.snappy')
    temp_fasta = Path(Path(parquet_file).name.replace('.parquet', '.fasta'))

    n_sequences = 0
    n_dropped_no_species = 0
    parquet_writer = None
    try:
        with open(temp_fasta, 'w') as f_out:
            for batch in stqdm(pf.iter_batches(batch_size=batch_size, columns=columns),
                               total=n_batches, desc="Reading & processing BOLD snapshot"):
                chunk_df = batch.to_pandas().fillna("")
                n_sequences += len(chunk_df)

                has_species = chunk_df['species'].str.strip() != ''
                n_dropped_no_species += int((~has_species).sum())
                chunk_df = chunk_df.loc[has_species]

                if chunk_df.empty:
                    continue

                # stream this batch's taxonomy columns straight into the output parquet file
                taxonomy_table = pa.Table.from_pandas(
                    chunk_df[taxonomy_columns].rename(columns={'processid': 'Accession', 'kingdom': 'superkingdom'}),
                    preserve_index=False
                )
                if parquet_writer is None:
                    # schema is only known once we see the first (non-empty) batch
                    parquet_writer = pq.ParquetWriter(taxonomy_file_snappy, taxonomy_table.schema, compression='snappy')
                parquet_writer.write_table(taxonomy_table)

                # vectorized fasta building (pandas .str ops) instead of a per-row Python loop --
                # same filter as before (strip gap chars, keep sequences longer than 50bp)
                sequences = chunk_df['nuc'].str.strip('-')
                keep = sequences.str.len() > 50
                fasta_lines = '>' + chunk_df.loc[keep, 'processid'] + '\n' + sequences[keep] + '\n'
                f_out.writelines(fasta_lines)
    finally:
        if parquet_writer is not None:
            parquet_writer.close()

    print(f'{datetime.datetime.now():%H:%M:%S}: Dropped {n_dropped_no_species:,} records with no species-level taxonomy.')
    st.write(f'{datetime.datetime.now():%H:%M:%S} : Collected {n_sequences:,} sequences from BOLD snapshot.')
    st.write(f'{datetime.datetime.now():%H:%M:%S} : Starting to create database.')

    # build the blast db from the fasta we just streamed to disk
    db_folder = db_folder.joinpath('db')
    command = f'makeblastdb -in {temp_fasta} -title db -dbtype nucl -out {db_folder}'
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        # os.system() swallowed failures silently -- surface them instead
        st.write(f'{datetime.datetime.now():%H:%M:%S} : makeblastdb failed with exit code {result.returncode}.')

    # delete tmp fasta now that the blast db has been built from it
    if temp_fasta.exists():
        os.remove(temp_fasta)

    # zip the finished db folder
    output = Path(output_path).joinpath(f'db_{db_name}')
    shutil.make_archive(output, 'zip', output)

    st.write(f'{datetime.datetime.now():%H:%M:%S} : Finished database creation.')

    return db_folder

### GENERAL FUNCTIONS

def zip_to_gz(fasta_file):
    # Determine the output .gz file path
    gz_path = os.path.splitext(fasta_file)[0] + '.gz'

    temp_dir = 'temp_unzip'
    os.makedirs(temp_dir, exist_ok=True)

    # Unzip the contents
    with zipfile.ZipFile(fasta_file, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    # Gzip the contents
    with open(gz_path, 'wb') as gz_file:
        for root, _, files in os.walk(temp_dir):
            for file in files:
                file_path = os.path.join(root, file)
                with open(file_path, 'rb') as f_in:
                    with gzip.open(gz_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

    # Clean up the temporary directory
    shutil.rmtree(temp_dir)
    os.remove(fasta_file)

    return gz_path

def open_folder(folder_path):
    # Get the current operating system
    current_os = platform.system()

    # Open the folder based on the OS
    try:
        if current_os == "Windows":
            subprocess.Popen(f'explorer "{folder_path}"')
        elif current_os == "Darwin":  # macOS
            subprocess.Popen(['open', folder_path])
        else:  # Linux
            subprocess.Popen(['xdg-open', folder_path])
    except Exception as e:
        st.write(f"Failed to open folder: {e}")

def open_file(path: Path):
    path = str(path)
    if platform.system() == "Darwin":
        subprocess.run(["open", path])
    elif platform.system() == "Windows":
        os.startfile(path)
    else:
        subprocess.run(["xdg-open", path])

def validation():

    all_dbs = [Path(i) for i in glob.glob('/Volumes/Coruscant/APSCALE_projects/APSCALE_databases_input/*.fasta.gz')]
    output_fasta = Path('/Volumes/Coruscant/APSCALE_projects/APSCALE_databases_input/combined_first_reads.fasta')
    max_seqs_per_file = 4
    with open(output_fasta, "w") as out_f:
        for file in all_dbs:
            name = file.name.replace(".fasta.gz", "")
            with gzip.open(file, "rt") as handle:
                for i, record in enumerate(SeqIO.parse(handle, "fasta"), start=1):
                    if i > max_seqs_per_file:
                        break
                    record.id = f"{name}_{i}"
                    record.name = ""
                    record.description = ""
                    SeqIO.write(record, out_f, "fasta")

    query = '/Volumes/Coruscant/APSCALE_projects/APSCALE_databases_input/combined_first_reads.fasta'
    all_dbs = glob.glob('/Volumes/Coruscant/APSCALE_projects/APSCALE_databases/*/db_taxonomy.parquet.snappy')
    for db in all_dbs:
        db_path = Path(db).parent
        out_path = Path(query).with_name(Path(db).parent.name)
        command = f'apscale_blast -q {query} -db {db_path} -o {out_path}'
        os.system(command)

### STREAMLIT FUNCTIONS

def create_dropdown(title, description, icon, filetype, run_function, working_dir, input_folder, weblink, downloadlink, example_file_name):

    with st.expander(f"{icon} {title}", expanded=False):
        st.markdown(f"**{description}**")
        st.write("")

        col1, col2 = st.columns(2)

        with col1:
            st.write('Weblinks:')
            if st.button(f"🔍 Learn more about {title}", key=f"learn_{title}", width='stretch'):
                webbrowser.open(weblink)
            if st.button(f"🔗 Download link to latest version", key=f"download_{title}", width='stretch'):
                webbrowser.open(downloadlink)
            st.info(f'Example file:\n\n{example_file_name}')

        with col2:
            st.write('')
            # load real files WITHOUT modifying filenames
            file_names = sorted(glob.glob(str(input_folder / f'*{filetype}*')))
            file_names_dict = {Path(i).name: Path(i) for i in file_names}
            selected_file_key = st.selectbox(f"Select a database file ({filetype}):", list(file_names_dict.keys()), key=f"{title}_db_file")

            # create database
            if selected_file_key != None:
                selected_file = file_names_dict[selected_file_key]
                if st.button(f"⚙️ Create {title} database", key=f"create_{title}", width='stretch', type='primary'):
                    if selected_file.suffix == '.zip':
                        st.write('{}: Converting .zip to .gz format.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
                        selected_file = zip_to_gz(selected_file)
                    db_path = run_function(working_dir, str(selected_file))
                    check_blast_db(db_path)
        st.write("___")

def create_dropdown_silva(title, description, icon, filetype, run_function, working_dir, input_folder, weblink, downloadlink, example_file_name, example_file_name_tax):

    with st.expander(f"{icon} {title}", expanded=False):
        st.markdown(f"**{description}**")
        st.write("")

        col1, col2 = st.columns(2)

        with col1:
            st.write('Weblinks:')
            if st.button(f"🔍 Learn more about {title}", key=f"learn_{title}", width='stretch'):
                webbrowser.open(weblink)
            if st.button(f"🔗 Download link to latest version", key=f"download_{title}", width='stretch'):
                webbrowser.open(downloadlink)
            if st.button(f"🔗 Download link to latest taxonomy", key=f"download_{title}_tax", width='stretch'):
                webbrowser.open('https://www.arb-silva.de/current-release/Exports/taxonomy')
            st.info(f'Example file:\n\n{example_file_name}\n\n')
            st.info(f'Example taxonomy file:\n\n{example_file_name_tax}\n\n')

        with col2:
            st.write('')
            # load real files WITHOUT modifying filenames
            file_names = sorted(glob.glob(str(input_folder / f'*{filetype}*')))
            file_names_dict = {Path(i).name: Path(i) for i in file_names}
            selected_file_key = st.selectbox(f"Select a database file ({filetype}):", list(file_names_dict.keys()), key=f"{title}_db_file")

            # load real files WITHOUT modifying filenames
            file_names_tax = sorted(glob.glob(str(input_folder / f'tax_slv_ssu_*.acc_taxid.gz')))
            file_names_dict_tax = {Path(i).name: Path(i) for i in file_names_tax}
            selected_tax_file_key = st.selectbox("Select a taxonomy file (.txt.gz):", list(file_names_dict_tax.keys()), key=f"{title}_tax_file")

            # create database
            if selected_file_key != None and selected_tax_file_key != None:

                selected_file = file_names_dict[selected_file_key]
                selected_tax_file = file_names_dict_tax[selected_tax_file_key]
                if st.button(f"⚙️ Create {title} database", key=f"create_{title}", width='stretch', type='primary'):
                    if selected_file.suffix == '.zip':
                        st.write('{}: Converting .zip to .gz format.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
                        selected_file = zip_to_gz(selected_file)
                    db_path = run_function(working_dir, str(selected_file), str(selected_tax_file))
                    check_blast_db(db_path)
        st.write("___")

def create_dropdown_custom(title, description, icon, filetype, run_function, working_dir, input_folder):

    with st.expander(f"{icon} {title}", expanded=False):
        st.markdown(f"**{description}**")
        st.write("")

        col1, col2 = st.columns(2)

        with col1:
            st.write('Weblinks:')
            if st.button(f"🔍 Learn more about {title}", key=f"learn_{title}", width='stretch'):
                webbrowser.open("https://github.com/TillMacher/apscale_blast")
            if st.button(f"📝 Create Template Excel File", key=f"excel_example_{title}", width='stretch'):
                row = [['JFKC01000052.336.3179',
                          'Bacteria',
                          'Pseudomonadota',
                          'Alphaproteobacteria',
                          'Rhodobacterales',
                          'Roseobacteraceae',
                          'Marivita',
                          'Marivita geojedonensis',
                           'ACGT']]
                template_df = pd.DataFrame(row, columns=['Accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'sequence'])
                template_xlsx = input_folder / 'custom_import_template.xlsx'
                template_df.to_excel(template_xlsx, index=False)
                open_file((template_xlsx))
            if st.button(f"📝 Create Template FASTA File", key=f"fasta_example_{title}", width='stretch'):
                template_fasta = input_folder / 'custom_import_template.fasta'
                f = open(template_fasta, "w")
                # NCBI accession number
                f.write(f'>Accession;;JFKC01000052.336.3179\n')
                f.write(f'ACGT\n')
                # Custom taxonomy
                f.write('>Taxonomy;;Bacteria;;Pseudomonadota;;Alphaproteobacteria;;Rhodobacterales;;Roseobacteraceae;;Marivita;;Marivita geojedonensis\n')
                f.write(f'ACGT\n')
                f.close()
                open_file(template_fasta)

        with col2:
            st.write('')
            # load real files WITHOUT modifying filenames
            file_names = sorted(glob.glob(str(input_folder / f'*{filetype}*')))
            file_names_dict = {Path(i).name: Path(i) for i in file_names}
            selected_file_key = st.selectbox("Select a database file:", list(file_names_dict.keys()), key=f"{title}_db_file")

            # create database
            if selected_file_key != None:
                selected_file = file_names_dict[selected_file_key]
                if st.button(f"⚙️ Create {title}", key=f"create_{title}", width='stretch', type='primary'):
                    db_path = run_custom_db(working_dir, str(selected_file), selected_file.suffix)
                    check_blast_db(db_path)
        st.write("___")

########################################################################################################################
########################################################################################################################
### GUI ###
st.title("🗃️ APSCALE-blast Database Creator")
st.markdown("A simple tool to build reference databases for APSCALE-blast.")
st.write("___")

########################################################################################################################
colored_header(
    label="Select an APSCALE working directory",
    description="Databases will be written to the APSCALE_databases folder.",
    color_name="blue-70")

col1, col2 = st.columns(2)
working_dir = False
input_folder = Path('.')
with col1:
    st.text_input('Please provide a PATH to your APSCALE project folder:', key='output_path')
with col2:
    output_path = Path(st.session_state['output_path'])
    st.write('Chosen working directory:')
    if output_path != Path('.'):
        working_dir = output_path / "APSCALE_databases"
        if working_dir.exists():
            st.success('The PATH exists and you are ready to create databases!')
            input_folder = output_path / 'APSCALE_databases_input'
            if not input_folder.exists():
                if st.button('Create APSCALE database input folder.'):
                    os.makedirs(input_folder, exist_ok=True)
        else:
            st.warning('The provided working directory does not exists!')
    else:
        st.warning('Please provide a PATH.')

if st.button("🔄 Refresh files and folders", width='stretch'):
    st.rerun()

########################################################################################################################
colored_header(
    label="NCBITaxa Input",
    description="Some Modules Require NCBI taxonomy.",
    color_name="blue-70")

# --- Email input ---
usermail = st.text_input(label="📧 Enter your email address for NCBI (required):", key="usermail", placeholder="your.name@email.com")
if usermail:
    Entrez.email = usermail

# --- Check NCBI taxonomy DB ---
db_available = check_ncbi()
if db_available:
    st.success("NCBI taxonomy database is available and ready to use ✅")
else:
    st.error("NCBI taxonomy database not found or corrupted ❌")
    st.info(
        "This database is required for taxonomic lookups.\n\n"
        "📦 Download size: ~150–300 MB\n"
        "💾 Disk usage after extraction: ~300–600 MB\n\n"
        "The download is performed once and reused automatically."
            )
    if st.button("⬇️ Download NCBI Taxonomy Database"):
        if not usermail:
            st.warning("Please provide your email address before downloading (NCBI requirement).")
        else:
            with st.spinner("Downloading taxonomy database... this may take a few minutes ⏳"):
                ncbi = NCBITaxa()
                ncbi.update_taxonomy_database()
            st.success("Download complete ✅")
            st.rerun()
# --- Accession → taxonomy file check ---
if working_dir != False and input_folder.exists():
    accession2taxonomy_parquet = path_to_apscale_blast / "accession_taxonomy.parquet.snappy"
    st.session_state['accession2taxonomy_parquet'] = accession2taxonomy_parquet
    if accession2taxonomy_parquet.exists():
        st.success("Accession-to-taxonomy mapping file found ✅")
    else:
        st.warning(
            "No accession-to-taxonomy mapping file found.\n\n"
            "All accession numbers will be resolved from scratch, which may take longer.")

########################################################################################################################
colored_header(
    label="Create your own reference database",
    description="Choose a curated reference source and build your custom BLAST database.",
    color_name="blue-70")

if working_dir == False:
    st.warning('Please select a working directory first.')
else:

    col1, col2 = st.columns(2)
    with col1:
        if st.button(label='📥 Open APSCALE dbcreator Input Folder', width='stretch'):
            open_folder(input_folder)
    with col2:
        if st.button(label='📤 Open APSCALE Database Folder', width='stretch'):
            open_folder(working_dir)

    if Path(working_dir).exists() and Path(input_folder).exists():
        create_dropdown("DiatBarcode",
                        "A curated diatom reference database for freshwater biomonitoring.",
                        "🟩",
                        ".xlsx",
                        run_diat_barcode,
                        working_dir,
                        input_folder,
                        "https://carrtel-collection.hub.inrae.fr/barcoding-databases/diat.barcode/database-download",
                        "https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/TOMBYZ",
                        "2024-05-28-Diat.barcode_release-version 12.4.xlsx"
                        )

        create_dropdown("MIDORI2",
                        "Extensive mitochondrial reference database for metazoan biodiversity.",
                        "🟦",
                        ".fasta",
                        run_midori2,
                        working_dir,
                        input_folder,
                        "https://www.reference-midori.info/",
                        "https://www.reference-midori.info/download.php",
                        "MIDORI2_UNIQ_NUC_GB266_srRNA_BLAST.fasta.zip"
                        )

        create_dropdown("PR2",
                        "Protist Ribosomal Reference database for 18S metabarcoding.",
                        "🟪",
                        ".fasta",
                        run_pr2,
                        working_dir,
                        input_folder,
                        "https://pr2-database.org/",
                        "https://github.com/pr2database/pr2database/releases/tag/v5.1.0.0",
                        "pr2_version_5.1.0_SSU_UTAX.fasta.gz"
                        )

        create_dropdown_silva("SILVA",
                        "Comprehensive ribosomal RNA database (16S/18S).",
                        "🟫",
                        ".fasta",
                        run_silva,
                        working_dir,
                        input_folder,
                        "https://www.arb-silva.de/",
                        "https://www.arb-silva.de/current-release/Exports",
                        "SILVA_144_SSURef_tax_silva_trunc.fasta.gz",
                        "tax_slv_ssu_144.acc_taxid.gz"
                        )

        create_dropdown("trnL",
                        "Plant chloroplast trnL reference database.",
                        "🌿",
                        ".fasta",
                        run_trnl,
                        working_dir,
                        input_folder,
                        "https://ucedna.com/reference-databases-for-metabarcoding",
                        "https://ucla.app.box.com/s/n4hbocub5pdoffdtcfv1tdjr7cl4cc1n",
                        "trnL.fasta.gz"
                        )

        create_dropdown("UNITE",
                        "Fungal ITS reference database.",
                        "🍄",
                        ".fasta",
                        run_unite,
                        working_dir,
                        input_folder,
                        "https://unite.ut.ee/",
                        "https://unite.ut.ee/repository.php",
                        "UNITE_sh_general_release_dynamic_s_all_eukaryotes_19.02.2025.fasta.gz"
                        )

        create_dropdown("BOLDsystems v5",
                        "BOLD DNA Barcode Reference Library.",
                        "🟧",
                        ".parquet",
                        run_BOLD,
                        working_dir,
                        input_folder,
                        "https://www.boldsystems.org/data/data-packages/",
                        "https://bench.boldsystems.org/index.php/datapackages/Latest",
                        "BOLD_Public.30-Jun-2026.parquet"
                        )

        create_dropdown_custom("Custom APSCALE database",
                        "Import your own sequences from .fasta or .xlsx files.",
                        "🧬",
                        "",
                        run_custom_db,
                        working_dir,
                        input_folder,
                        )