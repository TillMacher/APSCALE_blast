import argparse
import multiprocessing
import os
import sys
import datetime
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import subprocess
import importlib
from ete3 import NCBITaxa, Tree
# APSCALE blast
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter
# BOLDigger3
from boldigger3 import id_engine
from boldigger3 import add_metadata
from boldigger3 import select_top_hit
from boldigger3 import download_database

def organism_filter(name):
    try:
        name = int(name)
    except ValueError:
        pass
    if isinstance(name, int):
        try:
            # Initialize NCBITaxa database
            ncbi = NCBITaxa()
            res = ncbi.get_taxid_translator([name])
            organism = f'{res[name]} (taxid:{name})'
            return organism
        except:
            return f'{name} not found!'
    else:
        try:
            # Initialize NCBITaxa database
            ncbi = NCBITaxa()
            res = ncbi.get_name_translator([name])
            taxid = res[name][-1]
            organism = f'{name} (taxid:{taxid})'
            return organism
        except:
            return f'{name} not found!'

def issue_list_flag(merged_df, issue_list):

    # Merge on shared columns
    merged_with_issues = merged_df.merge(
        issue_list,
        on=["Species", "Family", "Phylum"],
        how="left",  # keep all rows from merged_df
        suffixes=("", "_issue")  # avoid name conflicts
    )

    merged_with_issues = merged_with_issues.fillna('')
    col = "Status"  # column to move
    merged_with_issues = merged_with_issues[[c for c in merged_with_issues.columns if c != col] + [col]]

    return merged_with_issues

def merge_reblast(taxonomy_df, taxonomy_reblast_df, db_name):
    # Example key column (adjust as needed)
    key_col = "unique ID"

    # Merge both tables
    merged = taxonomy_df.merge(
        taxonomy_reblast_df,
        on=key_col,
        suffixes=("", "_reblast"),
        how="outer"  # keep all sequences
    )

    # Compare similarity columns
    def choose_row(row):
        # Case 1: Only one of the two has data
        if pd.isna(row["Similarity"]):
            row["Similarity"] = row["Similarity_reblast"]
            row["Reblast result"] = True
        elif pd.isna(row["Similarity_reblast"]):
            row["Reblast result"] = ''
        # Case 2: Both exist — compare
        elif row["Similarity_reblast"] > row["Similarity"]:
            for col in taxonomy_reblast_df.columns:
                if col != key_col:
                    row[col] = row[f"{col}_reblast"]
            row["Reblast result"] = True
        else:
            row["Reblast result"] = ''
        return row

    merged = merged.apply(choose_row, axis=1)

    # Drop the reblast columns (clean up)
    merged = merged[[c for c in merged.columns if not c.endswith("_reblast")]]

    return merged

def run_boldigger3(query_fasta):

    # Collect input values
    thresholds = [97, 95, 90, 85, 75, 50]
    query_fasta = str(query_fasta)
    mode = 3 # exhaustive Search
    db = 2 # ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)

    # download the current metadata from BOLD
    metadata_download.main()

    # run the id engine
    id_engine.main(
        query_fasta,
        db,
        mode,
    )

    # add additional data via the metadata
    add_metadata.main(query_fasta)

    # select the top hit
    select_top_hit.main(query_fasta, thresholds)

def merge_boldigger3(taxonomy_df, boldigger_df, db_name):
    # Example key column (adjust as needed)
    key_col = "unique ID"

    boldigger_df.columns = [i.capitalize() for i in boldigger_df.columns.tolist()]
    boldigger_df = boldigger_df.rename(columns={"Id": "unique ID", "Pct_identity": "Similarity"})

    # Merge both tables
    merged = taxonomy_df.merge(
        boldigger_df,
        on=key_col,
        suffixes=("", "_reblast"),
        how="outer"  # keep all sequences
    )

    # Compare similarity columns
    def choose_row(row):
        # Case 1: Only one of the two has data
        if pd.isna(row["Similarity"]):
            row["Similarity"] = row["Similarity_reblast"]
            row["Reblast result"] = True
        elif pd.isna(row["Similarity_reblast"]):
            row["Reblast result"] = ''
        # Case 2: Both exist — compare
        elif row["Similarity_reblast"] > row["Similarity"]:
            for col in boldigger_df.columns:
                reblast_col = f"{col}_reblast"
                if col != key_col and reblast_col in merged.columns:
                    row[col] = row[reblast_col]
            row["Reblast result"] = True
        else:
            row["Reblast result"] = ''
        return row

    merged = merged.apply(choose_row, axis=1)
    merged = merged.fillna('')

    # Drop the reblast columns (clean up)
    merged = merged[[c for c in merged.columns if not c.endswith("_reblast")]]

    return merged

def check_blast_db(db_path):
    print('{}: Validating database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    try:
        subprocess.run(
            ["blastdbcmd", "-db", str(db_path), "-info"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )
        print('{}: Your database looks fine.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        return True
    except:
        print('{}: Database could not be verified!'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        return False

def get_package_versions(pkg):
    try:
        version = importlib.metadata.version(pkg)
        return version
    except importlib.metadata.PackageNotFoundError:
        return

import pandas as pd

def settings_to_df(args):
    settings = {
        # Main settings
        'database':         args.database,
        'query_fasta':      args.query_fasta,
        'out':              args.out,
        # BLASTn settings
        'n_cores':          args.n_cores,
        'task':             args.task,
        'subset_size':      args.subset_size,
        'max_target_seqs':  args.max_target_seqs,
        # Filter settings
        'thresholds':       args.thresholds,
        'filter':           args.filter,
        'masking':          args.masking,
        'rating_range':     args.rating_range,
        'sim_range':        args.sim_range,
        # Re-BLAST settings
        'reblast_db':       args.reblast_db,
        'reblast_sim':      args.reblast_sim,
        # Misc
        'blastn_exe':       args.blastn_exe,
    }

    df = pd.DataFrame(settings.items(), columns=['Setting', 'Value'])
    return df

def main():
    """
    APSCALE BLASTn suite
    Command-line tool to run and filter BLASTn searches.
    """

    apscale_blast_version = get_package_versions('apscale_blast')

    # Introductory message with usage examples
    message = f"""
    APSCALE blast command line tool - {apscale_blast_version}
    Example commands:
    $ apscale_blast -h
    
    Basic blast 
    $ apscale_blast -db ./MIDORI2_UNIQ_NUC_GB266_srRNA_BLAST -q ./12S_apscale_ESVs.fasta -o ./12S_apscale_ESV
    
    For large datasets or large databases (e.g., MIDORI COI) use the task "megablast"
    $ apscale_blast -db ./db_MIDORI2_UNIQ_NUC_GB266_CO1_BLAST -q ./COI_apscale_ESVs.fasta -o ./COI_apscale_ESV -task megablast
    
    Database creator (streamlit GUI)
    $ apscale_blast dbcreator
    
    """
    print(message)

    if sys.argv[1] == "dbcreator":
        # Determine theme
        theme = "dark"  # default
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            script_path = os.path.join(script_dir, 'c_dbcreator.py')

            print(f"Starting APSCALE-blast dbcreator in {theme} mode. Press CTRL + C to exit.")
            subprocess.run([
                'streamlit', 'run', script_path,
                '--theme.base', theme,
                '--server.address', 'localhost'
            ])
        except KeyboardInterrupt:
            print("Exiting...")
            sys.exit()

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=f'APSCALE BLAST {apscale_blast_version}')

    # === Main settings ===
    main_settings = parser.add_argument_group("Main settings")
    main_settings.add_argument('-database', '-db', type=str, required=True,
                               help='Path to the local reference database.')
    main_settings.add_argument('-query_fasta', '-q', type=str, required=True,
                               help='Path to the query FASTA file.')
    main_settings.add_argument('-out', '-o', type=str, default='./blastn',
                               help='Output directory. A new folder will be created here. [Default: ./blastn]')

    # === BLASTn settings ===
    blastn_settings = parser.add_argument_group("BLASTn settings")
    blastn_settings.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 2,
                                 help='Number of CPU cores to use. [Default: CPU count - 2]')
    blastn_settings.add_argument('-task', type=str, default='megablast',
                                 help='BLASTn task: blastn, megablast, or dc-megablast. [Default: megablast]')
    blastn_settings.add_argument('-subset_size', type=int, default=100,
                                 help='Number of sequences per FASTA subset. [Default: 100]')
    blastn_settings.add_argument('-max_target_seqs', type=int, default=20,
                                 help='Maximum number of hits retained per query. [Default: 20]')

    # === Filter settings ===
    filter_settings = parser.add_argument_group("Filter settings")
    filter_settings.add_argument('-thresholds', type=str, default='97,95,90,87,85',
                                 help='Taxonomic similarity thresholds (comma-separated). [Default: 97,95,90,87,85]')
    filter_settings.add_argument('-filter', type=str, default='2',
                                 help='Filtering mode: 1 = e-value → similarity, 2 = similarity → e-value, 3 = similarity [Default: 2]')
    filter_settings.add_argument('-masking', action='store_false',
                                 help='Disable masking. Masking is enabled by default.')
    filter_settings.add_argument('-rating_range', type=int, default=5,
                                 help='Range of allowed rating values. [Default: 5]')
    filter_settings.add_argument('-sim_range', type=int, default=0,
                                 help='Value that is subtracted from the highest similarity to retain hits. [Default: 0]')

    # === Re-BLAST settings ===
    reblast_settings = parser.add_argument_group("Re-BLAST settings")
    reblast_settings.add_argument('-reblast_db', '-db2', type=str,
                                  help='Path to an optional secondary database for re-BLAST. Choose "boldigger" for boldigger3 assignment.')
    reblast_settings.add_argument('-reblast_sim', type=int, default=98,
                                  help='Similarity threshold for re-BLASTing hits against db2. [Default: 98]')

    # === Misc settings ===
    misc_settings = parser.add_argument_group("Misc")
    misc_settings.add_argument('-blastn_exe', type=str, default='blastn',
                               help='Path to the BLASTn executable. [Default: blastn]')

    # Parse the arguments
    args = parser.parse_args()

    # Handle missing arguments interactively for both commands
    if not args.database and not args.query_fasta:
        args.database = input("Please enter PATH to database: ").strip('"')
        args.query_fasta = input("Please enter PATH to query fasta: ").strip('"')

        # Set output directory if default value is used
        if args.out == './':
            args.out = str(args.query_fasta).replace('.fasta', '')
            if not os.path.isdir(args.out):
                os.mkdir(Path(args.out))  # Create the output directory

    ## CHECK IF FILES ALREADY EXIST
    project_folder = Path(args.out)  # Use the output directory specified by the user

    # Convert db to Path
    database = Path(args.database.strip('"'))
    db_name = Path(database).stem

    # Verify db
    if ' ' in str(database):
        print('\nError: Your database PATH contains white spaces! Please move your database to a different folder!')
        print('')
        return
    db_path = database / 'db'
    db_check = check_blast_db(db_path)
    if db_check == False:
        print('')
        return

    if args.query_fasta:
        # Run the BLASTn function
        a_blastn(args.blastn_exe,
                 args.query_fasta.strip('"'),
                 database,
                 project_folder,
                 args.n_cores,
                 args.task,
                 args.subset_size,
                 args.max_target_seqs,
                 args.masking,
                                  )
    else:
        print('\nError: Please provide a fasta file!')

    # Handle the 'filter' command
    if not os.path.isfile(Path(project_folder).joinpath('log.txt')):
        print('\nError: Could not find the BLAST results folder!')
    else:
        print('')
        # Run the filter function
        if not args.filter:
            filter_mode = "2"
        else:
            filter_mode = str(args.filter)
        if filter_mode == "1":
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Using filter mode 1: similarity (-{args.sim_range}%) -> e-value -> rating (optional)')
        elif filter_mode == "2":
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Using filter mode 2: e-value -> similarity (-{args.sim_range}%) -> rating (optional)')
        else:
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Using filter mode 3: similarity (-{args.sim_range}%) -> rating (optional)')
        b_filter(project_folder, database, args.thresholds, args.n_cores, filter_mode, args.rating_range, args.sim_range)

    # Perform reblast if defined
    if args.reblast_db:
        print(f'\n{datetime.datetime.now().strftime("%H:%M:%S")}: Performing a re-blast of low-smiliarity sequences.')

        if not args.reblast_db:
            print(f'\n{datetime.datetime.now().strftime("%H:%M:%S")}: Please provide a secondary database.')
            sys.exit(1)

        # Define required folders and files
        reblast_db = Path(args.reblast_db)
        reblast_folder = project_folder / 'reblast'
        os.makedirs(reblast_folder, exist_ok=True)
        query_fasta = Path(args.query_fasta)
        low_similarity_query_fasta = reblast_folder / 'lowsim.fasta'
        name = project_folder.name
        blastn_filtered_xlsx = Path('{}/{}_taxonomy.xlsx'.format(project_folder, name))

        # Extract all sequences with less than 97% similarity (default) and write them to a new fasta
        reblast_sim = args.reblast_sim
        taxonomy_df = pd.read_excel(blastn_filtered_xlsx).fillna('')
        lowsim_df = taxonomy_df[(taxonomy_df['Similarity'] < reblast_sim) |(taxonomy_df['Species'] == '')]
        lowsim_hashes = lowsim_df['unique ID'].values.tolist()
        n_lowsim = len(lowsim_hashes)
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Found {n_lowsim} sequences for re-assignment (<{reblast_sim}% sim. or non-species-level).')

        with open(low_similarity_query_fasta, "w") as out_handle:
            for record in SeqIO.parse(query_fasta, "fasta"):
                if record.id in lowsim_hashes:
                    # Write to output FASTA
                    SeqIO.write(record, out_handle, "fasta")

        if reblast_db == Path('boldigger'):
            pass
            #run_boldigger3(low_similarity_query_fasta)
            #boldigger_results_xlsx = reblast_folder / 'boldigger3_data' / 'lowsim_identification_result.xlsx'
            #boldigger_df = pd.read_excel(boldigger_results_xlsx).fillna('')
            #merged_df = merge_boldigger3(taxonomy_df, boldigger_df, 'boldigger3')
            #col = "Status"  # column to move
            #merged_df = merged_df[[c for c in merged_df.columns if c != col] + [col]]
            #blastn_reblast_filtered_xlsx = Path('{}/{}_reblast_taxonomy.xlsx'.format(project_folder, name))
            #merged_df.to_excel(blastn_reblast_filtered_xlsx, index=False)
            #print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Finished re-assignment using BOLDigger3!')

        else:
            # Run the BLASTn function
            a_blastn(args.blastn_exe,
                     low_similarity_query_fasta,
                     reblast_db,
                     reblast_folder,
                     args.n_cores,
                     args.task,
                     args.subset_size,
                     args.max_target_seqs,
                     args.masking,
                                      )
            # Run the filter command
            print('')
            b_filter(reblast_folder, reblast_db, args.thresholds, args.n_cores, filter_mode, args.rating_range, args.sim_range)

            taxonomy_reblast_xlsx = reblast_folder / 'reblast_taxonomy.xlsx'
            taxonomy_reblast_df = pd.read_excel(taxonomy_reblast_xlsx).fillna('')
            merged_df = merge_reblast(taxonomy_df, taxonomy_reblast_df, db_name)
            col = "Status"  # column to move
            merged_df = merged_df[[c for c in merged_df.columns if c != col] + [col]]
            blastn_reblast_filtered_xlsx = Path('{}/{}_reblast_taxonomy.xlsx'.format(project_folder, name))
            merged_df.to_excel(blastn_reblast_filtered_xlsx, index=False)
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Finished re-blast!')

    # If defined, use an issue list
    issue_list_xlsx = False
    if issue_list_xlsx == True:
        issue_list_flag(taxonomy_df, issue_list_df)

    # save settings to df
    df = settings_to_df(args)
    settings_xlsx = Path('{}/{}_blastn_run_settings.xlsx'.format(project_folder, db_name))
    df.to_excel(settings_xlsx, index=False)

    ## Finish script
    print(f'\n{datetime.datetime.now().strftime("%H:%M:%S")}: Done! Have a nice day!')

# Run the main function if script is called directly
if __name__ == "__main__":
    main()