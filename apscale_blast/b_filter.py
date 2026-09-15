import argparse
import os
import sys
from pathlib import Path
import datetime
import subprocess
import shutil
import time
import glob
import multiprocessing
import pandas as pd
import pandas as pdE
import numpy as np
from joblib import Parallel, delayed
import threading
import pyarrow as pa
import pyarrow.parquet as pq
from ete3 import NCBITaxa, Tree
import requests
import xmltodict
from tqdm import tqdm

# Lock for thread-safe print statements
print_lock = threading.Lock()

# Ensure multiprocessing works when frozen (e.g. for executables)
multiprocessing.freeze_support()

def accession2taxonomy(df_1, taxid_df, col_names_2, db_name):
    """
    Convert accession IDs to taxonomy information.
    """

    # file = '/Volumes/Coruscant/APSCALE_projects/APSCALE_databases_input/blast_test/BOLD/subsets/subset_1_megablast.csv'
    # df_1 = pd.read_csv(file, header=None, sep=';;', names=col_names, engine='python').fillna('NAN')
    # taxid_df = pd.read_parquet('/Volumes/Coruscant/APSCALE_projects/APSCALE_databases/db_BOLD_Public_30-Jun-2026/db_taxonomy.parquet.snappy')

    df_2_list = []
    cols = taxid_df.columns.tolist()
    levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    for row in df_1.values.tolist():
        ID_name = row[0]
        accession = row[1]
        evalue = row[-1]
        similarity = row[-2]
        try:
            taxonomy = taxid_df[taxid_df['Accession'] == accession][levels].values.tolist()[0]
        except:
            try:
                taxonomy = taxid_df[taxid_df['Accession'].str.startswith(accession[:30])][levels].values.tolist()[0]
            except:
                taxonomy = ['Missing-Accession'] * 7

        if 'rating' in cols:
            rating = taxid_df[taxid_df['Accession'] == accession]['rating'].values.tolist()
            df_2_list.append([ID_name] + taxonomy + [similarity, evalue] + rating + [accession])
        else:
            df_2_list.append([ID_name] + taxonomy + [similarity, evalue, 0, accession])

    df_2 = pd.DataFrame(df_2_list, columns=col_names_2 + ['rating', 'accession'])
    return df_2

def filter_blastn_csvs(file, taxid_df, i, n_subsets, thresholds, db_name, filter_mode, rating_range, sim_range, n_jobs):
    """
    Filter BLASTn results from CSV files.
    """

    blastn_filtered_parquet = Path(file.replace('.csv', '_filtered.parquet.snappy'))
    if blastn_filtered_parquet.exists():
        print('{}: Skipping {} (already exists and is not empty).'.format(datetime.datetime.now().strftime('%H:%M:%S'),
                                                                          blastn_filtered_parquet.stem))
        return

    ## load blast results
    col_names = ['unique ID', 'Sequence ID', 'Similarity', 'evalue']

    if not os.path.isfile(file):
        print('{}: Skipping missing subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))
        return
    if os.path.getsize(file) == 0:
        print('{}: Skipping empty subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))
        return

    csv_df = pd.read_csv(file, header=None, sep=';;', names=col_names, engine='python').fillna('NAN')
    csv_df['Similarity'] = [float(i) for i in csv_df['Similarity'].values.tolist()]

    if len(csv_df) == 0:
        print('{}: Error during filtering for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))
        return

    thresholds = thresholds.split(',')
    # Ensure correct number of thresholds
    if len(thresholds) != 5:
        print('Please provide 5 comma-separated threshold values!')
        print('Using default values...')
        thresholds = ['97', '95', '90', '87', '85']

    species_threshold = int(thresholds[0])
    genus_threshold = int(thresholds[1])
    family_threshold = int(thresholds[2])
    order_threshold = int(thresholds[3])
    class_threshold = int(thresholds[4])

    ## filter hits
    ID_set = csv_df['unique ID'].drop_duplicates().values.tolist()
    col_names_2 = ['unique ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'evalue']
    taxonomy_df = pd.DataFrame()

    # loop through IDs
    # ID = "OTU_39"
    # for ID in (ID_set, desc=f'Subset {i + 1}/{n_subsets}', position=i % n_jobs, leave=True):
    for ID in ID_set:
        ## filter by evalue
        df_0 = csv_df.loc[csv_df['unique ID'] == ID]

        ##### 1) FILTERING BY SIMILARITY THEN BY EVALUE (or the other way round, let the user decide)
        filter_mode = int(filter_mode)
        sim_range = int(sim_range)
        if filter_mode == 1:
            max_sim = max(df_0['Similarity'])
            sim_cut = max_sim - sim_range
            df_1 = df_0.loc[df_0['Similarity'] >= sim_cut]
            min_e = min(df_1['evalue'])
            df_1 = df_1.loc[df_1['evalue'] == min_e]
        elif filter_mode == 2:
            min_e = min(df_0['evalue'])
            df_1 = df_0.loc[df_0['evalue'] == min_e]
            max_sim = max(df_1['Similarity'])
            sim_cut = max_sim - sim_range
            df_1 = df_1.loc[df_1['Similarity'] >= sim_cut]
        else:
            max_sim = max(df_0['Similarity'])
            sim_cut = max_sim - sim_range
            df_1 = df_0.loc[df_0['Similarity'] >= sim_cut]
            min_e = min(df_1['evalue'])

        ############################################################################################################
        ## convert fasta headers to taxonomy
        df_2 = accession2taxonomy(df_1, taxid_df, col_names_2, db_name)
        del df_0
        del df_1

        ## Filter out missing taxids
        if df_2['Species'].str.contains(r'\bUnknown taxid\b', na=False).any():
            df_2_reduced = df_2.loc[~df_2['Species'].str.contains(r'\bUnknown taxid\b', na=False)]
            if len(df_2_reduced) != 0:
                df_2 = df_2_reduced.copy()

        ############################################################################################################
        ## Filter by rating
        min_rating = int(max(df_2['rating'])) - int(rating_range)
        rating_dropped = [', '.join([f'{i} (r>{min_rating})' for i in df_2.loc[df_2['rating'] < min_rating]['Species'].unique()])]
        df_2 = df_2.loc[df_2['rating'] >= min_rating]
        dup_cols = [c for c in df_2.columns if c not in ('rating', 'accession')]
        group_sizes = df_2.groupby(dup_cols).size().rename("size")
        idx = df_2.groupby(dup_cols)['rating'].idxmax()
        df_2 = df_2.loc[idx].reset_index(drop=True)
        df_2 = df_2.merge(group_sizes, on=dup_cols, how="left")
        df_2 = df_2.drop(columns=['accession'])

        ############################################################################################################
        ##### 2) ROBUSTNESS TRIMMING BY SIMILARITY

        if max_sim >= species_threshold:
            pass
        elif max_sim < species_threshold and max_sim >= genus_threshold:
            df_2['Species'] = ''
        elif max_sim < genus_threshold and max_sim >= family_threshold:
            df_2['Species'] = ''
            df_2['Genus'] = ''
        elif max_sim < family_threshold and max_sim >= order_threshold:
            df_2['Species'] = ''
            df_2['Genus'] = ''
            df_2['Family'] = ''
        elif max_sim < order_threshold and max_sim >= class_threshold:
            df_2['Species'] = ''
            df_2['Genus'] = ''
            df_2['Family'] = ''
            df_2['Order'] = ''
        else:
            df_2['Species'] = ''
            df_2['Genus'] = ''
            df_2['Family'] = ''
            df_2['Order'] = ''
            df_2['Class'] = ''

        ##### 3) REMOVAL OF DUPLICATE HITS
        df_3 = df_2.drop_duplicates().copy()

        ##### 4) MORE THAN ONE TAXON REMAINING?
        if len(df_3) == 1:
            df_3['Flag'] = ''
            df_3['Ambiguous taxa'] = ''
            df_3 = df_3.drop(columns=['size'], errors='ignore')
            df_3['Rating drop'] = rating_dropped * len(df_3) if rating_dropped else ''
            taxonomy_df = pd.concat([taxonomy_df, df_3])

        ##### 5) SPECIES LEVEL REFERENCE?
        elif max_sim < species_threshold:
            # remove taxonomic levels until a single hit remains
            for level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'][::-1]:
                n_hits = len(df_3)
                df_3.loc[:, level] = [''] * n_hits
                df_3 = df_3.drop_duplicates()
                if len(df_3) == 1:
                    break
            df_3['Flag'] = ''
            df_3['Ambiguous taxa'] = ''
            df_3['Rating drop'] = rating_dropped * len(df_3) if rating_dropped else ''
            df_3 = df_3.drop(columns=['size'], errors='ignore')
            taxonomy_df = pd.concat([taxonomy_df, df_3])

        ##### 6) AMBIGOUS SPECIES WORKFLOW (FLAGGING SYSTEM)
        else:
            ##### 7) DOMINANT SPECIES PRESENT? (F1)
            df_2_dominant = df_2.loc[df_2['size'] == max(df_2['size'])].drop_duplicates()

            if len(df_2_dominant) == 1:
                df_3 = df_2_dominant.copy()
                df_3['Flag'] = ['F1 (Dominant species)'] * len(df_3)
                df_3["Ambiguous taxa"] = [", ".join(sorted([f"{s} (r={r}, n={n})" for s, r, n in
                                                            df_2.drop_duplicates(subset=["Species", "rating", "size"])[
                                                                ["Species", "rating", "size"]].values]))] * len(df_3)
                df_3['Rating drop'] = rating_dropped * len(df_3) if rating_dropped else ''
                df_3 = df_3.drop(columns=['size'], errors='ignore')
                taxonomy_df = pd.concat([taxonomy_df, df_3])
            else:
                ##### 8) TWO SPECIES OF ONE GENUS? (F2)
                df_3 = df_3.drop(columns=['size'], errors='ignore')
                n_genera = len(set(df_3['Genus']))
                if n_genera == 1 and len(df_3) == 2:
                    genus = df_3['Genus'].values.tolist()[0]
                    species = '/'.join([i.replace(genus + ' ', '') for i in sorted(df_3['Species'].values.tolist())])
                    df_3['Species'] = ['{} {}'.format(genus, species)]*len(df_3)
                    df_3['Flag'] = ['F2 (Two species of one genus)'] * len(df_3)
                    df_3['Ambiguous taxa'] = [', '.join(sorted([f"{s} (r={r})" for s, r in df_2[['Species', 'rating']].drop_duplicates().values]))] * len(df_3)
                    df_3['Rating drop'] = rating_dropped * len(df_3) if rating_dropped else ''
                    taxonomy_df = pd.concat([taxonomy_df, df_3])

                ##### 9) MULTIPLE SPECIES OF ONE GENUS? (F3)
                elif n_genera == 1 and len(df_3) != 2:
                    genus = df_3['Genus'].values.tolist()[0]
                    df_3['Species'] = [genus + ' sp.']*len(df_3)
                    df_3['Flag'] = ['F3 (Multiple species of one genus)'] * len(df_3)
                    df_3['Ambiguous taxa'] = [', '.join(sorted([f"{s} (r={r})" for s, r in df_2[['Species', 'rating']].drop_duplicates().values]))] * len(df_3)
                    df_3['Rating drop'] = rating_dropped * len(df_3) if rating_dropped else ''
                    taxonomy_df = pd.concat([taxonomy_df, df_3])

                ##### 10) TRIMMING TO MOST RECENT COMMON TAXON (F4)
                # df_3 = df_2.drop_duplicates().copy()
                else:
                    # remove taxonomic levels until a single hit remains
                    for level in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'][::-1]:
                        n_hits = len(df_3)
                        if df_3[level].nunique() == 1:
                            break
                        df_3.loc[:, level] = ['']*n_hits
                    df_3 = df_3.copy().head(1)
                    df_3['Flag'] = ['F4 (Trimming to MRCA)'] * len(df_3)
                    df_3['Ambiguous taxa'] = [', '.join(sorted([f"{s} (r={r})" for s, r in df_2[['Species', 'rating']].drop_duplicates().values]))] * len(df_3)
                    df_3['Rating drop'] = rating_dropped * len(df_3) if rating_dropped else ''
                    taxonomy_df = pd.concat([taxonomy_df, df_3])

    # export dataframe
    taxonomy_df.columns = ['unique ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'E-value', 'Rating', 'Flag', 'Ambiguous taxa', 'Rating drop']

    # Fix type problems before parquet export
    taxonomy_df['Rating'] = pd.to_numeric(taxonomy_df['Rating'], errors='coerce').astype('Int64')

    # export
    taxonomy_df.to_parquet(blastn_filtered_parquet, compression='snappy')

    # final print
    print('{}: Finished filtering of subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))

def main(blastn_folder, blastn_database, thresholds, n_cores, filter_mode, rating_range, sim_range):
    """ Filter results according to Macher et al., 2023 (Fish Mock Community paper) """

    print('{}: Starting to filter blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))
    print('{}: Your database: {}'.format(datetime.datetime.now().strftime('%H:%M:%S'),  Path(blastn_database).stem))

    ## load blast results
    print('{}: Importing Accession IDs'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    csv_files = glob.glob('{}/subsets/*.csv'.format(blastn_folder))
    all_ids = []
    for file in csv_files:
        csv_df = pd.read_csv(file, header=None, sep=';;', engine='python').fillna('NAN')
        ids = csv_df[1].unique().tolist()
        all_ids.extend(ids)
    all_ids = list(set(all_ids))

    ## load taxid table
    print('{}: Importing reference taxonomy'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    taxid_table = Path(blastn_database).joinpath('db_taxonomy.parquet.snappy')
    taxid_df = pd.read_parquet(taxid_table).fillna('')
    # filter taxid_df for present ids
    taxid_df = taxid_df.loc[taxid_df['Accession'].isin(all_ids)].reset_index(drop=True)

    ## collect name of database
    db_name = Path(blastn_database).stem

    # PARALLEL FILTER COMMAND
    n_subsets = len(csv_files)
    # file = csv_files[2]
    # i = 0
    print('{}: Filtering and flagging of {} files...'.format(datetime.datetime.now().strftime('%H:%M:%S'),  len(csv_files)))
    Parallel(n_jobs=n_cores, backend='loky')(delayed(filter_blastn_csvs)(file, taxid_df, i, n_subsets, thresholds, db_name, filter_mode, rating_range, sim_range, n_cores) for i, file in enumerate(csv_files))

    ## also already define the no match row
    NoMatch = ['NoMatch'] * 7 + [0, 1, '', '']

    # Get a list of all the xlsx files
    parquet_files = glob.glob('{}/subsets/*_filtered.parquet.snappy'.format(blastn_folder))

    if len(parquet_files) == 0:
        print('## Warning ##')
        print('{}: Error during filtering of the blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))
        print('{}: Please check your database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        print('{}: Errors usually occur during unzipping of the database file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        print('## Warning ##')
    else:

        # Create a list to hold all the individual DataFrames
        df_list = []

        # Loop through the list of xlsx files
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Collecting {len(parquet_files)} filtered subsets.')
        for file in tqdm(parquet_files):
            # Read each xlsx file into a DataFrame
            df = pd.read_parquet(file).fillna('')
            # Append the DataFrame to the list
            df_list.append(df)

        # Concatenate all the DataFrames in the list into one DataFrame
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Concatenating subsets.')
        merged_df = pd.concat(df_list, ignore_index=True)
        name = Path(blastn_folder).name
        blastn_filtered_xlsx = Path('{}/{}_taxonomy.xlsx'.format(blastn_folder, name))
        blastn_filtered_parquet = Path('{}/{}_taxonomy.parquet.snappy'.format(blastn_folder, name))

        ## add OTUs without hit
        # Drop duplicates in the DataFrame
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Dropping duplicates.')
        merged_df = merged_df.drop_duplicates()
        output_df_list = []

        # Read the IDs from the file
        ID_list = Path(blastn_folder).joinpath('IDs.txt')
        IDs = [i.rstrip() for i in ID_list.open()]

        # Check if each ID is already in the DataFrame
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Constructing taxonomy table.')
        for ID in tqdm(IDs):
            if ID not in merged_df['unique ID'].values.tolist():
                # Create a new row with the ID and other relevant information
                row = [ID] + NoMatch
                output_df_list.append(row)
            else:
                row = merged_df.loc[merged_df['unique ID'] == ID].values.tolist()[0]
                output_df_list.append(row)

        ## sort table
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Saving taxonomy table...')

        merged_df.columns.tolist()
        output_df = pd.DataFrame(output_df_list, columns=merged_df.columns.tolist())
        output_df['Status'] = db_name
        try:
            output_df.to_excel(blastn_filtered_xlsx, sheet_name='Taxonomy table', index=False)
        except:
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Unable to write to Excel.')
        try:
            output_df.to_parquet(blastn_filtered_parquet, compression='snappy')
        except:
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Unable to write to Parquet.')

        print('{}: Finished to filter blast results for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), blastn_folder))

if __name__ == '__main__':
    main()
