from pygbif import occurrences as occ
from geopy.distance import geodesic
import math

species = "Lutra lutra"
lat = 49.7722
lon = 6.6492
radius_km = 50
after_year = 2020

def gbif_validation(species, lat, lon, radius_km, after_year):

    # Bounding box
    lat_delta = radius_km / 111.32
    lon_delta = radius_km / (111.32 * math.cos(math.radians(lat)))

    result = occ.search(
        scientificName=species,
        decimalLatitude=f"{lat-lat_delta},{lat+lat_delta}",
        decimalLongitude=f"{lon-lon_delta},{lon+lon_delta}",
        limit=300
    )

    # Keep only records within 50 km
    nearby = []

    for record in result["results"]:
        try:
            if "decimalLatitude" not in record or "decimalLongitude" not in record:
                continue

            if species != record['species']:
                continue

            distance = geodesic(
                (lat, lon),
                (record["decimalLatitude"], record["decimalLongitude"])
                ).km

            if distance <= radius_km:
                nearby.append(record)
        except:
            pass

    print(f"{species} (max=300)\n"
            f"  - after {after_year} = {len(nearby)}"
          )
    return [species, len(nearby)]

import pandas as pd
df = pd.read_excel('/Users/tillmacher/Desktop/TTT_projects/Lippe_eRNA_tele02_vertebrates/TaXon_tables/Lippe_eRNA_tele02_taxon_table_cons_NCsub_vertebrates.xlsx')
species_lst = df['Species'].fillna('').drop_duplicates().values.tolist()
for species in species_lst:
    gbif_validation(species, lat, lon, radius_km, after_year)
