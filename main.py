import geojson

from pycart import cartogram

import geopandas as gpd

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from libpysal.weights import W
from dataprep.clean import clean_country


def parse_geojson(fname, is_pop=False):
    with open(fname) as f:
        gj = geojson.load(f)
    features = gj['features']

    keys = list(features[0]['properties'].keys())

    if is_pop:
        code_type = keys[1]
        name_type = keys[2]

        return features, code_type, name_type
    else:
        return features, keys


def init_geojson(map_type):
    cua = "./data/UK/Dec2020/Counties_and_Unitary_Authorities_(December_2020)_UK_BGC.geojson"
    countries = "./data/UK/Dec2020/Countries_(December_2020)_UK_BGC.geojson"
    lad = "./data/UK/Dec2020/Local_Authority_Districts_(December_2020)_UK_BGC.geojson"
    regions = "./data/UK/Dec2020/Regions_(December_2020)_EN_BGC.geojson"
    wales_name = None
    removal = None

    if map_type.lower() == "countries":
        features, code_type, name_type = parse_geojson(countries, True)
    elif map_type.lower() == "regions":
        features, code_type, name_type = parse_geojson(regions, True)
    elif map_type.lower() == "cua":
        features, code_type, name_type = parse_geojson(cua, True)
        wales_name = "CTYUA20NMW"
        removal = ["E10000021"]
    elif map_type.lower() == "lad":
        features, code_type, name_type = parse_geojson(lad, True)
        wales_name = "LAD20NMW"
        removal = ["E07000150", "E07000151", "E07000152", "E07000153", "E07000154", "E07000155", "E07000156"]
    else:
        raise NameError("Map not supported, use: 'countries' | 'regions' | 'cua' | 'lad'")

    features_df = gpd.GeoDataFrame.from_features(features)

    if wales_name is not None and removal is not None:
        to_remove = features_df[features_df[code_type].isin(removal)].index
        features_df.drop(['OBJECTID', name_type, wales_name], axis=1, inplace=True)
        features_df.drop(to_remove, axis=0, inplace=True)

        temp_lad = "./data/UK/Dec2021/Local_Authority_Districts_(December_2021)_GB_BGC.geojson"
        temp_features, temp_code_type, temp_name_type = parse_geojson(temp_lad, True)
        temp_features_df = gpd.GeoDataFrame.from_features(temp_features)
        temp_features_df.drop(['OBJECTID', 'LAD21NMW', 'GlobalID', temp_name_type], axis=1, inplace=True)
        temp_features_df.rename(columns={temp_code_type: code_type}, inplace=True)

        addition = ["E06000061", "E06000062"]
        to_add = temp_features_df[temp_features_df[code_type].isin(addition)]
        features_df = gpd.GeoDataFrame(pd.concat([features_df, to_add], ignore_index=True))
    else:
        features_df.drop(['OBJECTID', name_type], axis=1, inplace=True)

    if map_type.lower() != "lad":
        features_df.drop(['GlobalID'], axis=1, inplace=True)

    return features_df, features, code_type


def to_int(x):
    int_x = None
    try:
        int_x = int(x.replace(',', ''))
    except AttributeError:
        int_x = int(x)
    finally:
        return int_x


def parse_pop(fname):
    pop_df = pd.read_csv(fname)
    pop_df['Population'] = pop_df['Population'].apply(to_int)

    return pop_df


def get_sub_pop(pop, features_df, code_type):
    sub_pop = pd.DataFrame()
    pop.rename(columns={'Code': code_type}, inplace=True)

    features = features_df[code_type].to_list()
    for feature in features:
        sub_pop = pd.concat([sub_pop, pop[pop[code_type] == feature]])

    sub_pop.reset_index(drop=True, inplace=True)
    return sub_pop


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n}, {a:.2f}, {b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))

    return new_cmap


def make_gdf(places_df, pop_df):
    sub_pop = get_sub_pop(pop_df, places_df, code_type)
    geo_pop = places_df.merge(sub_pop, on=code_type)

    gdf = gpd.GeoDataFrame(geo_pop).set_crs('EPSG:3857')

    return gdf


def parse_world_pop(fname):
    world_df = pd.read_csv(fname)

    world_df['Country'] = world_df[['Country Name']]
    world_df['ISO'] = world_df[['Country Code']]
    world_df['Population'] = world_df[['2021']]


    world_df = world_df[['Country', 'ISO', 'Population']]

    world_df = world_df.astype({'Country': str, 'ISO': str, 'Population': float})

    return world_df


def init_world_geojson():
    pop = "./data/World/API_SP.POP.TOTL_DS2_en_csv_v2_5358404.csv"
    geo = "./data/World/World_Countries_(Generalized).geojson"

    world_pop = parse_world_pop(pop)
    world_f, world_k = parse_geojson(geo)

    features_df = gpd.GeoDataFrame.from_features(world_f)

    features_df['Country'] = features_df['COUNTRY']
    features_df = features_df[['Country', 'ISO', 'SHAPE_Length', 'SHAPE_Area', 'geometry']]

    # # Drop all irrelevant rows
    # to_drop = [
    #     'WLD',
    #     'IBT',
    #     'LMY',
    #     'MIC',
    #     'IBD',
    #     'EAR',
    #     'LMC',
    #     'UMC',
    #     'EAS',
    #     'LTE',
    #     'EAP',
    #     'TEA',
    #     'SAS',
    #     'TSA',
    #     'IDA',
    #     'OED',
    #     'HIC',
    #     'IDX',
    #     'SSF',
    #     'TSS',
    #     'SSA',
    #     'PST',
    #     'LDC',
    #     'PRE',
    #     'FCS',
    #     'ECS',
    #     'HPC',
    #     'LIC',
    #     'AFE',
    #     'LCN',
    #     'TLA',
    #     'IDB',
    #     'LAC',
    #     'MEA',
    #     'AFW',
    #     'TEC',
    #     'ARB',
    #     'EUU',
    #     'MNA',
    #     'TMN',
    #     'ECA',
    #     'NAC',
    #     'EMU',
    #     'CEB',
    #     'SST',
    #     'OSS',
    #     'CSS',
    #     'PSS',
    #     'INX',
    # ]
    #
    # dropping = world_pop['ISO'].isin(to_drop)
    # world_pop = world_pop[~dropping].reset_index(drop=True)

    world_pop = clean_country(world_pop, "ISO", input_format="alpha-3", output_format="alpha-2")
    world_pop.drop(columns=["ISO"], inplace=True)
    world_pop.rename(columns={"ISO_clean": "ISO"}, inplace=True)

    return world_pop, features_df


def combine_geo_val(geo, val):
    combined = val.merge(geo, on='ISO')
    combined = combined.drop(columns={'Country_y'})
    combined = combined.rename(columns={'Country_x': 'Country'})

    combined = gpd.GeoDataFrame(combined).set_crs('EPSG:3857')

    return combined


if __name__ == '__main__':
    pop = "./data/UK/Dec2020/ukpopestimatesdec2020.csv"
    world_pop = "./data/World/API_SP.POP.TOTL_DS2_en_csv_v2_5358404.csv"
    world_geo = "./data/World/World_Countries_(Generalized).geojson"

    # world_df = parse_world_pop(world_pop)
    #
    # world_f, world_k = parse_geojson(world_geo)
    #
    # world_f_df = gpd.GeoDataFrame.from_features(world_f)

    world_df, world_features = init_world_geojson()

    combined = combine_geo_val(world_features, world_df)

    cmap = plt.get_cmap('Reds')
    new_cmap = truncate_colormap(cmap, 0.2, 0.9)

    fig, ax = plt.subplots(1)
    ax.axis('equal')
    ax.axis('off')

    pop_df = parse_pop(pop)

    # 'countries' | 'regions' | 'cua' | 'lad'
    try:
        places_df, features, code_type = init_geojson("cua")
    except NameError as e:
        print(e)
        exit(-1)

    gdf = make_gdf(places_df, pop_df)

    # gdf.plot(color='w', ax=ax, zorder=0, edgecolor='0', linewidth=0.1, legend=False)

    # cart = cartogram.Cartogram(gdf, value_field="Population", id_field="Name", geometry_field="geometry")
    # dorling = cart.dorling(iterations=100, stop=None)
    # non_con = cart.non_contiguous(position='centroid', size_value=1.0)

    # gdf.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    # non_con.plot(color='r', ax=ax, edgecolor='0', linewidth=0.1, legend=False)
    # dorling.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)

    combined_cart = cartogram.Cartogram(combined, "Population", id_field="ISO", geometry_field="geometry")
    combined_dorling = combined_cart.dorling(iterations=5, stop=None)
    # combined_non_con = combined_cart.non_contiguous(size_value=5.0)

    combined.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    # combined_non_con.plot(color='r', ax=ax, edgecolor='0', linewidth=0.1, legend=False)
    combined_dorling.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)

    # Plot Figure
    plt.savefig("./out/dorling_world_new5.png", dpi=1200)
