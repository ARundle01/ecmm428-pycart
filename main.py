import geojson
import shapely.ops
from shapely import geometry

import cartogram

# import os
# os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from libpysal.weights import W


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
    cua = "./data/Dec2020/Counties_and_Unitary_Authorities_(December_2020)_UK_BGC.geojson"
    countries = "./data/Dec2020/Countries_(December_2020)_UK_BGC.geojson"
    lad = "./data/Dec2020/Local_Authority_Districts_(December_2020)_UK_BGC.geojson"
    regions = "./data/Dec2020/Regions_(December_2020)_EN_BGC.geojson"
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

        temp_lad = "./data/Dec2021/Local_Authority_Districts_(December_2021)_GB_BGC.geojson"
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


if __name__ == '__main__':
    pop = "./data/Dec2020/ukpopestimatesdec2020.csv"

    cmap = plt.get_cmap('Reds')
    new_cmap = truncate_colormap(cmap, 0.2, 0.9)

    fig, ax = plt.subplots(1, figsize=(4, 4))
    ax.axis('equal')
    ax.axis('off')

    pop_df = parse_pop(pop)

    # 'countries' | 'regions' | 'cua' | 'lad'
    try:
        places_df, features, code_type = init_geojson("lad")
    except NameError as e:
        print(e)
        exit(-1)

    gdf = make_gdf(places_df, pop_df)

    # gdf.plot(color='w', ax=ax, zorder=0, edgecolor='0', linewidth=0.1, legend=False)

    cart = cartogram.Cartogram(gdf, "Population", id_field="Name")
    # dorling = cart.dorling(iterations=100, stop=None)
    # non_con = cart.non_contiguous(position='centroid', size_value=1.0)
    fftrho = cart.diffusion()

    # gdf.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)

    # diffusion.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)

    # places_df.plot(ax=ax, edgecolor='0', linewidth=0.1)

    # gdf.plot(column='Population', cmap=new_cmap, ax=ax, edgecolor='0', linewidth=0.1, legend=True)
    # ax.set_title('Population of LADs', fontdict={'fontsize': '15', 'fontweight': '3'})

    # Plot Figure
    # plt.savefig("./out/fishnet_256_density_2.png", dpi=1200)


