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

    # squares_features, squares_keys = parse_geojson("./data/Test/square_test.geojson")
    # squares_df = gpd.GeoDataFrame.from_features(squares_features)
    #
    # squares_pop = parse_pop("./data/Test/squares_pop.csv")
    #
    # squares_geo_pop = squares_df.merge(squares_pop, on="name")
    # squares_geo_pop.drop(["shape"], axis=1, inplace=True)

    # squares_cart = cartogram.Cartogram(squares_geo_pop, "Population", id_field="name")
    # squares_noncon = squares_cart.non_contiguous(position="centroid", size_value=1.0)
    # squares_dorling = squares_cart.dorling(iterations=256)

    # squares_geo_pop.plot(color='w', ax=ax, zorder=0, edgecolor='0', linewidth=0.5, legend=False)
    # squares_noncon.plot(color='r', ax=ax, edgecolor='0', linewidth=0.1, legend=False)
    # squares_dorling.plot(color='blue', ax=ax, edgecolor='0', linewidth=0.1, legend=False)
    # plt.savefig("squares_dorling.png", dpi=750)

    # gdf.plot(color='w', ax=ax, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    # gdf.to_file('out.shp')

    # londonless_geo_pop = geo_pop[geo_pop["Name"] != "LONDON"].reset_index(drop=True)

    # cart = cartogram.Cartogram(geo_pop, "Population", id_field="Name")
    # borders, islands, regions, current_region, neighbours, xrepel, yrepel, xattract, yattract, repel_dist, attract_dist = cart.dorling(iterations=1, stop=77)
    # regions = cart.dorling(iterations=100, stop=None)
    # diffusion = cart.diffusion()

    # for i in range(1, 101):
    #     fig, ax = plt.subplots(1, figsize=(4, 4))
    #     ax.axis('equal')
    #     ax.axis('off')
    #     plt.xlim((-9.492274504549895, 2.6818561740215383))
    #     plt.ylim((49.31499355, 61.41058145))
    #
    #     print(f"Running {i} iterations:\n")
    #     regions, displacement = cart.dorling(iterations=i)
    #     geo_pop.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    #     regions.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    #
    #     plt.savefig(f"./dorling_out/dorling_iters_{i}.png", dpi=1200)
    #
    #     print("\n")
    #     plt.cla()

    # centroids = regions["geometry"].to_numpy()
    # centroids = np.vstack(centroids)
    # plt.plot(centroids[:, 0], centroids[:, 1], '.', color='r', markersize=0.2)

    # borders_W = W.from_adjlist(borders)
    #
    # for k, neighs in borders_W.neighbors.items():
    #     if k not in islands:
    #         origin = centroids[k]
    #         for neigh in neighs:
    #             segment = centroids[[k, neigh]]
    #             plt.plot(segment[:, 0], segment[:, 1], '-', linewidth=0.2)

    # no_island_centroids = regions.drop(islands, axis=0, inplace=False)
    # no_island_centroids = np.vstack(no_island_centroids["geometry"].to_numpy())
    # plt.plot(no_island_centroids[:,0], no_island_centroids[:,1], '.', color='r', markersize=0.2)

    # if islands is not None:
    #     regions = regions.drop(islands, axis=0, inplace=False)
    #
    # radius_regions = gpd.GeoDataFrame(
    #     data=regions.drop(columns=["geometry", "radius"]),
    #     geometry=regions.apply(lambda x: x["geometry"].buffer(x["radius"]), axis=1)
    # )
    #
    # neighbours_idx = neighbours.index.values.tolist()

    # non_con = cart.non_contiguous(position='centroid', size_value=1.0)
    # dorling = cart.old_dorling(iterations=1000)

    # borders = cartogram.shared_borders(geo_pop)
    # knn_borders = cartogram.knn_borders(geo_pop)

    # wq, wk, WK = borders.get_borders(geo_pop)
    # neighbours, weights = borders.get_borders(geo_pop)

    # for neighbour, weight in zip(neighbours, weights):
    #     if len(neighbours[neighbour]) != len(weights[weight]):
    #         print(f"Neighbour {neighbour} and Weight {weight} are mismatched")

    # WQ, WK = borders.get_borders(geo_pop)
    # W = borders.get_borders(geo_pop)

    # test_geo_pop = geo_pop.copy()
    #
    # for island in wq.islands:
    #     test_geo_pop.drop([island], inplace=True)
    #
    # test_geo_pop.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)

    # if 0 in knn_borders["focal"].values:
    #     if 3 in knn_borders[knn_borders["focal"] == 0]["neighbor"].values:
    #         print("it is in here")

    # unique_borders = knn_borders.merge(
    #     borders,
    #     how='left',
    #     indicator=True
    # )

    # geo_pop.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)
    minx, miny, maxx, maxy = gdf.total_bounds
    gdf['density'] = gdf['Population'] / gdf.area
    # gdf.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)

    # x, y = (minx, miny)
    # geom_array = []
    # densities = []

    w = maxx - minx
    h = maxy - miny

    n_cells_x = 256
    n_cells_y = 256
    cell_x = w / n_cells_x
    cell_y = h / n_cells_y

    x_coords = np.arange(minx, maxx, cell_x)
    y_coords = np.arange(miny, maxy, cell_y)

    geometries = [geometry.box(x, y, x + cell_x, y + cell_y) for y in y_coords for x in x_coords]
    fishnet = gpd.GeoDataFrame(geometry=geometries, crs=gdf.crs)

    mean_density = gdf['density'].mean()

    sindex = gdf.sindex

    densities = []

    for i, cell in fishnet.iterrows():
        cell_geom = cell['geometry']
        possible_matches_index = list(sindex.intersection(cell_geom.bounds))
        possible_matches = gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(cell_geom)]
        if not precise_matches.empty:
            density = precise_matches['density'].max()
        else:
            density = mean_density
        densities.append(density)

    fishnet['density'] = densities

    # cell_x_num = 0
    # cell_y_num = 0
    #
    # while y < maxy:
    #     while x < maxx:
    #         print(f"Cell at (x, y): ({cell_x_num}, {cell_y_num})")
    #         geom = geometry.Polygon([(x,y), (x, y+cell_y), (x+cell_x, y+cell_y), (x+cell_x, y), (x, y)])
    #         geom_array.append(geom)
    #
    #         cell_x_num += 1
    #
    #         # Check if cell intersects with gdf
    #         cell = gpd.GeoDataFrame(geometry=[geom]).set_crs('EPSG:3857')
    #         intersection = gpd.overlay(cell, gdf, how='intersection')
    #
    #         if len(intersection) > 0:
    #             print("Adding intersection density")
    #             density = intersection['density'].iloc[0]
    #         else:
    #             print("Adding mean density")
    #             density = mean_density
    #
    #         densities.append(density)
    #         x += cell_x
    #     x = minx
    #     y += cell_y
    #     cell_y_num += 1
    #     cell_x_num = 0
    #
    # fishnet = gpd.GeoDataFrame(geom_array, columns=['geometry']).set_crs('EPSG:3857')
    # fishnet['density'] = densities
    fishnet.plot(facecolor='none', ax=ax, alpha=0.5, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)

    # joined = gpd.sjoin(fishnet, gdf, op='intersects')
    # joined.plot(color='w', ax=ax, zorder=0, edgecolor='0', linewidth=0.1, legend=False)

    # fishnet['density'] = joined['density'].fillna(mean_density)

    OrRd = plt.get_cmap('OrRd')
    trunced_OrRd = truncate_colormap(OrRd, 0.2)

    fishnet.plot(column='density', cmap=trunced_OrRd, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=True)

    # density = joined.groupby(['index_right']).agg({'density': 'mean'})

    # non_con.plot(color='r', ax=ax, edgecolor='0', linewidth=0.1, legend=False)
    # dorling.plot(color='blue', ax=ax, edgecolor='0', linewidth=0.1, legend=False)

    # New Dorling Method Testing
    # radius_regions.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    # radius_regions.loc[neighbours_idx].plot(color='b', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    # radius_regions.loc[[current_region]].plot(color='r', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)
    # regions.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)

    # diffusion.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)

    # places_df.plot(ax=ax, edgecolor='0', linewidth=0.1)

    # geo_pop.plot(column='Population', cmap=new_cmap, ax=ax, edgecolor='0', linewidth=0.1, legend=True)
    # ax.set_title('Population of LADs', fontdict={'fontsize': '15', 'fontweight': '3'})

    # Plot Figure
    plt.savefig("./out/fishnet_256_density_2.png", dpi=1200)


