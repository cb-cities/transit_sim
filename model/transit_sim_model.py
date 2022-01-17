import sys
import itertools
import numpy as np
import pandas as pd 
import geopandas as gpd
from shapely.wkt import loads
from shapely.geometry import Point

### shortest path
sys.path.insert(0, '/Users/bingyu')
from sp import interface

class Network():
    def __init__(self, all_nodes, all_links):
        
        ### nodes and links dataframe
        self.all_nodes = all_nodes
        self.all_links = all_links
        ### dictionary map of station name and ID
        self.station_nm_id_dict = {getattr(row, 'route_stop_id'): getattr(
            row, 'node_id') for row in self.all_nodes.itertuples()}
        self.station_id_nm_dict = {getattr(row, 'node_id'): getattr(
            row, 'route_stop_id') for row in self.all_nodes.itertuples()}
        
        ### create graph for shortest path calculations
        self.network_g = interface.from_dataframe(self.all_links, 'start_nid', 'end_nid', 'initial_weight')

        ### create geometry dictionary for visualization
        self.station_locations = {getattr(row, 'route_stop_id'): getattr(row, 'geometry') for row in self.all_nodes.itertuples()}
        self.link_sections = {'{}-{}'.format(getattr(row, 'route_stop_id'), getattr(row, 'next_route_stop_id')):
            getattr(row, 'geometry').interpolate(0.2, 0.3) for row in self.all_links.itertuples()}

class Trains():
    def __init__(self):
        self.schedule_df = None
    
    def add_schedule(self, schedule_table):
        ### convert the cleaned GTFS data (schedule_table) into class attribute (self.schedule)
        schedule_list = []
        for row in schedule_table.itertuples():
            schedule_list.append((getattr(row, 'trip_id'),
                                  getattr(row, 'arrival_time'), getattr(row, 'departure_time'), 'stop', 
                                  getattr(row, 'route_stop_id')))
            schedule_list.append((getattr(row, 'trip_id'),
                                  getattr(row, 'departure_time'), getattr(row, 'next_arrival_time'), 'on_link',
                                  '{}-{}'.format(getattr(row, 'route_stop_id'), 
                                                 getattr(row, 'next_route_stop_id'))))
        self.schedule_df = pd.DataFrame(schedule_list, columns=['trip_id', 'time', 'next_time', 'status', 'location'])
    
    def add_network(self, schedule_table):
        
        ### process network
        all_links = schedule_table.drop_duplicates(subset=['route_stop_id', 'next_route_stop_id'])
        all_links = all_links[['route_stop_id', 'next_route_stop_id', 
                               'stop_lon', 'stop_lat', 'next_stop_lon', 'next_stop_lat']].copy()

        ### create nodes
        all_nodes = pd.DataFrame(np.vstack(
            [all_links[['route_stop_id', 'stop_lon', 'stop_lat']].values,
             all_links[['next_route_stop_id', 'next_stop_lon', 'next_stop_lat']].values]),
            columns=['route_stop_id', 'stop_lon', 'stop_lat'])
        all_nodes = all_nodes.drop_duplicates(subset=['route_stop_id'])
        all_nodes['stop_id'] = all_nodes['route_stop_id'].apply(lambda x: x.split('-')[1])
        all_nodes['type'] = 'platform'
        ### station nodes
        virtual_nodes = all_nodes.groupby('stop_id').agg({'stop_lon': np.mean, 'stop_lat': np.mean}).reset_index(drop=False)
        virtual_nodes['stop_lon'] *= 0.999999
        virtual_nodes['route_stop_id'] = virtual_nodes['stop_id'].apply(lambda x: 'all_{}'.format(x))
        virtual_nodes['type'] = 'station'
        all_nodes = pd.concat([all_nodes, virtual_nodes[all_nodes.columns]])
        all_nodes['node_id'] = np.arange(all_nodes.shape[0])
        all_nodes = gpd.GeoDataFrame(
            all_nodes, crs='epsg:4326', 
            geometry=[Point(xy) for xy in zip(all_nodes.stop_lon, all_nodes.stop_lat)])
        station_nm_id_dict = {getattr(row, 'route_stop_id'): getattr(
            row, 'node_id') for row in all_nodes.itertuples()}
        station_id_nm_dict = {getattr(row, 'node_id'): getattr(
            row, 'route_stop_id') for row in all_nodes.itertuples()}

        ### add transfer links
        transfer_links = []
        for stop_id, grp in all_nodes.groupby('stop_id'):
            for (stop1, stop2) in list(itertools.permutations(grp.to_dict('records'), 2)):
                transfer_links.append([stop1['route_stop_id'], stop2['route_stop_id'], 
                                       stop1['stop_lon'], stop1['stop_lat'], stop2['stop_lon'], stop2['stop_lat']])
        transfer_links_df = pd.DataFrame(transfer_links, columns=['route_stop_id', 'next_route_stop_id',
                                                                 'stop_lon', 'stop_lat', 'next_stop_lon', 'next_stop_lat'])
        all_links = pd.concat([all_links, transfer_links_df])

        ### map stop names to node_ids
        all_links['start_nid'] = all_links['route_stop_id'].map(station_nm_id_dict)
        all_links['end_nid'] = all_links['next_route_stop_id'].map(station_nm_id_dict)
        all_links['initial_weight'] = 1.0
        all_links['geometry'] = all_links.apply(
            lambda x: 'LINESTRING({} {}, {} {})'.format(
                x['stop_lon'], x['stop_lat'], x['next_stop_lon'], x['next_stop_lat']
            ), axis=1)
        all_links = all_links[['route_stop_id', 'next_route_stop_id', 'start_nid', 'end_nid',
                               'initial_weight', 'geometry']]
        all_links = gpd.GeoDataFrame(all_links, crs='epsg:4326', geometry=all_links['geometry'].map(loads))
        return all_nodes, all_links
    
    def schedule_and_network_from_gtfs(self, stop_times_file, trips_file, stops_file):
        ### read GTFS files
        stop_times_table = pd.read_csv(stop_times_file)
        trips_table = pd.read_csv(trips_file)
        stops_table = pd.read_csv(stops_file)
        
        ### merge the tables
        schedule_table = stop_times_table[['trip_id', 'arrival_time', 'departure_time', 'stop_id']]
        ### assign a route code to individual trains
        schedule_table = pd.merge(schedule_table, trips_table[['trip_id', 'route_id']],
                                   how='left', on='trip_id')
        ### assign a route code to individual stops
        schedule_table['route_stop_id'] = schedule_table.apply(lambda x:
                                                               '{}-{}'.format(x['route_id'], x['stop_id']), axis=1)
        ### assign locations to individual stations
        schedule_table = pd.merge(schedule_table, 
                                    stops_table[['stop_id', 'stop_lon', 'stop_lat']], how='left', on='stop_id')
        ### shift line for better plotting
        route_seq_dict = dict()
        seq_id = 0
        for route_id, _ in schedule_table.sort_values(by='route_id', ascending=True).groupby('route_id'):
            route_seq_dict[route_id] = seq_id
            seq_id += 1
        schedule_table = gpd.GeoDataFrame(schedule_table, crs='epsg:4326', 
                                           geometry=[Point(xy) for xy in zip(schedule_table.stop_lon, 
                                                                             schedule_table.stop_lat)])
        schedule_table = schedule_table.to_crs(3857)
        schedule_table['stop_x'] = schedule_table.geometry.centroid.x + 10*schedule_table['route_id'].map(route_seq_dict)
        schedule_table['stop_y'] = schedule_table.geometry.centroid.y + 10*schedule_table['route_id'].map(route_seq_dict)
        schedule_table['geometry'] = [Point(xy) for xy in zip(schedule_table.stop_x, schedule_table.stop_y)]
        schedule_table = schedule_table.to_crs(4326)
        schedule_table['stop_lon'] = schedule_table.geometry.centroid.x
        schedule_table['stop_lat'] = schedule_table.geometry.centroid.y
        
        ### convert arrival and departure time to seconds since midnight
        schedule_table['arrival_time'] = schedule_table['arrival_time'].apply(
            lambda x: 3600*int(x.split(':')[0]) + 60*int(x.split(':')[1]) +
            int(x.split(':')[2]))
        schedule_table['departure_time'] = schedule_table['departure_time'].apply(
            lambda x: 3600*int(x.split(':')[0]) + 60*int(x.split(':')[1]) +
            int(x.split(':')[2]))
        ### add 30 seconds dwell time at stop if train arrival time = train departure time in GTFS
        schedule_table['departure_time'] = np.where(schedule_table['arrival_time']==schedule_table['departure_time'],
                                                schedule_table['departure_time']+30, schedule_table['departure_time'])
        
        ### link to next stops
        schedule_table = schedule_table.sort_values(by=['trip_id', 'arrival_time'], ascending=True)
        schedule_table['next_route_stop_id'] = schedule_table['route_stop_id'].shift(-1)
        schedule_table['next_stop_lon'] = schedule_table['stop_lon'].shift(-1)
        schedule_table['next_stop_lat'] = schedule_table['stop_lat'].shift(-1)
        schedule_table['next_trip_id'] = schedule_table['trip_id'].shift(-1)
        schedule_table['next_arrival_time'] = schedule_table['arrival_time'].shift(-1)
        schedule_table = schedule_table[schedule_table['trip_id']==schedule_table['next_trip_id']]
        
        ### create schedule and network
        self.add_schedule(schedule_table)
        all_nodes, all_links = self.add_network(schedule_table)
        return all_nodes, all_links
        
    def update_location_occupancy(self, t):
        self.schedule_df['current_location'] = np.where(
            self.schedule_df['time']>t, 'future', np.where(
            self.schedule_df['next_time']>t, 'current', 'past'))
        
    def get_all_train_positions(self, network):
        ### get train position
        train_positions = self.schedule_df[self.schedule_df['current_location']=='current'].copy()
        train_positions['geometry'] = np.where(train_positions['status']=='stop',
                                              train_positions['location'].map(network.station_locations),
                                              train_positions['location'].map(network.link_sections))
        train_positions = gpd.GeoDataFrame(train_positions, crs='epsg:4326', geometry=train_positions['geometry'])
        return train_positions
        
class Travelers():
    def __init__(self):
        self.travelers_df = None
        self.travelers_paths = dict() ### graph path
        self.travelers_key_stops = dict() ### boarding, alignthing, transfering key stops
       
    def random_od(self, all_nodes=None, num_travelers=1):
        
        #traveler_origins = np.random.choice(all_nodes['node_id'], num_travelers)
        #traveler_destins = np.random.choice(all_nodes['node_id'], num_travelers)
        traveler_origins = [9]
        traveler_destins = [22]
        self.travelers_df = pd.DataFrame({
            'origin_nid': traveler_origins, 'destin_nid': traveler_destins})
        self.travelers_df = self.travelers_df[self.travelers_df['origin_nid'] != self.travelers_df['destin_nid']].copy()
        self.travelers_df['traveler_id'] = np.arange(self.travelers_df.shape[0])
        self.travelers_df['departure_time'] = np.random.randint(25980, 26980, self.travelers_df.shape[0])
    
    def set_initial_status(self, station_id_nm_dict):
        ### initialize traveler_df
        self.travelers_df['traveler_status'] = 'pretrip'
        self.travelers_df['association'] = None 
        self.travelers_df['next_station'] = None
    
    def find_routes(self, network_g, station_id_nm_dict):
        ### find paths using Dijkstra's algorithm
        for traveler in self.travelers_df.itertuples():
            traveler_origin = getattr(traveler, 'origin_nid')
            traveler_destin = getattr(traveler, 'destin_nid')
            sp = network_g.dijkstra(traveler_origin, traveler_destin)
            sp_dist = sp.distance(traveler_destin)

            if sp_dist > 10e7:
                sp.clear()
                traveler_path = []
                key_stops = []
                print(traveler)
            else:
                sp_path = sp.route(traveler_destin)
                traveler_path = [
                    station_id_nm_dict[start_nid] for (start_nid, end_nid) in sp_path] + [
                    station_id_nm_dict[traveler_destin]]
                sp.clear()

                key_stops = []
                current_route = None ### current line number/name
                for stop in traveler_path[:-1]: ### always append the last element, so no need to iterate through
                    new_route = stop.split('-')[0]
                    if new_route != current_route:
                        key_stops.append(stop)
                        current_route = new_route
                key_stops.append(traveler_path[-1])

            self.travelers_paths[getattr(traveler, 'traveler_id')] = traveler_path
            self.travelers_key_stops[getattr(traveler, 'traveler_id')] = {
                xy[0]:xy[1] for xy in zip(key_stops, key_stops[1:])}
            
    def find_next_station(self, x):
        try:
            next_station = self.travelers_key_stops[x['traveler_id']][x['association']]
        except KeyError:
            next_station = None
        return next_station

    def traveler_update(self, network, trains, t):

        ### get current train locations
        train_locations = trains.schedule_df.loc[trains.schedule_df['current_location']=='current']
        ### we are only interested in trains stop at platforms, as it is when travelers board or alight
        stop_trains = train_locations.loc[train_locations['status']=='stop']
        ### lookup dictionary 1: from platform to trip_id
        stop_train_locations_dict = {getattr(train, 'location'): 
                                     getattr(train, 'trip_id') for train in stop_trains.itertuples()}
        ### lookup dictionary 2: from trip_id to platform
        stop_trip_ids_dict = {getattr(train, 'trip_id'): 
                               getattr(train, 'location') for train in stop_trains.itertuples()}
        
        ### load departure travelers
        ### (1) change status from "pretrip" to "platform"
        ### (2) change association from "None" to "origin_nid"
        ### (3) change next_stop from "None" to next key stop
        departure_travelers = (
            self.travelers_df['departure_time']<t) & (
            self.travelers_df['traveler_status']=='pretrip')
        self.travelers_df.loc[departure_travelers, 'traveler_status'] = 'platform'
        self.travelers_df.loc[departure_travelers, 'association'] = self.travelers_df.loc[
            departure_travelers, 'origin_nid'].map(network.station_id_nm_dict)
        self.travelers_df.loc[departure_travelers, 'next_station'] = self.travelers_df.loc[
            departure_travelers].apply(lambda x: self.find_next_station(x) , axis=1)
        
        ### aboard: travelers ready to aboard
        ### conditions: status is "platform"
        board_travelers = self.travelers_df['traveler_status']=='platform'
        ### (1) change association from platform to trip_id
        self.travelers_df.loc[board_travelers, 'association'] = self.travelers_df.loc[board_travelers, 
                                                                     'association'].replace(stop_train_locations_dict)
        ### (2) change status from "platform" to "train"
        self.travelers_df['traveler_status'] = np.where(
            board_travelers, 'train', self.travelers_df['traveler_status'])
        ### (3) change next stop according to key routes
        self.travelers_df.loc[board_travelers, 'next_station'] = self.travelers_df.loc[
            board_travelers].apply(lambda x: self.find_next_station(x) , axis=1)

        ### alight: travelers ready to get off the train
        ### conditions: (1) status is "train"
        alight_travelers = self.travelers_df['traveler_status']=='train'
        ### conditions: (2) train location is alighting location
        train_locations = self.travelers_df.loc[alight_travelers, 'association'].replace(stop_trip_ids_dict)
        alight_travelers = alight_travelers & (train_locations==self.travelers_df.loc[alight_travelers, 'next_station'])
        ### (1) change association from trip_id to platform
        self.travelers_df.loc[alight_travelers, 'association'] = self.travelers_df.loc[alight_travelers, 
                                                                     'association'].replace(stop_trip_ids_dict)
        ### (2) change status from "train" to "platform"
        self.travelers_df['traveler_status'] = np.where(alight_travelers, 'platform', self.travelers_df['traveler_status'])
        ### (3) change next stop according to key routes
        self.travelers_df.loc[alight_travelers, 'next_station'] = self.travelers_df.loc[
            alight_travelers, 'next_station'].apply(lambda x: self.find_next_station(x) , axis=1)
    
    def get_all_traveler_positions(self, station_locations=None, train_positions=None):
        
        ### get traveler locations
        ### group by station or train location
        traveler_locations = self.travelers_df.groupby(
            ['traveler_status', 'association']).size().to_frame(
            name='num_travelers').reset_index(drop=False)
        traveler_locations['association'] = traveler_locations['association'].astype(str)
        traveler_locations['geometry'] = None
        ### assign geometry to stations
        traveler_locations['geometry'] = traveler_locations['geometry'].map(station_locations)
        ### assign geometry to travelers on trains
        train_positions['trip_id'] = train_positions['trip_id'].astype(str)
        traveler_locations = pd.merge(
            traveler_locations, train_positions[['trip_id', 'geometry']], 
            how='left', left_on='association', right_on='trip_id', suffixes=['_platform', '_train'])
        traveler_locations['geometry'] = np.where(pd.isnull(traveler_locations['geometry_platform']), 
                                                    traveler_locations['geometry_train'],
                                                    traveler_locations['geometry_platform'])
        traveler_locations = traveler_locations[['association', 'num_travelers', 'geometry']]
        traveler_locations = gpd.GeoDataFrame(traveler_locations, crs='epsg:4326', 
                                                 geometry=traveler_locations['geometry'])
        return traveler_locations