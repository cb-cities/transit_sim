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

### crs
project_crs=None

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
        self.station_id_route_dict = {getattr(row, 'node_id'): getattr(
            row, 'route_stop_id').split('-')[0] for row in self.all_nodes.itertuples()}
        
        ### create graph for shortest path calculations
        self.network_g = interface.from_dataframe(self.all_links, 'start_nid', 'end_nid', 'initial_weight')

        ### create geometry dictionary for visualization
        self.station_locations = {getattr(row, 'route_stop_id'): getattr(row, 'geometry') for row in self.all_nodes.itertuples()}
        #self.link_sections = {'{}-{}'.format(getattr(row, 'route_stop_id'), getattr(row, 'next_route_stop_id')):
        #    getattr(row, 'geometry').interpolate(0.2, 0.3) for row in self.all_links.itertuples()}
        self.station_cx = dict()
        self.station_cy = dict()
        self.link_cx = dict()
        self.link_cy = dict()
        self.link_ux = dict()
        self.link_uy = dict()
        for row in self.all_nodes.itertuples():
            self.station_cx[getattr(row, 'route_stop_id')] = getattr(row, 'geometry').x
            self.station_cy[getattr(row, 'route_stop_id')] = getattr(row, 'geometry').y
        for row in self.all_links.itertuples():
            link_nm = '{}-{}'.format(getattr(row, 'route_stop_id'), getattr(row, 'next_route_stop_id'))
            link_geom = getattr(row, 'geometry')
            self.link_cx[link_nm] = link_geom.interpolate(0.2, 0.3).x
            self.link_cy[link_nm] = link_geom.interpolate(0.2, 0.3).y
            self.link_ux[link_nm] = link_geom.coords[-1][0] - link_geom.coords[0][0]
            self.link_uy[link_nm] = link_geom.coords[-1][1] - link_geom.coords[0][1]

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
        ### add final destination of the trip
        for row in schedule_table.groupby('trip_id').tail(1).itertuples():
            schedule_list.append((getattr(row, 'trip_id'),
                                  getattr(row, 'next_arrival_time'), getattr(row, 'next_arrival_time')+30, 'stop', 
                                  getattr(row, 'next_route_stop_id')))
        self.schedule_df = pd.DataFrame(schedule_list, columns=['trip_id', 'time', 'next_time', 'status', 'location'])
        self.schedule_df = self.schedule_df.sort_values(by=['trip_id', 'time'], ascending=True)
    
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
        virtual_nodes['route_stop_id'] = virtual_nodes['stop_id'].apply(lambda x: 'all-{}'.format(x))
        virtual_nodes['type'] = 'station'
        all_nodes = pd.concat([all_nodes, virtual_nodes[all_nodes.columns]])
        all_nodes['node_id'] = np.arange(all_nodes.shape[0])
        all_nodes = gpd.GeoDataFrame(
            all_nodes, crs=project_crs, 
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
        self.schedule_df['location_id'] = self.schedule_df['location'].map(station_nm_id_dict).fillna(-1)
        all_links['start_nid'] = all_links['route_stop_id'].map(station_nm_id_dict)
        all_links['end_nid'] = all_links['next_route_stop_id'].map(station_nm_id_dict)
        all_links['initial_weight'] = 1.0
        all_links['geometry'] = all_links.apply(
            lambda x: 'LINESTRING({} {}, {} {})'.format(
                x['stop_lon'], x['stop_lat'], x['next_stop_lon'], x['next_stop_lat']
            ), axis=1)
        all_links = all_links[['route_stop_id', 'next_route_stop_id', 'start_nid', 'end_nid',
                               'initial_weight', 'geometry']]
        all_links = gpd.GeoDataFrame(all_links, crs=project_crs, geometry=all_links['geometry'].map(loads))
        return all_nodes, all_links
    
    def schedule_and_network_from_gtfs(self, stop_times_file, trips_file, stops_file, service_id):
        ### read GTFS files
        stop_times_table = pd.read_csv(stop_times_file)
        trips_table = pd.read_csv(trips_file)
        stops_table = pd.read_csv(stops_file)
        
        ### check format
        stops_table = stops_table.drop_duplicates(subset=['stop_id'])
        
        ### merge the tables
        schedule_table = stop_times_table[['trip_id', 'arrival_time', 'departure_time', 'stop_id']]
        ### assign a route code to individual trains
        trips_table = trips_table[trips_table['service_id']==service_id]
        schedule_table = pd.merge(schedule_table, trips_table[['trip_id', 'route_id']],
                                   how='right', on='trip_id')
        ### assign a route-stop code to individual stops
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
        schedule_table = gpd.GeoDataFrame(schedule_table, crs=project_crs, 
                                           geometry=[Point(xy) for xy in zip(schedule_table.stop_lon, 
                                                                             schedule_table.stop_lat)])
        if project_crs is not None:
            schedule_table = schedule_table.to_crs(3857)
        ### calculate shift factor
        minx, miny, maxx, maxy = schedule_table.total_bounds
        shift_factor = (maxx - minx)/10000
        schedule_table['stop_x'] = schedule_table.geometry.x + shift_factor*schedule_table['route_id'].map(route_seq_dict)
        schedule_table['stop_y'] = schedule_table.geometry.y + shift_factor*schedule_table['route_id'].map(route_seq_dict)
        schedule_table['geometry'] = [Point(xy) for xy in zip(schedule_table.stop_x, schedule_table.stop_y)]
        if project_crs is not None:
            schedule_table = schedule_table.to_crs(4326)
        schedule_table['stop_lon'] = schedule_table.geometry.x
        schedule_table['stop_lat'] = schedule_table.geometry.y
        
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
        train_positions = self.schedule_df.loc[self.schedule_df['current_location']=='current', 
                                           ['trip_id', 'status', 'location']].copy()
        train_positions['cx'] = np.where(train_positions['status']=='stop',
                                              train_positions['location'].map(network.station_cx),
                                              train_positions['location'].map(network.link_cx))
        train_positions['cy'] = np.where(train_positions['status']=='stop',
                                              train_positions['location'].map(network.station_cy),
                                              train_positions['location'].map(network.link_cy))
        train_positions['ux'] = np.where(train_positions['status']=='stop',
                                              0, train_positions['location'].map(network.link_ux))
        train_positions['uy'] = np.where(train_positions['status']=='stop',
                                              0, train_positions['location'].map(network.link_uy))
        return train_positions
        
class Travelers():
    def __init__(self):
        self.travelers_df = None
        # self.travelers_paths = dict() ### graph path
        # self.travelers_key_stops_list = []
        self.travelers_key_stops = None ### boarding, alignthing, transfering key stops
       
    def random_od(self, all_nodes=None, num_travelers=1):
        
        od_nodes = all_nodes['route_stop_id'].str.split('-').str[0]=='all'
        traveler_origins = np.random.choice(all_nodes.loc[od_nodes, 'node_id'], num_travelers)
        traveler_destins = np.random.choice(all_nodes.loc[od_nodes, 'node_id'], num_travelers)
        #traveler_origins = [241]
        #traveler_destins = [267]
        self.travelers_df = pd.DataFrame({
            'origin_nid': traveler_origins, 'destin_nid': traveler_destins})
        self.travelers_df = self.travelers_df[self.travelers_df['origin_nid'] != self.travelers_df['destin_nid']].copy()
        self.travelers_df['traveler_id'] = np.arange(self.travelers_df.shape[0])
        self.travelers_df['departure_time'] = np.random.randint(3600*6, 3600*10, self.travelers_df.shape[0])
        #self.travelers_df['departure_time'] = 26664-120
    
    def set_initial_status(self, station_id_nm_dict):
        ### initialize traveler_df
        self.travelers_df['traveler_status'] = 0 ### {0: 'pretrip', 1: 'walking', 2: 'platform', 3: 'train', 4: 'arrival'}
        self.travelers_df['update_time'] = 0
        self.travelers_df['association'] = -11
        self.travelers_df['next_station_id'] = -111
        self.travelers_df['init_boarding_time'] = 1e7
        self.travelers_df['final_alighting_time'] = 1e7
    
    def find_routes(self, network_g, station_id_nm_dict, station_id_route_dict):
        ### Group by origin and destinations
        travelers_key_stops_list = []
        
        for (traveler_origin, traveler_destin), travelers in self.travelers_df.groupby(['origin_nid', 'destin_nid']):
            traveler_origin = int(traveler_origin)
            traveler_destin = int(traveler_destin)
            ### find paths using Dijkstra's algorithm
            sp = network_g.dijkstra(traveler_origin, traveler_destin)
            sp_dist = sp.distance(traveler_destin)

            if sp_dist > 10e7:
                sp.clear()
                traveler_path = []
                key_stops = []
                print(traveler)
            else:
                sp_path = sp.route(traveler_destin)
                ### only record when changing line
                travelers_key_stops_list += [(traveler_origin, traveler_destin, nid) for (
                    start_nid, end_nid) in sp_path for nid in [start_nid, end_nid] 
                             if station_id_route_dict[start_nid]!=station_id_route_dict[end_nid]]
                sp.clear()
        
        self.travelers_key_stops = pd.DataFrame(travelers_key_stops_list, 
                                                columns=['origin_nid', 'destin_nid', 'current_stop_id'])
        self.travelers_key_stops['next_stop_id'] = self.travelers_key_stops.groupby(
            ['origin_nid', 'destin_nid'])['current_stop_id'].shift(-1)
        self.travelers_key_stops = self.travelers_key_stops[~pd.isnull(self.travelers_key_stops['next_stop_id'])]
        self.travelers_key_stops = self.travelers_df.merge(
            self.travelers_key_stops, how='left', on =['origin_nid', 'destin_nid'])[['traveler_id', 'current_stop_id', 'next_stop_id']]
        self.travelers_key_stops['next_stop_id'] = self.travelers_key_stops['next_stop_id'].astype(int)
            
    # def find_next_station(self, x):
    #     try:
    #         next_station = self.travelers_key_stops[x['traveler_id']][x['association']]
    #     except KeyError:
    #         next_station = None
    #     return next_station
    
    # @profile
    def traveler_update(self, network, trains, t, time_step_size=20, train_capacity=1460, transfer_time=120, exit_walking_time=120):

        ### get current train locations
        train_locations = trains.schedule_df.loc[trains.schedule_df['current_location']=='current']
        ### we are only interested in trains stop at platforms, as it is when travelers board or alight
        stop_trains = train_locations.loc[train_locations['status']=='stop'].copy()
        stop_trains['location_id'] = stop_trains['location_id'].astype(int)
        ### lookup dictionary 1: from platform_id to trip_id
        stop_train_locations_dict = {getattr(train, 'location_id'): 
                                     getattr(train, 'trip_id') for train in stop_trains.itertuples()}
        ### lookup dictionary 2: from trip_id to platform_id
        stop_trip_ids_dict = {getattr(train, 'trip_id'): 
                               getattr(train, 'location_id') for train in stop_trains.itertuples()}
        ### lookup dictionary 3: platform to stopped train total capacity
        ### hard coded capacity!!!
        stop_trains = stop_trains.merge(
            self.travelers_df[self.travelers_df['traveler_status']==3].groupby(
                'association').size().to_frame('train_occupancy'), left_on='trip_id', right_index=True, how='left')
        stop_trains['train_occupancy'] = stop_trains['train_occupancy'].fillna(0)
        stop_train_capacities_dict = {getattr(train, 'location_id'): 
                                     (train_capacity-getattr(train, 'train_occupancy')) for train in stop_trains.itertuples()}
        # if stop_trains.shape[0]>0:
        #     display(stop_trains)
        
        ### load departure travelers
        ### (1) change status from "pretrip" to "walking"
        ### (2) change association from "None" to "origin_nid"
        ### (3) change next_stop from "None" to next key stop
        departure_travelers = (
            self.travelers_df['departure_time']<=t) & (
            self.travelers_df['departure_time']>t-time_step_size) & (
            self.travelers_df['traveler_status']==0)
        self.travelers_df.loc[departure_travelers, 'traveler_status'] = 1
        self.travelers_df.loc[departure_travelers, 'update_time'] = t
        self.travelers_df.loc[departure_travelers, 'association'] = self.travelers_df.loc[departure_travelers, 'origin_nid'].values
        self.travelers_df.loc[departure_travelers, 'next_station_id'] = self.travelers_df.loc[departure_travelers].set_index(
            ['traveler_id', 'association']).join(self.travelers_key_stops[self.travelers_key_stops['traveler_id'].isin(
                self.travelers_df.loc[departure_travelers, 'traveler_id'])].set_index(['traveler_id', 'current_stop_id']), 
            how='left', on=['traveler_id', 'association'])['next_stop_id'].values
        
        ### transfer
        ### conditions: status is "walking" and transfer time is 2 minutes
        walking_travelers = (
            self.travelers_df['traveler_status']==1) & (
            self.travelers_df['update_time']<=t-transfer_time)
        self.travelers_df.loc[walking_travelers, 'traveler_status'] = 2
        self.travelers_df.loc[walking_travelers, 'update_time'] = t
        self.travelers_df.loc[walking_travelers, 'association'] = self.travelers_df.loc[walking_travelers, 'next_station_id']
        self.travelers_df.loc[walking_travelers, 'next_station_id'] = self.travelers_df.loc[walking_travelers].merge(
            self.travelers_key_stops[self.travelers_key_stops['traveler_id'].isin(self.travelers_df.loc[walking_travelers, 'traveler_id'])], how='left', left_on=['traveler_id', 'association'], right_on=['traveler_id', 'current_stop_id'])['next_stop_id'].values
        
        ### arrival
        arrival_travelers = self.travelers_df['traveler_status']==2
        self.travelers_df['arrival_tmp'] = False
        self.travelers_df.loc[arrival_travelers, 'arrival_tmp'] = self.travelers_df.loc[
            arrival_travelers, 'association'] == self.travelers_df.loc[arrival_travelers, 'destin_nid']
        arrival_travelers = arrival_travelers.values & self.travelers_df['arrival_tmp'].values
        self.travelers_df.loc[arrival_travelers, 'traveler_status'] = 4
        self.travelers_df.loc[arrival_travelers, 'update_time'] = t + exit_walking_time-transfer_time
        self.travelers_df.loc[arrival_travelers, 'final_alighting_time'] = t
        #self.travelers_df.loc[arrival_travelers, 'association'] = None
        self.travelers_df.loc[arrival_travelers, 'next_station_id'] = None
        
        ### aboard: travelers ready to aboard
        ### check if next stop is covered
        ### condition 1: status is "platform"
        board_travelers = self.travelers_df['traveler_status']==2
        ### condition 2: a train at the platform
        self.travelers_df['aboard_tmp'] = False
        self.travelers_df.loc[board_travelers, 'aboard_tmp'] = self.travelers_df.loc[board_travelers, 'association'].isin(stop_train_locations_dict.keys())
        board_travelers = board_travelers.values & self.travelers_df['aboard_tmp'].values
        ### condition 3: within capacity -- note it is hard-coded capacity
        self.travelers_df['remain_cap'] = 0
        self.travelers_df.loc[board_travelers, 'remain_cap'] = self.travelers_df.loc[
            board_travelers, 'association'].map(stop_train_capacities_dict).fillna(-10)
        self.travelers_df['order'] = 1e7
        self.travelers_df.loc[self.travelers_df['remain_cap']>0, 'order'] = self.travelers_df.loc[
            self.travelers_df['remain_cap']>0].sort_values(
            by='update_time', ascending=True).groupby('association').cumcount()
        board_travelers = board_travelers & (self.travelers_df['order'] < self.travelers_df['remain_cap']).values
        ### (1) change next stop according to key routes
        self.travelers_df.loc[board_travelers, 'next_station_id'] = self.travelers_df.loc[board_travelers].merge(
            self.travelers_key_stops[self.travelers_key_stops['traveler_id'].isin(self.travelers_df.loc[board_travelers, 'traveler_id'])],  
            how='left', left_on=['traveler_id', 'association'], right_on=['traveler_id', 'current_stop_id'])['next_stop_id'].values
        ### (2) change association from platform to trip_id
        new_board = self.travelers_df.loc[board_travelers].copy()
        new_board['boarding_platform'] = new_board['association']
        new_board['boarding_train'] = new_board['association'].replace(stop_train_locations_dict)
        self.travelers_df.loc[board_travelers, 'association'] = self.travelers_df.loc[board_travelers, 
                                                                     'association'].replace(stop_train_locations_dict)
        ### (3) change status from "platform" to "train"
        self.travelers_df.loc[board_travelers, 'traveler_status'] = 3
        #self.travelers_df['traveler_status'] = np.where(board_travelers, 'train', self.travelers_df['traveler_status'])
        new_board['prev_time'] = new_board['update_time']
        self.travelers_df.loc[board_travelers, 'update_time'] = t
        self.travelers_df.loc[board_travelers, 'init_boarding_time'] = np.minimum(t, self.travelers_df.loc[board_travelers, 'init_boarding_time'])
        new_board['boarding_time'] = t
        new_board['waiting_time'] = new_board['boarding_time'] - new_board['prev_time']
        new_board = new_board[['traveler_id', 'boarding_platform', 'boarding_train', 
                               'prev_time', 'boarding_time', 'waiting_time']]

        ### alight: travelers ready to get off the train
        ### conditions: (1) status is "train"
        alight_travelers = self.travelers_df['traveler_status']== 3
        ### conditions: (2) train location is alighting location
        self.travelers_df['train_location_id'] = -1
        self.travelers_df.loc[alight_travelers, 'train_location_id'] = self.travelers_df.loc[
            alight_travelers, 'association'].replace(stop_trip_ids_dict)
        # if (t>=22100) and (t<22700): print(stop_trip_ids_dict)
        alight_travelers = alight_travelers.values & (
            self.travelers_df['train_location_id']==self.travelers_df['next_station_id']).values
        ### (1) change association from trip_id to platform
        self.travelers_df.loc[alight_travelers, 'association'] = self.travelers_df.loc[alight_travelers, 'next_station_id']
        ### (2) change status from "train" to "walking"
        self.travelers_df.loc[alight_travelers, 'traveler_status'] = 1
        #self.travelers_df['traveler_status'] = np.where(alight_travelers, 'walking', self.travelers_df['traveler_status'])
        self.travelers_df.loc[alight_travelers, 'update_time'] = t
        ### (3) change next stop according to key routes
        # self.travelers_df.loc[alight_travelers, 'next_station'] = self.travelers_df.loc[
        #     alight_travelers].apply(lambda x: self.find_next_station(x), axis=1)
        # self.travelers_df.loc[alight_travelers, 'next_station'] = self.travelers_df.loc[alight_travelers].merge(
        #     self.travelers_key_stops[self.travelers_key_stops['traveler_id'].isin(self.travelers_df.loc[alight_travelers, 'traveler_id'])], 
        #     how='left', left_on=['traveler_id', 'association'], right_on=['traveler_id', 'current_stop'])['next_stop']
        self.travelers_df.loc[alight_travelers, 'next_station_id'] = self.travelers_df.loc[alight_travelers].set_index(['traveler_id', 'association']).join(
            self.travelers_key_stops[self.travelers_key_stops['traveler_id'].isin(self.travelers_df.loc[alight_travelers, 'traveler_id'])].set_index(
                ['traveler_id', 'current_stop_id']), 
            how='left', on=['traveler_id', 'association'])['next_stop_id'].values
        
        return new_board
    
    def get_all_traveler_positions(self, train_positions):
        
        ### get traveler locations
        ### group by station or train location
        traveler_locations = self.travelers_df.groupby(
            ['traveler_status', 'association']).size().to_frame(
            name='num_travelers').reset_index(drop=False)
        
        travelers_on_trains = traveler_locations[traveler_locations['traveler_status']==3]
        travelers_on_trains['association'] = travelers_on_trains['association'].astype(int)
        travelers_on_trains = travelers_on_trains.merge(train_positions[['trip_id', 'cx', 'cy']], 
                                                        how='left', left_on='association', right_on='trip_id')
        return traveler_locations