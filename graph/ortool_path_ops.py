

from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
import numpy as np
import os
from graph.visualize_utils import visualize_graph

def create_ortool_data_model(coordinate_list,edge_pairs,edge_density, num_vehicles=20,cutoff_length=10,relax_choice=True):
    """Stores the data for the problem."""
    max_Value = 99999999
    adj_matrix = np.ones([len(coordinate_list)+1,len(coordinate_list)+1])*max_Value#super big to avoid this connections
    #adj matrix must be an integer matrix
    for k in range(len(edge_pairs)):
        node1, node2 = edge_pairs[k]
        current_prob = edge_density[k]
        #edge pairs already begins with 1
        adj_matrix[node1,node2]= int(current_prob*100)
        adj_matrix[node2,node1]= int(current_prob*100)

    print("number of edges in the data: ",len(np.argwhere(adj_matrix!=max_Value)))
    print("distance mean %f"%np.mean(adj_matrix[adj_matrix!=max_Value]))

    #adj_matrix[adj_matrix==99999]=sys.maxsize
    adj_matrix[0,:]=0#use 0 as a pseudo search starting point
    adj_matrix[:,0]=0
    adj_matrix[np.arange(1,len(coordinate_list)+1),np.arange(1,len(coordinate_list)+1)]=0



    #continue add edges for remained nodes: 10 times of the actual distances.
    if relax_choice:#necessary to be able to search useful candidates
        max_distance = np.max(adj_matrix[adj_matrix!=max_Value])
        for i in range(1,len(coordinate_list)+1):
            for j in range(1,len(coordinate_list)+1):
                if adj_matrix[i,j]==max_Value:
                    tmp_coord1 = coordinate_list[int(i)-1]
                    tmp_coord2 = coordinate_list[int(j)-1]
                    cur_distance = 0
                    for kk in range(3):
                        cur_distance+=(tmp_coord1[kk]-tmp_coord2[kk])**2
                    cur_distance = np.sqrt(cur_distance)
                    if cur_distance>=cutoff_length:
                        adj_matrix[i,j] = int(2*max_distance+cur_distance*200)#big penalty for them
                    else:
                        adj_matrix[i,j] = int(max_distance+cur_distance*100)
        print("updated distance mean %f"%np.mean(adj_matrix[adj_matrix!=max_Value]))
    data = {}
    data['distance_matrix'] = adj_matrix#final_adj_matrix
    data['num_vehicles'] = num_vehicles
    print("number vehicles: ",num_vehicles,"adj matrix shape: ",len(adj_matrix))
    data['depot'] = 0
    return data

def print_solution_vrp(data, manager, routing, solution):
    """Prints solution on console."""
    print(f'Objective: {solution.ObjectiveValue()}')
    max_route_distance = 0
    All_Route_List = []
    for vehicle_id in range(data['num_vehicles']):
        cur_Route_List =[]
        index = routing.Start(vehicle_id)
        plan_output = 'Route for vehicle {}:\n'.format(vehicle_id)
        route_distance = 0
        while not routing.IsEnd(index):
            plan_output += ' {} -> '.format(manager.IndexToNode(index))
            previous_index = index
            cur_Route_List.append(manager.IndexToNode(index))
            index = solution.Value(routing.NextVar(index))
            route_distance += routing.GetArcCostForVehicle(
                previous_index, index, vehicle_id)
        plan_output += '{}\n'.format(manager.IndexToNode(index))
        cur_Route_List.append(manager.IndexToNode(index))
        plan_output += 'Distance of the route: {}m\n'.format(route_distance)
        print(plan_output)
        All_Route_List.append(cur_Route_List)
        max_route_distance = max(route_distance, max_route_distance)
    print('Maximum of the route distances: {}m'.format(max_route_distance))
    return All_Route_List

def build_travel_path(coordinate_list,travel_list):
    final_coord =[]
    final_pair = []
    for k in range(len(travel_list)):
        travel_id = int(travel_list[k])
        if travel_id==0:
            continue
        travel_id -=1
        final_coord.append(coordinate_list[travel_id])
    for k in range(len(final_coord)-1):
        final_pair.append([k+1,k+2])
    return final_coord,final_pair
def ortools_build_path(save_path,coordinate_list,edge_pairs,edge_density,cutoff_length,relax_choice=True):

    max_Value = 99999999

    #number_vehicles = int(np.ceil(len(coordinate_list) / 100))
    number_vehicles = int(np.ceil(len(coordinate_list) / 20))#find 100 did not work so well for very big structures for new model
    if number_vehicles<=5:
        number_vehicles=5
    ortool_data = create_ortool_data_model(coordinate_list, edge_pairs, edge_density, number_vehicles,cutoff_length,relax_choice=relax_choice)
    distance_path = os.path.join(save_path, "distance.npy" )
    np.save(distance_path, np.array(ortool_data['distance_matrix']))
    distance_matrix = ortool_data['distance_matrix']
    drop_penalty= int(2*100*np.max(edge_density)+cutoff_length*200)#int(np.max(distance_matrix[distance_matrix!=max_Value]))
    print("drop penalty %d"%drop_penalty)
    # Create the routing index manager.
    manager = pywrapcp.RoutingIndexManager(len(ortool_data['distance_matrix']),
                                               ortool_data['num_vehicles'], ortool_data['depot'])
    routing = pywrapcp.RoutingModel(manager)

   # for vehicle_id in range(ortool_data['num_vehicles']):
    #    routing.ConsiderEmptyRouteCostsForVehicle(True, vehicle_id)

    def distance_callback(from_index, to_index):
        """Returns the distance between the two nodes."""
        # Convert from routing variable Index to distance matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return ortool_data['distance_matrix'][from_node][to_node]
        # Create Routing Model.

    print("start routing")

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)
    print("add capacity")
    # add dimension for vrp program
    dimension_name = 'Distance'
    routing.AddDimension(
        transit_callback_index,
        0,  # no slack
        max_Value,  # vehicle maximum travel distance
        True,  # start cumul to zero
        dimension_name)
    print("add penalty")
    # Allow to drop nodes.
    penalty = drop_penalty
    for node in range(1, len(ortool_data['distance_matrix'])):
        routing.AddDisjunction([manager.NodeToIndex(node)], penalty)
    # Setting first solution heuristic.
    print("add params")
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
         routing_enums_pb2.FirstSolutionStrategy.AUTOMATIC)
    # search_parameters.first_solution_strategy = (
    # routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)
    search_parameters.local_search_metaheuristic = (
        routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)
    search_parameters.time_limit.seconds = max(3600,int(1000*len(coordinate_list)/400)) # few hours search
    search_parameters.solution_limit = 1000
    search_parameters.log_search = True
    # Solve the problem.
    solution = routing.SolveWithParameters(search_parameters)
    travel_list = print_solution_vrp(ortool_data, manager, routing, solution)
    for j, cur_travel_list in enumerate(travel_list):
        travel_coord, travel_edge = build_travel_path(coordinate_list, cur_travel_list)
        visualize_graph(save_path,"travel_path%d"%j,travel_coord,travel_edge)
        tmp_save_path = os.path.join(save_path,"travel_record_%d.txt"%j)
        with open(tmp_save_path,'w') as file:
            for travel_id in cur_travel_list:
                file.write("%d\n"%(travel_id-1))
