

def tool_extract_node_arrays(cfi, nodes):
    """Extracts node arrays from a dictionary of nodes. 

    Args:
        cfi: Comfit instance.
    
    Returns:
        A dictionary of node arrays.
    """

    node_arrays = {}

    # Check if there are nodes to be treated, if not return an empty dictionary
    if not nodes:
        return node_arrays
    
    node_arrays['velocity_given'] = 'velocity' in nodes[0].keys()

    x_coordinates = []
    y_coordinates = []
    z_coordinates = []

    x_coordinates_positive = []
    y_coordinates_positive = []
    z_coordinates_positive = []

    x_coordinates_negative = []
    y_coordinates_negative = []
    z_coordinates_negative = []

    node_arrays['charge_given'] = 'charge' in nodes[0].keys()
    node_arrays['velocity_given'] = 'velocity' in nodes[0].keys()
    node_arrays['Burgers_vector_given'] = 'Burgers_vector' in nodes[0].keys()
    node_arrays['tangent_vector_given'] = 'tangent_vector' in nodes[0].keys()

    for node in nodes:
        x_coordinates.append(node['position'][0])
        y_coordinates.append(node['position'][1] if cfi.dim > 1 else None)
        z_coordinates.append(node['position'][2] if cfi.dim > 2 else None)

    node_arrays['x_coordinates'] = x_coordinates
    node_arrays['y_coordinates'] = y_coordinates if cfi.dim > 1 else None
    node_arrays['z_coordinates'] = z_coordinates if cfi.dim > 2 else None

    if node_arrays['charge_given']:

        for node in nodes:

            if node['charge'] > 0:
                x_coordinates_positive.append(node['position'][0])
                y_coordinates_positive.append(node['position'][1] if cfi.dim > 1 else None)
                z_coordinates_positive.append(node['position'][2] if cfi.dim > 2 else None)

            else:
                x_coordinates_negative.append(node['position'][0])
                y_coordinates_negative.append(node['position'][1] if cfi.dim > 1 else None)
                z_coordinates_negative.append(node['position'][2] if cfi.dim > 2 else None)
            
        node_arrays['x_coordinates_positive'] = x_coordinates_positive
        node_arrays['y_coordinates_positive'] = y_coordinates_positive if cfi.dim > 1 else None
        node_arrays['z_coordinates_positive'] = z_coordinates_positive if cfi.dim > 2 else None

        node_arrays['x_coordinates_negative'] = x_coordinates_negative
        node_arrays['y_coordinates_negative'] = y_coordinates_negative if cfi.dim > 1 else None
        node_arrays['z_coordinates_negative'] = z_coordinates_negative if cfi.dim > 2 else None

    for property in ['velocity', 'tangent_vector', 'Burgers_vector']:
        node_arrays[property + '_given'] = property in nodes[0].keys()

        if node_arrays[property + '_given']:
            x_coordinates = []
            y_coordinates = []
            z_coordinates = []

            for node in nodes:
                x_coordinates.append(node[property][0])
                y_coordinates.append(node[property][1] if cfi.dim > 1 else None)
                z_coordinates.append(node[property][2] if cfi.dim > 2 else None)
            
            node_arrays[property + '_x_coordinates'] = x_coordinates
            node_arrays[property + '_y_coordinates'] = y_coordinates if cfi.dim > 1 else None
            node_arrays[property + '_z_coordinates'] = z_coordinates if cfi.dim > 2 else None

    return node_arrays
        
        
