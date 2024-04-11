import networkx as nx
import matplotlib.pyplot as plt

def visualize_connections(component_map):
    """
    Visualizes the connections between components in a component map.

    Parameters:
    - component_map (ComponentMap): The component map containing the components and connections.

    Returns:
    - None
    """
    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes for each component
    for component in component_map.components.values():
        G.add_node(component.name)

    # Add edges for each connection
    for component_name, ports in component_map.connections.items():
        for port_name, (connected_component_name, _) in ports.items():
            if port_name in component_map.components[component_name].output_ports:
                G.add_edge(component_name, connected_component_name)

    # Draw the graph
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=500, font_size=10, edge_color='gray', arrows=True)

    # Show the plot
    plt.show()
