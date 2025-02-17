import matplotlib.pyplot as plt

from Edge_weight_and_feedback_score_module import get_links_with_positive_edge_weight


def get_initially_perturbred_region_for_perturbation(links,
                                                     profile_before_perturbation,
                                                     profile_after_perturbation,
                                                     perturbed_nodes:"list of node names",
                                                     edge_weight_threshold):
    linkswithpositiveedgeweight_edgeweight_map = get_links_with_positive_edge_weight(links, profile_before_perturbation, profile_after_perturbation)
    stack_of_nodes = perturbed_nodes.copy()
    initially_perturbed_region = set(perturbed_nodes)
    while stack_of_nodes:
        node_to_check = stack_of_nodes.pop(0)
        for link, edgeweight in linkswithpositiveedgeweight_edgeweight_map.items():
            if link[0] == node_to_check:
                if edgeweight > edge_weight_threshold:
                    target_node = link[-1]
                    if target_node not in initially_perturbed_region:
                        initially_perturbed_region.add(target_node)
                        stack_of_nodes.append(target_node)
    
    initially_perturbed_region.difference_update(perturbed_nodes)
    return initially_perturbed_region


def find_downstream_nodes_and_max_edge_weights(edges_dict, start_nodes):
    """
    Finds downstream nodes from the given start nodes and prints the maximum edge weight used at each step.
    
    :param edges_dict: A dictionary where keys are edges (source, sign, target), and values are edge weights.
                       Key format: (source, sign, target)
    :param start_nodes: A set of source node names to start the traversal from.
    """
    current_nodes = set(start_nodes)  # Nodes to explore in the current step
    visited_edges = set()  # Track visited edges to avoid reprocessing

    Max_weights = []
    
    while current_nodes:
        downstream_nodes = set()  # Nodes discovered in this step
        max_weight = float('-inf')  # Maximum weight in the current step
        
        for (source, sign, target), weight in edges_dict.items():
            if source in current_nodes and (source, sign, target) not in visited_edges:
                downstream_nodes.add(target)  # Add connected downstream node
                visited_edges.add((source, sign, target))  # Mark edge as visited
                max_weight = max(max_weight, weight)  # Update the maximum weight
        
        if downstream_nodes:
            Max_weights.append(max_weight)
        else:
            break
        
        # Move to the next step
        current_nodes = downstream_nodes
    
    return Max_weights

# def trace_max_edge_weights_for_initially_perturbed_region(links, profile_before_perturbation,
#                                                           profile_after_perturbation1, profile_after_perturbation2,
#                                                           perturbed_nodes1, perturbed_nodes2):
#     linkswithpositiveedgeweight_edgeweight_map1 = get_links_with_positive_edge_weight(links, profile_before_perturbation, profile_after_perturbation1)
#     linkswithpositiveedgeweight_edgeweight_map2 = get_links_with_positive_edge_weight(links, profile_before_perturbation, profile_after_perturbation2)

#     max_weights1 = find_downstream_nodes_and_max_edge_weights(linkswithpositiveedgeweight_edgeweight_map1, perturbed_nodes1)
#     max_weights2 = find_downstream_nodes_and_max_edge_weights(linkswithpositiveedgeweight_edgeweight_map2, perturbed_nodes2)

#     name_of_perturbation1 = "perturb node(s) {}".format(','.join(perturbed_nodes1))
#     name_of_perturbation2 = "perturb node(s) {}".format(','.join(perturbed_nodes2))
#     # Create x values (indices)
#     x1 = range(len(max_weights1))
#     x2 = range(len(max_weights2))

#     # Find the segment with the largest decrease
#     max_drop1 = 0
#     max_segment1 = (0, 1)  # Indices of the segment with the largest drop
#     for i in range(1, len(max_weights1)):
#         drop1 = max_weights1[i - 1] - max_weights1[i]
#         if drop1 > max_drop1:
#             max_drop1 = drop1
#             max_segment1 = (i - 1, i)
    
#     max_drop2 = 0
#     max_segment2 = (0, 1)  # Indices of the segment with the largest drop
#     for i in range(1, len(max_weights2)):
#         drop2 = max_weights2[i - 1] - max_weights2[i]
#         if drop2 > max_drop2:
#             max_drop2 = drop2
#             max_segment2 = (i - 1, i)

#     # Create the plot
#     plt.figure(figsize=(8, 6))

#     # Plot the first line graph
#     plt.plot(x1, max_weights1, marker='o', linestyle='--', color='blue', label=name_of_perturbation1)

#     # Highlight the segment with the largest drop
#     start_idx, end_idx = max_segment1
#     plt.plot(
#         [start_idx, end_idx],
#         [max_weights1[start_idx], max_weights1[end_idx]],
#         color="red",
#         linewidth=3,
#         label=f"Largest drop in {name_of_perturbation1}: {max_drop1:.2f}",
#     )

def trace_max_edge_weights_for_initially_perturbed_region(links, profile_before_perturbation,
                                                          profile_after_perturbation,
                                                          perturbed_nodes):
    linkswithpositiveedgeweight_edgeweight_map = get_links_with_positive_edge_weight(links, profile_before_perturbation, profile_after_perturbation)
    
    max_weights = find_downstream_nodes_and_max_edge_weights(linkswithpositiveedgeweight_edgeweight_map, perturbed_nodes)

    name_of_perturbation = "Initially perturbed node(s) by {}".format(','.join(perturbed_nodes))

    # Create x values (indices)
    x = range(len(max_weights))

    max_weights_len_is_longer_than_1 = len(max_weights) > 1

    # Find the segment with the largest decrease
    if max_weights_len_is_longer_than_1:
        max_drop = 0
        max_segment = (0, 1)  # Indices of the segment with the largest drop
        for i in range(1, len(max_weights)):
            drop = max_weights[i - 1] - max_weights[i]
            if drop > max_drop:
                max_drop = drop
                max_segment = (i - 1, i)


    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot the first line graph
    plt.plot(x, max_weights, marker='o', linestyle='--', color='blue', label=name_of_perturbation)

    # Highlight the segment with the largest drop
    if max_weights_len_is_longer_than_1:
        start_idx, end_idx = max_segment
        plt.plot(
            [start_idx, end_idx],
            [max_weights[start_idx], max_weights[end_idx]],
            color="red",
            linewidth=3,
            label=f"The most negative derivative in {name_of_perturbation}: {-1*max_drop:.2f}",
        )
    
    # Add text labels for each point
    for i in range(len(x)):
        plt.text(
            x[i], max_weights[i],
            f"{max_weights[i]:.3f}", 
            fontsize=10, ha="center", color="black"
        )

    # Add labels and title
    plt.xlabel("Depth", fontsize=12)
    plt.ylabel("Max edge weight", fontsize=12)
    plt.title("Comparison of Max edge weights \nfor {}".format(name_of_perturbation), fontsize=16)

    # Add legend
    plt.legend(fontsize=12)

    # Add grid
    plt.grid(True, linestyle='--', alpha=0.7)

    # Show the plot
    plt.show()


    

# def show_edgeweights_of_edges_connecting_phenotype_and_regulators(phenotype_nodes:set,
#                                                                 link_edgeweight_map):
#     print("regulator\tsign\tphenotype_node\tedge_weight")
#     for link, edge_weight in link_edgeweight_map.items():
#         target_node = link[-1]
#         if target_node in phenotype_nodes:
#             regulator_node = link[0]
#             sign = link[1]
#             print("{:>9}\t{:>4}\t{:>14}\t{:>11.5f}".format(regulator_node, str(sign), target_node, edge_weight))

def visualize_largest_decrease_in_edge_weight_across_depth(phenotype_nodes:set, link_edgeweight_map, combi_target=""):
    values = []
    for link, edge_weight in link_edgeweight_map.items():
        target_node = link[-1]
        regulator = link[0]
        if target_node in phenotype_nodes:
            if regulator not in phenotype_nodes:
                values.append(edge_weight)

    # Sort the values in descending order
    sorted_values = sorted(values, reverse=True)

    # Create x values (indices)
    x = range(len(sorted_values))

    len_of_sorted_values_is_longer_than_1 = len(sorted_values) > 1

    # Find the segment with the largest decrease
    if len_of_sorted_values_is_longer_than_1:
        max_drop = 0
        max_segment = (0, 1)  # Indices of the segment with the largest drop
        for i in range(1, len(sorted_values)):
            drop = sorted_values[i - 1] - sorted_values[i]
            if drop > max_drop:
                max_drop = drop
                max_segment = (i - 1, i)

    # Plot the line graph
    plt.figure(figsize=(8, 6))
    plt.plot(x, sorted_values, marker="o", linestyle="--", color="blue", label="Edge weights")

    # Highlight the segment with the largest drop
    if len_of_sorted_values_is_longer_than_1:
        start_idx, end_idx = max_segment
        plt.plot(
            [start_idx, end_idx],
            [sorted_values[start_idx], sorted_values[end_idx]],
            color="red",
            linewidth=3,
            label=f"The most negative derivative: {-1*max_drop:.3f}",
        )

    # Add text labels for each point
    for i in range(len(x)):
        plt.text(
            x[i], sorted_values[i],
            f"{sorted_values[i]:.3f}", 
            fontsize=10, ha="center", color="black"
        )

    # Add labels and title
    plt.xlabel("Edges connecting regulators and phenotypes, \nsorted in descending order by edge weight values.", fontsize=12)
    plt.ylabel("Edge weight", fontsize=12)
    plt.title("Edge weights between regulator and phenotype {}".format(combi_target), fontsize=16)

    # Add legend
    plt.legend(fontsize=12)

    # Add grid
    plt.grid(True, linestyle="--", alpha=0.7)

    # Show the plot
    plt.show()

def select_significant_regulators_for_phenotype_nodes(phenotype_nodes:set,
                                                      link_edgeweight_map,
                                                      edgeweight_threshold):
    significant_regulators = set()
    for link, edge_weight in link_edgeweight_map.items():
        target_node = link[-1]
        regulator_node = link[0]
        if target_node in phenotype_nodes:
            if regulator_node not in phenotype_nodes:
                if edge_weight > edgeweight_threshold:
                    significant_regulators.add(regulator_node)

    return significant_regulators
        



def calculate_PPR_seq_and_related_info(feedback_collector_to_use, 
                                     links, 
                                     links_edgeweight_map,
                                     fixed_nodes,
                                     initially_perturbed_region, 
                                     regulators_destination, 
                                     verbose=True):
    current_PPR = initially_perturbed_region.copy() #PPR=perturbation-propagated region
    PPR_seq = []
    feedback_seq = []
    Ftb_seq = []
    score_seq = []

    while True:
        PPR_seq.append(current_PPR)
        if current_PPR.issuperset(regulators_destination):
            break
    
        #print("PPR in this step is ",current_PPR)
        feedback_to_add, Ftb = feedback_collector_to_use.select_feedback_to_be_added_to_current_PPR(current_PPR,links, links_edgeweight_map, fixed_nodes, verbose)
        _, feedback_score = feedback_collector_to_use.analyze_cycle(feedback_to_add)
        feedback_seq.append(feedback_to_add)
        Ftb_seq.append(Ftb)
        score_seq.append(feedback_score)
        nodes_in_feedback = feedback_collector_to_use.get_nodes_in_cycle(feedback_to_add)
        current_PPR = current_PPR.union(nodes_in_feedback)
    
    return PPR_seq, feedback_seq, Ftb_seq, score_seq