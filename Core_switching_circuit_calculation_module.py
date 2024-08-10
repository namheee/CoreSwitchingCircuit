
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






    

def show_edgeweights_of_edges_connecting_phenotype_and_regulators(phenotype_nodes:set,
                                                                link_edgeweight_map):
    print("regulator\tsign\tphenotype_node\tedge_weight")
    for link, edge_weight in link_edgeweight_map.items():
        target_node = link[-1]
        if target_node in phenotype_nodes:
            regulator_node = link[0]
            sign = link[1]
            print("{:>9}\t{:>4}\t{:>14}\t{:>11.5f}".format(regulator_node, str(sign), target_node, edge_weight))

def select_significant_regulators_for_phenotype_nodes(phenotype_nodes:set,
                                                      link_edgeweight_map,
                                                      edgeweight_threshold):
    significant_regulators = set()
    for link, edge_weight in link_edgeweight_map.items():
        target_node = link[-1]
        if target_node in phenotype_nodes:
            if edge_weight > edgeweight_threshold:
                regulator_node = link[0]
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