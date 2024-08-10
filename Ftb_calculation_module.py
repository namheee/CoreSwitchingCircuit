def calculate_Ftb_of_feedback(perturbation_propagated_region, 
                              feedback_to_calculate,
                              links, links_positiveedgeweight_map):
    nodes_in_feedback = set()
    for link in feedback_to_calculate:
        nodes_in_feedback.add(link[0])
        nodes_in_feedback.add(link[-1])
    
    nodes_outside_part = nodes_in_feedback.difference(perturbation_propagated_region)

    node_incomingedges_map = {}
    for node in nodes_outside_part:
        node_incomingedges_map[node] = []
        for link in links:
            if link[-1] == node:
                node_incomingedges_map[node].append(link)
    
    node_feedbacklink_map = {}
    node_unfrustratedlink_map = {}
    node_frustratedlink_map = {}
    node_fromunknownlink_map = {}
    for node, incominglinks in node_incomingedges_map.items():
        node_feedbacklink_map[node] = []
        node_unfrustratedlink_map[node] = []
        node_frustratedlink_map[node] = []
        node_fromunknownlink_map[node] = []
        for link in incominglinks:
            if link in feedback_to_calculate:
                node_feedbacklink_map[node].append(link)
                continue
            if link[0] not in perturbation_propagated_region:
                node_fromunknownlink_map[node].append(link)
                #not used
                continue
            if link in links_positiveedgeweight_map:
                node_unfrustratedlink_map[node].append(link)
            else:
                node_frustratedlink_map[node].append(link)
    
    node_Ntb_map = {}
    for node, incominglinks in node_incomingedges_map.items():
        union_feedback_unfrustrated = set(node_feedbacklink_map[node]+node_unfrustratedlink_map[node])
        tendency = (len(union_feedback_unfrustrated)-len(node_frustratedlink_map[node]))/len(incominglinks)
        node_Ntb_map[node] = 1-tendency
    
    Ftb = sum(node_Ntb_map.values())

    return Ftb