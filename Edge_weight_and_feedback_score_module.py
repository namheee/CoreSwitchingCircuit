import bisect

import matplotlib.pyplot as plt
import numpy as np

from Ftb_calculation_module import calculate_Ftb_of_feedback
from Cycle_analysis import Find_cycles_containing_the_node

def get_links_with_positive_edge_weight(links, profile_before_perturb, profile_after_perturb):
    """return a dictionary, 
    whose key is link with positive edge weight 
    and value is the edge weight of the link"""
    average_node_activity_difference = Average_node_activity_difference(profile_before_perturb, profile_after_perturb)
    return average_node_activity_difference.filter_links_with_edge_weight_positive(links)

def plot_histogram_and_pdf_of_positive_edge_weights(links_positiveedgeweight_map, bins=10):
    """
    Plot a histogram and probability density function (PDF) for values in a given dictionary.

    Parameters:
    - links_positiveedgeweight_map (dict): A dictionary where values are floating-point numbers between 0 and 1.
    - bins (int): Number of bins for the histogram (default is 10).

    The function displays:
    1. A histogram showing the distribution of values.
    2. A probability density function (PDF) estimated using a normalized histogram.
    """
    # Extract values from the dictionary
    values = np.array(list(links_positiveedgeweight_map.values()))

    # Compute percentiles (90%, 80%, ..., 10%)
    percentiles = [np.percentile(values, p) for p in range(90, 0, -10)]

    # Create histogram data (density=True normalizes the histogram to represent a probability distribution)
    hist_values, bin_edges = np.histogram(values, bins=bins, density=True)

    # Compute bin centers for smooth PDF plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Plot the histogram
    plt.figure(figsize=(8, 6))
    plt.bar(bin_centers, hist_values, width=(bin_edges[1] - bin_edges[0]), alpha=0.6, color="blue", label="Histogram")

    # Plot the probability density function (PDF)
    plt.plot(bin_centers, hist_values, marker="o", linestyle="-", color="red", label="Estimated PDF")

    # Annotate percentiles on the x-axis
    for p, x_val in zip(range(90, 0, -10), percentiles):
        plt.axvline(x_val, color="gray", linestyle="--", alpha=0.7)
        plt.text(x_val, max(hist_values) * 0.9, f"{p}%", fontsize=10, rotation=90, verticalalignment="top")

    # Add labels and title
    plt.xlabel("edge weight", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title("Histogram and Estimated PDF\nof positive edge weights", fontsize=14)

    # Add legend
    plt.legend(fontsize=12)

    # Show grid
    plt.grid(True, linestyle="--", alpha=0.7)

    # Display the plot
    plt.show()

class Average_node_activity_difference:
    def __init__(self, profile_before_perturb, profile_after_perturb):
        self.profile_before_perturb = profile_before_perturb
        self.profile_after_perturb = profile_after_perturb
        
        self.node_averagenodeactivitydiff_map = {}
        for node, state_after in self.profile_after_perturb.node_averagenodeactivity_map.items():
            state_before = self.profile_before_perturb.node_averagenodeactivity_map[node]
            self.node_averagenodeactivitydiff_map[node] = state_after - state_before
    
    def __repr__(self):
        return "{} - {}".format(self.profile_after_perturb, self.profile_before_perturb)

    def filter_links_with_edge_weight_positive(self, links):
        links_with_edge_weight_positive = {}
        for link in links:
            regulator = link[0]
            target = link[2]
            modality = link[1]
            
            regulator_diff_state = self.node_averagenodeactivitydiff_map[regulator]
            target_diff_state = self.node_averagenodeactivitydiff_map[target]
            diff_mul =  regulator_diff_state * target_diff_state

            if diff_mul >0 and modality == '+':
                links_with_edge_weight_positive[link] = abs(diff_mul)
            elif diff_mul <0 and modality == '-':
                links_with_edge_weight_positive[link] = abs(diff_mul)
        
        return links_with_edge_weight_positive
    

def _convert_links_in_cycle(cycle, links_map):
    """links_map is dictionary 
    mapping link with form (source, target) to link with form (source, sign, target)
    
    when the links in cycle has form of (source,target),
    convert the form of each link to form (source, sign, target)"""
    cycle_new = []
    for link in cycle:
        link_with_sign = links_map[link]
        cycle_new.append(link_with_sign)
    return tuple(cycle_new)    

# def get_feedbacks_having_feedback_score_higher_than_threshold(links_edgeweight_map, max_len=None, feedback_score_threshold=0):
#     """links_edgeweight_map is dictionary 
#     whose key is link having positive edge weight and value is edge weight of the link.

#     find all cycles whose length is less than or equal to max_len.
#     if max_len is None, find all cycles.

#     among found feedbacks, 
#     select feedbacks having feedback score higher than threshold and return them"""
#     links_map = {(x[0],x[-1]):x for x in links_edgeweight_map}
#     links_to_check = list(links_map)
#     all_nodes = set()
#     for link in links_to_check:
#         all_nodes.add(link[0])
#         all_nodes.add(link[-1])

#     feedback_collector = Feedbacks_over_feedback_score_threshold(links_edgeweight_map.copy(), feedback_score_threshold)
    
#     while all_nodes:
#         node_selected = all_nodes.pop()

#         calculator = Find_cycles_containing_the_node(node_selected, links_to_check)
#         cycles = calculator.find_cycles(algorithm="simple", max_len=max_len, return_node_form=False)
#         for cycle in cycles:
#             cycle_with_signs = _convert_links_in_cycle(cycle, links_map)
#             feedback_collector.put_new_feedback(cycle_with_signs)

#         links_to_check = [link for link in links_to_check if node_selected not in link]
#         all_nodes = set()
#         for link in links_to_check:
#             all_nodes.add(link[0])
#             all_nodes.add(link[-1])
    
#     return feedback_collector

def get_feedbacks_having_feedback_score_higher_than_threshold(links_edgeweight_map, max_len=None, feedback_score_threshold=0):
    """links_edgeweight_map is dictionary 
    whose key is link having positive edge weight and value is edge weight of the link.

    find all cycles whose length is less than or equal to max_len.
    if max_len is None, find all cycles.

    among found feedbacks, 
    select feedbacks having feedback score higher than threshold and return them
    
    
    이것은 리비전 동안만 사용하기 위해 pickle 기능을 넣은 함수.
    특정 max_len 값에 대해, 나온 결과를 저장하고, 다음 리비전에서는 저장된 결과를 불러와서 사용할 수 있도록 함."""

    import pickle
    import os
    import time
    time_now = time.time()
    pickle_file_name = "feedbacks_with_max_len_{}.pickle".format(str(max_len))
    pickle_address = os.path.join(r"D:\LG화학 project\논문 작성\241223 GPB 1차 revision\tmp_data_erasable_feedback_collectors_pickles", pickle_file_name)
    if os.path.exists(pickle_address):
        print("pickle file found. loading feedbacks from the pickle file")
        flag=True
        with open(pickle_address, 'rb') as f:
            node_cycles_map = pickle.load(f)
    else:
        print("no pickle file found. calculating feedbacks from scratch")
        flag=False
        node_cycles_map = {}

    links_map = {(x[0],x[-1]):x for x in links_edgeweight_map}
    links_to_check = list(links_map)
    all_nodes = set()
    for link in links_to_check:
        all_nodes.add(link[0])
        all_nodes.add(link[-1])

    feedback_collector = Feedbacks_over_feedback_score_threshold(links_edgeweight_map.copy(), feedback_score_threshold)
    
    while all_nodes:
        node_selected = all_nodes.pop()
        print("calculating feedbacks containing node: ", node_selected)
        print("the time taken until now: ", time.time()-time_now)

        
        if flag == True:
            if node_selected in node_cycles_map:
                cycles = node_cycles_map[node_selected]
            else:
                continue
        else:
            calculator = Find_cycles_containing_the_node(node_selected, links_to_check)
            cycles = calculator.find_cycles(algorithm="simple", max_len=max_len, return_node_form=False)
            node_cycles_map[node_selected] = cycles
        
        for cycle in cycles:
            cycle_with_signs = _convert_links_in_cycle(cycle, links_map)
            feedback_collector.put_new_feedback(cycle_with_signs)

        links_to_check = [link for link in links_to_check if node_selected not in link]
        all_nodes = set()
        for link in links_to_check:
            all_nodes.add(link[0])
            all_nodes.add(link[-1])
    
    if not flag:
        with open(pickle_address, 'wb') as f:
            pickle.dump(node_cycles_map, f)
    print("time taken for feedback calculation: ", time.time()-time_now)
    return feedback_collector

class Feedbacks_over_feedback_score_threshold:
    def __init__(self, link_edgeweight_map, feedback_score_threshold):
        self.feedback_score_threshold = feedback_score_threshold
        self.feedbacks = []
        self.feedback_scores = []
        #feedback score of self.feedbacks[i] is self.feedback_scores[i]
        self.link_edgeweight_map = link_edgeweight_map
    
    def put_new_feedback(self, new_feedback):
        """for a new feedback, check whether the feedback score of the feedback 
        is bigger than self.feedback_score_threshold.

        if the feedback has feedback score higher than the threshold,
        put the feedback to self.feedbacks, maintaining order according to feedback scores
        (the smaller the feedback score, the earlier the position in the self.feedbacks)
        """
        cycle_type, feedback_score = self.analyze_cycle(new_feedback)

        if cycle_type == -1:
            #negative feedback! ignore
            return
        else:
            #positive feedback
            if feedback_score < self.feedback_score_threshold:
                #discard feedback having feedback score less than the threshold
                return
            index_to_insert = bisect.bisect_right(self.feedback_scores, feedback_score)
            self.feedbacks.insert(index_to_insert, new_feedback)
            self.feedback_scores.insert(index_to_insert, feedback_score)
    
    def select_feedback_to_be_added_to_current_PPR(self, current_PPR:set,
                                                   links, links_edgeweight_map, 
                                                   except_nodes=set(), verbose=False):
        index_Ftb_map = {}
        for index, feedback in enumerate(self.feedbacks):
            nodes_in_feedback = self.get_nodes_in_cycle(feedback)
            if nodes_in_feedback.intersection(except_nodes):
                #ignore feedbacks containing nodes in exceot_nodes
                continue
            if current_PPR.intersection(nodes_in_feedback):
                #selected feedback should be overlapped to the current_PPR
                nodes_not_in_currnet_PPR = nodes_in_feedback.difference(current_PPR)
                if nodes_not_in_currnet_PPR:
                    #if nodes_not_in_currnet_PPR is empty, 
                    #the feedback is already contained the current_PPR
                    Ftb_of_feedback = calculate_Ftb_of_feedback(current_PPR, feedback, links, links_edgeweight_map)
                    index_Ftb_map[index] = Ftb_of_feedback
        
        if verbose:
            print("Ftbs are: ")
            for index, Ftb in index_Ftb_map.items():
                print(self.feedbacks[index], Ftb)
        
        if len(index_Ftb_map) == 0:
            print("\n\n********************************************")
            print("current PPR: {}".format(current_PPR))
            print("No extendable feedback can be found in the current PPR.")
            print("Please relax the feedback search criteria and recalculate.\n")
            print("Alternatively, the significant regulators may not be included in the feedback.")
            print("In this case, perform network reduction to modify the network so that the regulators of the phenotype node are included in the feedback.")
            print("********************************************")
            raise
        
        Ftb_min = min(index_Ftb_map.values())
        indexes_with_min_Ftb = []
        for index, Ftb in index_Ftb_map.items():
            if Ftb == Ftb_min:
                indexes_with_min_Ftb.append(index)
        
        if len(indexes_with_min_Ftb) >= 2:
            #print("more than two feedbacks have same min Ftb")
            #print("select feedback with maximum feedback score")
            maximum_feedback_score = -999
            index_with_max_feedback_score = None
            for index in indexes_with_min_Ftb:
                feedback = self.feedbacks[index]
                _, feedback_score = self.analyze_cycle(feedback)
                if feedback_score > maximum_feedback_score:
                    maximum_feedback_score = feedback_score
                    index_with_max_feedback_score = index
            
            feedback_selected = self.feedbacks[index_with_max_feedback_score]
        else:
            feedback_selected = self.feedbacks[indexes_with_min_Ftb[0]]
            
        return feedback_selected, Ftb_min
    
    def analyze_cycle(self, feedback):
        """return cycle_type and feedback_score
        if the cycle_type is 1, then the feedback is positive feedback.
        if the cycle_type is -1, then the feedback is negative feedback."""
        cycle_type = 1
        edge_weight_sum = 0
        for link in feedback:
            modality = link[1]
            if modality == '-':
                cycle_type *= -1
            elif modality == '+':
                pass
            else:
                raise("{} link has insane modality symbol".format(link))
            edge_weight_sum += self.link_edgeweight_map[link]
        feedback_score = edge_weight_sum/len(feedback)

        return cycle_type, feedback_score
    
    def get_nodes_in_cycle(self, feedback):
        """get set of nodes composing the feedback"""
        nodes_in_feedback = set()
        for link in feedback:
            nodes_in_feedback.add(link[0])
            nodes_in_feedback.add(link[-1])
        
        return nodes_in_feedback