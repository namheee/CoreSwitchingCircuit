import sys

from Data_reading_module import *
from Core_switching_circuit_calculation_module import *
from Edge_weight_and_feedback_score_module import *

def calculate_core_switching_circuit(argv):
    command_file_address = argv[1]
    command_file_obj = Command_file(command_file_address)
    links = read_links(command_file_obj.net_structure_address)
    
    target_nodes_of_drug = set(command_file_obj.drug_target_nodes)
    target_nodes_of_comb = set(command_file_obj.combi_target_nodes)
    perturbed_nodes = target_nodes_of_drug.union(target_nodes_of_comb)

    phenotype_nodes = set(command_file_obj.phenotype_nodes)

    profile_nominal = Average_node_activity_profile_for_perturbation("nominal")
    profile_nominal.read_from_file(command_file_obj.average_node_activity_files["nominal"])
    profile_drug_only = Average_node_activity_profile_for_perturbation("drug only")
    profile_drug_only.read_from_file(command_file_obj.average_node_activity_files["drug only"])
    profile_comb_target_only = Average_node_activity_profile_for_perturbation("combination target only")
    profile_comb_target_only.read_from_file(command_file_obj.average_node_activity_files["combination target only"])
    profile_drug_and_comb_target = Average_node_activity_profile_for_perturbation("drug and combination target")
    profile_drug_and_comb_target.read_from_file(command_file_obj.average_node_activity_files["drug and combination target"])

    links_positiveedgeweight_map = get_links_with_positive_edge_weight(links, profile_nominal, profile_drug_and_comb_target)

    print("The graph below shows the max weights of the edges forming the downstream of the drug target when only the drug target was processed, and the max weights of the edges forming the downstream of the combination target when only the combination target was processed.")
    trace_max_edge_weights_for_initially_perturbed_region(links, profile_nominal,
                                                          profile_drug_only, profile_comb_target_only,
                                                          list(target_nodes_of_drug), list(target_nodes_of_comb))
    while True:
        print("Please set the edge weight threshold to determine initially perturbed region.\nThe range is between 0 and 1. Once you enter a value,\nthe expected initially perturbed region will be displayed.")
        print("I recommend setting the edge weight threshold for the initially perturbed region to satisfy the section where the drop in max weight between the two line graphs of the graph is large across the Depth.")
        threshold_for_initially_perturbed_region = float(input())
        print("Regarding {}, the expected initially perturbed retion is as follows.".format(','.join(list(target_nodes_of_drug))))
        initially_perturbed_region_for_drug_only = get_initially_perturbred_region_for_perturbation(links, 
                                                                                                    profile_nominal, 
                                                                                                    profile_drug_only, 
                                                                                                    list(target_nodes_of_drug), 
                                                                                                    threshold_for_initially_perturbed_region)
        print(initially_perturbed_region_for_drug_only)
        print("\n")
        print("Regarding {}, the expected initially perturbed retion is as follows.".format(','.join(list(target_nodes_of_comb))))
        initially_perturbed_region_for_combi_target_only = get_initially_perturbred_region_for_perturbation(links, 
                                                                                                    profile_nominal, 
                                                                                                    profile_comb_target_only, 
                                                                                                    list(target_nodes_of_comb), 
                                                                                                    threshold_for_initially_perturbed_region)
        print(initially_perturbed_region_for_combi_target_only)

        print("To confirm this value, enter 'y'. To input a different value, enter 'n'.")
        if input().lower() == 'y':
            break
    initially_perturbed_region = initially_perturbed_region_for_drug_only.union(initially_perturbed_region_for_combi_target_only)

    print("\n\nPlease set the edge weight threshold to select significant regulators for the phenotype.\nBelow are the edge weights between the phenotype and regulators. \nRegulators with an edge weight equal to or greater than the edge weight threshold will be selected as significant regulators.")
    show_edgeweights_of_edges_connecting_phenotype_and_regulators(phenotype_nodes, links_positiveedgeweight_map)
    threshold = float(input(":"))
    significant_regulators = select_significant_regulators_for_phenotype_nodes(phenotype_nodes, 
                                                                               links_positiveedgeweight_map,
                                                                               threshold)
    print("The selected significant regulators are as follows.")
    print(significant_regulators)

    
    print("\n\nPlease set the maximum feedback length that can be calculated. \nChoosing a value that is too large may result in a long feedback search time.")
    max_len = int(input())
    print("\n\nPlease set the feedback score threshold for the core switching circuit search. The range is from 0 to 1. \nChoosing a value that is too low may result in insufficient memory for the computation.")
    threshold_for_feedback_score = float(input())

    feedback_collector = get_feedbacks_having_feedback_score_higher_than_threshold(links_positiveedgeweight_map, max_len, threshold_for_feedback_score)

    
    PPR_seq, feedback_seq, Ftb_seq, score_seq = calculate_PPR_seq_and_related_info(feedback_collector, 
                                     links, 
                                     links_positiveedgeweight_map,
                                     perturbed_nodes,
                                     initially_perturbed_region, 
                                     significant_regulators, 
                                     verbose=False)
    
    with open(command_file_obj.save_address, 'w') as f:
        f.write("initially perturbed region\n{}\n\n".format(PPR_seq[0]))
        for step in range(len(feedback_seq)):
            f.write("{}th step\n".format(step))
            f.write("PPR:{}\n".format(PPR_seq[step+1]))
            f.write("added feedback:{}\n".format(feedback_seq[step]))
            f.write("Ftb of the feedback:{}\n".format(Ftb_seq[step]))
            f.write("feedback score of the feedback:{}\n\n".format(score_seq[step]))

    





if __name__ == "__main__":
    calculate_core_switching_circuit(sys.argv)