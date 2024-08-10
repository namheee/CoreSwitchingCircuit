import json

class Command_file:
    def __init__(self, command_file_address):
        with open(command_file_address, 'r', encoding='UTF-8') as command_json:
            command_info = json.load(command_json)
        self.net_structure_address = command_info['network structure file address']
        self.drug_target_nodes = command_info['node names targeted by drug']
        self.combi_target_nodes = command_info['node names targeted by combination target']
        self.phenotype_nodes = command_info['phenotype node names']
        self.average_node_activity_files = command_info['average node activity profile address']
        self.save_address = command_info['save address']
    

        


def read_links(address):
    """read tsv file containing links information.
    it should have three columns, which are source, target, modality"""
    links = []
    with open(address, 'r') as f:
        f.readline() #column line
        for line in f:
            line=line.strip()
            line_splited = line.split('\t')
            source_node = line_splited[0]
            target_node = line_splited[1]
            sign = line_splited[2]
            if sign == "activation":
                links.append((source_node, '+', target_node))
            elif sign == "inhibition":
                links.append((source_node, '-', target_node))
            else:
                print(line)
                raise
    
    return links


class Average_node_activity_profile_for_perturbation:
    def __init__(self, perturbation):
        self.perturbation = perturbation
        self.node_averagenodeactivity_map = {}

    def __repr__(self):
        return self.perturbation

    def read_from_file(self, file_address):
        with open(file_address, 'r') as f:
            f.readline()
            #first line is column
            for line in f:
                line_splited = line.split('\t')
                node = line_splited[0]
                average_node_activity = float(line_splited[1])

                self.node_averagenodeactivity_map[node] = average_node_activity
