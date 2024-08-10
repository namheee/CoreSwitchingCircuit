# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 16:05:05 2021

@author: jwKim
"""
import numpy as np
import itertools
import SCC_decomposition

class Find_cycles_containing_the_node:
    #use Johnson's algorithm
    def __init__(self, s_node_contained_to_cycle, l_t_links):
        self.l_nodes = SCC_decomposition.get_SCC_containing_the_node(s_node_contained_to_cycle, l_t_links)
        self.links_of_SCC = self._extract_links_of_subnet(self.l_nodes, l_t_links)
        self.dict_i_node_t_i_targets = {}
        #self.dict_i_node_t_i_targets = {index of node:(index of target nodes1,)}
        self.dict_i_count_for_each_node = {}
        self.dict_i_count_for_each_node_saved = {}
        self._refine_topology_information(l_t_links)
        
        self.array_blocked = np.zeros((len(self.l_nodes),), dtype=bool)
        self.l_connectable_to_the_node = []#has same length to l_flow. if l_connectable_to_0[i] = True: l_flow[:i+1] can reach to 0(root node)
        
        self.i_index_of_the_node = self.l_nodes.index(s_node_contained_to_cycle)
        self.l_flow = []
        self.l_l_i_cycles = []
        self.dict_blocked_link = {i:set([]) for i in range(len(self.l_nodes))}#if dic_blocked_link[i] contains j, then j->i link is blocked
        self.i_max_len = len(self.l_nodes)
        
    @staticmethod
    def _extract_links_of_subnet(nodes_of_subnet, links):
        links_extracted = []
        for link in links:
            if (link[0] in nodes_of_subnet) and (link[-1] in nodes_of_subnet):
                links_extracted.append(link)
        return links_extracted
        
    def _refine_topology_information(self, l_t_links):
        for i,s_node in enumerate(self.l_nodes):
            self.dict_i_node_t_i_targets[i] = tuple([self.l_nodes.index(t_link[1]) for t_link in l_t_links if (t_link[0] == s_node) and (t_link[1] in self.l_nodes)])
            self.dict_i_count_for_each_node[i] = len(self.dict_i_node_t_i_targets[i]) -1
            self.dict_i_count_for_each_node_saved[i] = len(self.dict_i_node_t_i_targets[i]) -1
            
    def _extend_flow(self, i_node):
        self.l_flow.append(i_node)
        self.l_connectable_to_the_node.append(False)
        self.array_blocked[i_node] = True
    
    def _block_links(self, i_node):
        for i_node_end in self.dict_i_node_t_i_targets[i_node]:
            self.dict_blocked_link[i_node_end].add(i_node)

    def _unblock(self, i_node):
        if self.array_blocked[i_node]:#i is blocked:
            self.array_blocked[i_node] = False #reset of blocking
        if self.dict_blocked_link[i_node]:#dic_blocked_link[i_node] is set
            for i_node_start in list(self.dict_blocked_link[i_node]):#i_node_start->i_node links are blocked
                self.dict_blocked_link[i_node].discard(i_node_start)# reset i_node_start->i_node links
                self._unblock(i_node_start)
    
    def _passable_case(self, i_next_edge):
        i_node_next = self.dict_i_node_t_i_targets[self.l_flow[-1]][i_next_edge]
        self.dict_i_count_for_each_node[self.l_flow[-1]] -= 1
        if i_node_next == self.i_index_of_the_node:#it is a cycle containing 0
            self.l_connectable_to_the_node[-1] = True#l_flow[-1] node is connected to 0 node
            self.l_l_i_cycles.append(self.l_flow.copy())     
        elif (not self.array_blocked[i_node_next]) and (self.l_flow[-1] not in self.dict_blocked_link[i_node_next]):# not passed this node yet
            self._extend_flow(i_node_next)
    
    def _go_to_next_node(self, i_next_edge):
        """Find the next node from the current self.l_flow.
        If that node is the starting node, add the simple cycle to self.l_l_i_cycles and terminate.
        If the node is not blocked, extend self.l_flow to include that node and terminate.
        If the node is blocked and is not the starting node, terminate without making any changes to self.l_flow."""
        i_node_next = self.dict_i_node_t_i_targets[self.l_flow[-1]][i_next_edge]
        self.dict_i_count_for_each_node[self.l_flow[-1]] -= 1
        if i_node_next == self.i_index_of_the_node:#it is a cycle containing 0
            self.l_l_i_cycles.append(self.l_flow.copy())
        elif len(self.l_flow) >= self.i_max_len:
            return
        elif not self.array_blocked[i_node_next]:
            self._extend_flow(i_node_next)
    
    def _impassable_case(self):
        """Key point of this algorithm: 
        As the flow retraces its steps, if it determines that the flow cannot return to the starting point 
        (i.e., it cannot form a cycle, with b_end = False), 
        it blocks off all the links to ensure that this path will not be revisited."""
        i_end = self.l_flow.pop(-1)
        b_end = self.l_connectable_to_the_node.pop(-1)
        if self.l_flow:
            self.dict_i_count_for_each_node[i_end] = self.dict_i_count_for_each_node_saved[i_end]#reset the passage of links starting from this node
            if b_end:
                self._unblock(i_end)
                self.l_connectable_to_the_node[-1] = True
            else:
                self._block_links(i_end)
    
    def _back_from_node(self):
        i_end = self.l_flow.pop(-1)
        if self.l_flow:
            self.dict_i_count_for_each_node[i_end] = self.dict_i_count_for_each_node_saved[i_end]
            self.array_blocked[i_end] = False #reset of blocking
            
    def restore_link_form_from_node_form_feedback(self, l_i_cycle):
        cycle_copy = l_i_cycle.copy()
        cycle_copy.append(l_i_cycle[0])
        links_in_cycle = []
        for i in range(len(l_i_cycle)):
            node_from = self.l_nodes[cycle_copy[i]]
            node_to = self.l_nodes[cycle_copy[i+1]]
            links_this_part = []#같은 source, target 에 대해 2개 이상의 link 가 있을 경우를 대비하여.
            for link in self.links_of_SCC:
                if (link[0] == node_from) and (link[-1] == node_to):
                    links_this_part.append(link)
            links_in_cycle.append(links_this_part)
        
        cycles_link_form = []
        for link_form_feedback in itertools.product(*tuple(links_in_cycle)):
            cycles_link_form.append(link_form_feedback)
        
        return cycles_link_form
        
    
    def find_cycles(self, algorithm="Johnson", max_len=None, return_node_form=True):
        """max_len specifies the maximum length of simple cycles to search for on the network.
        For example, max_len=1 will only search for self-loops.
        If set to None, it means there is no limit.
        This is to prevent excessive computational load.
        Currently, the max_len function only works with the simple algorithm, and improvements are needed."""
        if max_len is None:
            self.i_max_len = len(self.l_nodes)
        elif max_len<1 or (type(max_len) != type(1)):
            raise(ValueError("max_len should be integer and bigger then 0, or should be None"))
        elif max_len != None:
            self.i_max_len = max_len
        
        if len(self.l_nodes) == 1:
            if self.i_index_of_the_node in self.dict_i_node_t_i_targets[self.i_index_of_the_node]:#self loop
                self.l_l_i_cycles.append([self.i_index_of_the_node])
        else:
            if algorithm =="Johnson":
                self._extend_flow(self.i_index_of_the_node)
                while self.l_flow:
                    i_next_edge = self.dict_i_count_for_each_node[self.l_flow[-1]]#마지막 node의 target nodes 중 갈 node의 선택.
                    if i_next_edge >=0:
                        self._passable_case(i_next_edge)
                    else:
                        self._impassable_case()
            elif algorithm == "simple":
                self._extend_flow(self.i_index_of_the_node)
                while self.l_flow:                    
                    i_next_edge = self.dict_i_count_for_each_node[self.l_flow[-1]]#마지막 node의 target nodes 중 갈 node의 선택
                    if i_next_edge >=0:
                        self._go_to_next_node(i_next_edge)
                    else:
                        self._back_from_node()
        if return_node_form:
            return [[self.l_nodes[i] for i in l_cycle] for l_cycle in self.l_l_i_cycles]
        else:#return as link form
            cycles = []
            for l_i_cycle in self.l_l_i_cycles:
                cycles.extend(self.restore_link_form_from_node_form_feedback(l_i_cycle))
            return cycles
        

