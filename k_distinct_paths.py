'''
Implementation of algorithm to find k distinct paths on a kinetic transition network (KTN).
The paths are distinct in the sense that they feature unique rate-limiting steps.
Input: A KTN of minima and transition state energies
Stage 1: compute edge costs
Stage 2: Dijkstra algorithm to find fastest path
Stage 3: iteratively block edge corresponding to rate-limiting step and use Marchetti-Spaccamela algorithm to find
         k-1 next-"fastest" paths
Output: k distinct (low-energy) pathways on the KTN
Daniel J. Sharpe (djs244)
Jul 2018
'''

import numpy as np
import heapq
import read_mindata

class K_Distinct_Paths(object):

    k_B = 1.9872036*10**(-3) # Boltzmann constant / kcal K^{-1} mol^{-1}

    def __init__(self):
        ### PARAMS - update as desired
        self.k = 25           # no. of paths to find
        self.M = float("inf") # a large number for effective blocking of edges
        self.T = 298.         # temperature / K (used to calculate inverse Boltzmann weights)
        self.s = 1            # start node (=1 for example min.data file)
        self.t = 17           # end node (=3345 for example min.data file, =17 for toy problem)
        self.costfunc = "noe_ts"    # choose a function for calculating the edge weights based on TS energies. Options:
                                    # noe_ts: Noe's scheme based on inverse Boltzmann weighting of TS energies
                                    # noe_b: Noe's scheme based on inverse Boltzmann weighting of "edge barriers"
                                    # evans: Evans' scheme based on a log-weighted adjacency matrix
        self.write_epath = False    # write Epath.<PATH> files Y/N
        self.write_mdf = False      # write min.data.fastest.<PATH> files Y/N
        self.write_mdf_all = True   # write min.data.fastest.all file Y/N

    '''Implementation of Dijkstra's algorithm using a min-priority queue. Used to find shortest path and to find
    corresponding shortest path tree for use in Marchetti-Spaccamella algorithm'''
    def dijkstra_ktn(self, G, min_energies=None, ts_energies=None):
        n_nodes = len(G)
        if self.write_mdf_all: self.in_mdf = {i: False for i in range(1,n_nodes+1)}
        self.dist = [float("inf")]*n_nodes
        self.dist[self.s-1] = 0. # start vertex has zero distance
        self.prev = [-1]*n_nodes # parents are initially null
        pq1 = K_Distinct_Paths.Priority_Queue()
        Q = pq1.pq # vertex set with specified priorities
        for i in range(1,n_nodes+1): pq1.add_with_priority(i,self.dist[i-1])
        i = 0
        while Q:
#            print ">>>>> ITERATION %i <<<<<" % i
            try:
                dist_u, u = pq1.extract_min()
            except pq1.QueueError:
                break
            for v, dist_v in G[u].iteritems():
                alt = self.dist[u-1] + dist_v[0]
                if alt < self.dist[v-1]:
                    self.dist[v-1] = alt
                    self.prev[v-1] = u
                    pq1.decrease_priority(v,alt)
            i += 1
        print ">>> FASTEST PATH ACCORDING TO DIJKSTRA'S ALGORITHM\n"
        if self.prev[self.t-1] == -1: quit("\nTerminating: No shortest path to end node %i" % self.t)
        P = self.trace_path(G)
        T = K_Distinct_Paths.get_shortest_path_tree(self.prev)
        self.process_writing(P, G, min_energies, ts_energies, 1)
        self.marchetti_ktn(G, P, T, min_energies, ts_energies)

    '''Given a flat list of the parent for each node, return a tree'''
    @staticmethod
    def get_shortest_path_tree(parents):
        tree = {i: [] for i in range(1,len(parents)+1)} # each array entry in the dict is a list of children of that node
        for node, parent in enumerate(parents):
            try:
                tree[parent].append(node+1)
            except KeyError:
                if parent == -1: pass # source node has no parent
                else: raise
        return tree

    '''Given a shortest path between start and end nodes, find the rate-limiting edge'''
    @staticmethod
    def find_rle(P, mode="absolute", G=None, min_energies=None, ts_energies=None):
        if mode == "absolute": # return edge with largest weight
            return P.index(max(P,key=lambda x: x[1]))
        elif mode == "barrier": # return edge corresponding to the TS representing the highest individual energy barrier
            P_barr = [[i[0],0.] for i in P]
            for step_no in range(len(P)-1):
                v1, v2 = P[step_no][0], P[step_no+1][0]
                ts_idx = G[v1][v2][1]
                barr_en = ts_energies[ts_idx-1] - min_energies[v1-1]
                P_barr[step_no+1][1] = barr_en
            return P_barr.index(max(P_barr,key=lambda x: x[1]))

    '''Algorithm of Frigioni, Marchetti-Spaccamela & Nanni for determining shortest paths (given an original shortest path tree)
       when an edge on the original shortest path is modified to have increased cost'''
    def marchetti_ktn(self, G, P, sp_tree, min_energies, ts_energies):
        print "\n>>> FINDING THE (k-1) NEXT DISTINCT PATHWAYS BY MARCHETTI-SPACCAMELLA ALGORITHM\n"
        for path_i in range(self.k-1):
            print "\n>>> MARCHETTI ITERATION # %i\n" % (path_i+1)
            rate_lim_edge_idx = K_Distinct_Paths.find_rle(P,"absolute")
#            rate_lim_edge_idx = K_Distinct_Paths.find_rle(P,"barrier",G,min_energies,ts_energies)
            rate_lim_edge = (P[rate_lim_edge_idx-1][0],P[rate_lim_edge_idx][0])
            y, x = rate_lim_edge[0], rate_lim_edge[1]
            print "The rate-limiting edge is:", rate_lim_edge, "with associated cost:", P[rate_lim_edge_idx][1], \
                  "\tThis edge will be blocked."
            G[y][x][0] += self.M # block edge
            G[x][y][0] += self.M # block edge
            pq_m = K_Distinct_Paths.Priority_Queue()
            pq_k = K_Distinct_Paths.Priority_Queue()
            # the priority queue starts with only the owner of the rate-limiting edge
            pq_m.add_with_priority(y,G[y][x][0])
            red_vertices = self.marchetti_colouring(G, pq_m, sp_tree)
            print "\nThe vertices that are coloured red are:\n", red_vertices
            self.marchetti_process_reds(G, pq_k, red_vertices)
            print "\n>>> PATH # %i\n" % (path_i+2)
            P = self.trace_path(G)
            sp_tree = K_Distinct_Paths.get_shortest_path_tree(self.prev)
            self.process_writing(P, G, min_energies, ts_energies, path_i+2)
        if self.write_mdf_all: read_mindata.write_mindatafastestall(self.in_mdf)

    '''Find vertices that are coloured "red"'''
    def marchetti_colouring(self, G, pq_m, sp_tree):
        marchetti_colours = {}
        M = pq_m.pq
        while M:
            nonred_shorter = False
            try:
                dist_z, z = pq_m.extract_min(blocked=True)
            except pq_m.QueueError:
                pass
            for q, dist_q in G[z].iteritems():
                qz_edgecost = G[z][q][0]
                try:
                    if marchetti_colours[q] != 2 and self.dist[q-1] + qz_edgecost == self.dist[z-1]:
                        # there is a nonred neighbour that retains the original shortest path distance to q
                        marchetti_colours[z] = 1 # vertex coloured "pink"
                        self.prev[z-1] = q
                        nonred_shorter = True
                except KeyError:
                    pass
            if not nonred_shorter:
                marchetti_colours[z] = 2 # vertex coloured "red"
                children = sp_tree[z]
                while children: # trace the shortest path tree to find the children of z and add to queue
                    child = children.pop()
                    pq_m.add_with_priority(child,self.dist[child-1])
                    children += sp_tree[child]
        red_vertices = [key for key in marchetti_colours.iterkeys() if marchetti_colours[key]==2]
        return red_vertices

    def marchetti_process_reds(self, G, pq_q, red_vertices):
        for z in red_vertices:
            nonred_nbr = False
            u = None # the best nonred neighbour
            f_level_z_best = float("inf")
            for q, dist_q in G[z].iteritems():
                if q not in red_vertices:
                    nonred_nbr = True
                    qz_edgecost = G[z][q][0]
                    f_level_z = self.dist[q-1] + qz_edgecost
                    if u is None or f_level_z < f_level_z_best: # found a new best nonred neighbour u
                        u = q
                        uz_edgecost = qz_edgecost
                        f_level_z_best = f_level_z
            if not nonred_nbr: # disclude z from shortest path tree
                self.dist[z-1] = float("inf")
                self.prev[z-1] = -1 # null
            else:
                self.dist[z-1] = self.dist[u-1] + uz_edgecost
                self.prev[z-1] = u
                pq_q.add_with_priority(z,self.dist[z-1])
        Q = pq_q.pq
        while Q:
            try:
                dist_z, z = pq_q.extract_min()
            except pq_q.QueueError:
                break
            for h, dist_h in G[z].iteritems():
                if h in red_vertices:
                    hz_edgecost = G[z][h][0]
                    f_level_z = self.dist[z-1] + hz_edgecost
                    if f_level_z < self.dist[h-1]:
                        self.dist[h-1] = f_level_z
                        self.prev[h-1] = z
                        pq_q.decrease_priority(h,f_level_z)
        if self.prev[self.t-1] == -1: quit("\nTerminating: No shortest path to end node %i" % self.t)

    def trace_path(self, G):
        P = [] # shortest path
        P.append((self.t,0.))
        cum_dist = self.dist[self.t-1]
        print "total dist:", cum_dist
        print "node\tcum dist\t\tcurr dist\t\tthis edge\t\t\tcurr edge\t\tcoll dist"
        print self.t, "\t", cum_dist, "\t", 0.
        parent = self.prev[self.t-1]
        old_parent = self.t
        curr_dist, prev_dist = 0., 0.
        coll_dist, curr_edge = 0., 0.
        while parent != -1:
            curr_dist += cum_dist - self.dist[parent-1]
            cum_dist = self.dist[parent-1]
            this_edge = curr_dist - prev_dist
            curr_edge = G[parent][old_parent][0]
            coll_dist += curr_edge
            P.append((parent,this_edge))
            print parent, "\t", cum_dist, "\t", curr_dist, "\t", this_edge, "\t\t", curr_edge, "\t", coll_dist
            old_parent = parent
            parent = self.prev[parent-1]
            prev_dist = curr_dist
        print "%i steps along fastest path" % len(P)
        return P

    def process_writing(self, path, G, min_energies, ts_energies, path_no):
        if self.write_epath: read_mindata.write_epath(path, G, min_energies, ts_energies, path_no)
        if self.write_mdf_all:
            for step in path:
                if self.in_mdf[step[0]] == False: self.in_mdf[step[0]] = True
            if path_no == self.k: read_mindata.write_mindatafastestall(self.in_mdf)

    def calc_edge_costs(self, E_TS, E_min1, E_min2):
        if self.costfunc == "noe_ts": 
            ts_cost1 = np.exp(E_TS/(K_Distinct_Paths.k_B*self.T))
            return ts_cost1, ts_cost1
        elif self.costfunc == "noe_b":
            B_TS = E_TS - np.max(np.array([E_min1,E_min2]))
            ts_cost1 = np.exp(B_TS/(K_Distinct_Paths.k_B*self.T))
#            return ts_cost1, ts_cost1
            return B_TS, B_TS
        elif self.costfunc == "evans":
            B_TS_1, B_TS_2 = E_TS - E_min1, E_TS - E_min2
            ts_cost1 = np.exp(B_TS_1/(K_Distinct_Paths.k_B*self.T))
            ts_cost2 = np.exp(B_TS_2/(K_Distinct_Paths.k_B*self.T))
#            return ts_cost1, ts_cost2
            return B_TS_1, B_TS_2
        else: quit("Invalid selection for cost function")
#        return np.exp(50*E_TS/self.T)

    def build_graph(self, min_energies, ts_energies, ts_conns):
        construct_graph1 = K_Distinct_Paths.Construct_Graph()
        for i in range(1,np.shape(min_energies)[0]+1): construct_graph1.add_vertex(i)
        for i in range(1,np.shape(ts_energies)[0]+1):
            ts_cost1, ts_cost2 = self.calc_edge_costs(ts_energies[i-1],min_energies[ts_conns[i-1,0]-1], \
                                                       min_energies[ts_conns[i-1,1]-1])
            construct_graph1.add_edge(ts_conns[i-1,0],ts_conns[i-1,1],ts_cost1,i)
            construct_graph1.add_edge(ts_conns[i-1,1],ts_conns[i-1,0],ts_cost2,i) # "bidirected" adjacency list
        return construct_graph1.G

    class Construct_Graph(object):

        def __init__(self):
            self.G = {}

        def add_vertex(self, v):
            self.G[v] = {}

        def add_edge(self, v1, v2, e_cost, label=None):
            if label==None: self.G[v1][v2] = e_cost
            else: self.G[v1][v2] = [e_cost, label]

    class Priority_Queue(object):

        def __init__(self):
            self.pq = []
            self.v_log = {} # mapping of elements v to position in priority queue

        def add_with_priority(self, v, priority):
            self.v_log[v] = priority
            heapq.heappush(self.pq,[priority,v])

        def decrease_priority(self, v, new_priority):
            self.add_with_priority(v, new_priority)

        # return highest priority element
        def extract_min(self, blocked=False):
            while self.pq:
#                print "highest priority elem is:", self.pq[0], "lowest priority elem is:", self.pq[-1]
                vpop = heapq.heappop(self.pq)
                try:
                    # encountered for a blocked edge or irrelevant (initialised Dijkstra) entry:
                    if self.v_log[vpop[1]] == float("inf"):
                        # with the blocked flag, entries with cost=inf may not be irrelevant, but may correspond to blocked edges,
                        # which we want to return
                        if blocked: return vpop
                        else: continue
                    if abs(vpop[0] - self.v_log[vpop[1]]) <= 1.0E-06: # the popped entry is the up-to-date value
                        del self.v_log[vpop[1]]
                        return vpop
                except KeyError: # occurs if there is no up-to-date entry for the vertex
                    continue
            raise self.QueueError("During pop attempt: Priority queue is empty")

        class QueueError(Exception):
            pass

if __name__ == "__main__":
    '''
    min_energies, ts_energies, ts_conns = read_mindata.get_data()
    k_distinct_paths1 = K_Distinct_Paths()
    G = k_distinct_paths1.build_graph(min_energies, ts_energies, ts_conns)
    k_distinct_paths1.dijkstra_ktn(G, min_energies, ts_energies)
    '''
    # A TOY TEST PROBLEM
    G = {1: {2: [11,1], 3: [4,2], 7: [6,3], 4: [2,4], 5: [5,5]}, 2: {6: [6,6], 3: [3,7]}, 3: {6: [4,8], 7: [2,9], 4: [1,10]}, \
         4: {5: [1,11], 7: [3,12], 10: [5,13], 8: [5,14], 9: [6,15]}, 5: {9: [4,16]}, 6: {7: [1,17], 12: [5,18], 10: [3,19]}, \
         7: {10: [3,20]}, 8: {12: [7,21], 13: [5,22], 11: [5,23], 9: [3,24]}, 9: {11: [5,25], 14: [5,26]}, \
         10: {12: [2,27], 13: [4,28], 11: [6,29]}, 11: {13: [4,30], 16: [1,31], 14: [2,32]}, 12: {13: [4,33], 15: [5,34]}, \
         13: {15: [2,35], 16: [3,36], 17: [4,37]}, 14: {16: [2,38], 17: [8,39]}, 15: {16: [4,40], 17: [6,41]}, \
         16: {17: [2,42]}, 17: {}}
    # make "bidirected" graph
    for v in G.iterkeys():
        for u in G[v].iterkeys():
            e_cost, ts_idx = G[v][u]
            if u < v: continue
            G[u][v] = [e_cost, ts_idx]
    k_distinct_paths1 = K_Distinct_Paths()
    k_distinct_paths1.dijkstra_ktn(G)
#    '''
