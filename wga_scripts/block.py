from tuples import *
from block_tools import walk_length, find_vertex_on_path_till_b

class CollinearBlock:
    carrying_path: list # vertex index: int
    carrying_path_orientations: list # orientation: int
    collinear_walks: list # genome index: int,
                            # start index from genome's path: int,
                            # end index from genome's path: int,
                            # orientaion with respect to the path: int
    scores: list
    carrying_path_length_so_far: int
    match_carrying: list # a list of the same length as block.collinear_walks, 
                        # containing lists of tuples 
                        # (position on carrying path - from 0 to len(carrying_path)-1, 
                        #  position on walk.genome - between walk.start and walk.end)
    def __init__(self, seed_idx, seed_occurrences, carrying_seed_orientation):
        self.carrying_path = [seed_idx] # vertex index
        self.carrying_path_orientations = [carrying_seed_orientation] # orientation
        self.collinear_walks = seed_occurrences # CollinearWalks of length 1; carrying_seed included
        self.scores = [Score(0, 0) for w in self.collinear_walks] # ['q1', 'q3']
        self.carrying_path_length_so_far = 0
        self.match_carrying = [[(0, w.start)] for w in self.collinear_walks]

    
    def remove_doubled_matches(self):
        '''
        Function removes from self.match_carrying the elements
        having the same first elements (describing the position on carrying path).
        The first and last elements are not removed.
        So, if the 1st and 2nd elements have the same first element, 
        the 2nd one is removed.
        For other elements, the one with smaller index is removed.
        '''
        for w_idx in range(len(self.collinear_walks)):
            matches = self.match_carrying[w_idx].copy()
            for i in range(len(matches)-1, 0, -1):
                prev_match_carrying, _ = self.match_carrying[w_idx][i-1]
                match_carrying, _ = self.match_carrying[w_idx][i]
                if prev_match_carrying==match_carrying:
                    if i==1: # we need a match at the beginning of the walk
                        self.match_carrying[w_idx].pop(i)
                    else:
                        self.match_carrying[w_idx].pop(i-1)
                
class BlockExtensions:
    extensions: dict # list of extensions of type CollinearWalk
    coverage: dict # {vertex index: coverage}, where coverage - number of occurrences in extensions
    shortest_walk: dict  # {vertex index: (distance, *shortest_walk)}
                    # distance - length of the shortest walk from w_0
                    # shortest_walk (walk number, start, end)
    def __init__(self, collinear_walks, graph, w0_idx, w0_orientation, PARAM_b):
        self.extensions = {}
        self.coverage = {}
        self.shortest_walk = {}
        
        for walk_nr, walk in enumerate(collinear_walks):
            g_idx = walk.genome
            g_path = graph.genomes[g_idx].path
            g_len = len(g_path)
            if walk.orient==1 and walk.end==g_len-1: # if we've reached the end of the genome
                continue
            elif walk.orient==-1 and walk.start==0:
                continue
            
            if walk.orient==1:
                to_search_from = walk.end # we start searching at the most distal position of the walk
                to_end_search = g_len - 1
            else:
                to_search_from = walk.start
                to_end_search = 0
            proximal = to_search_from + walk.orient # proximal - the end of extension proximal to the collinear walk
            
            if g_path[proximal].used==True:
                break
            w0_nr_on_path = find_vertex_on_path_till_b(graph, walk, 
                                                       w0_idx, w0_orientation,
                                                       proximal=to_search_from,
                                                       distal=to_end_search,
                                                       PARAM_b=PARAM_b) # taking the closest position to the collinear walk <--- maybe should take the last one?
            p_length = 0
            i = proximal
            if w0_nr_on_path is None: # only increase coverage
                while True:
                    if g_path[i].used==True:
                        i -= walk.orient
                        break
                    v_idx = g_path[i].vertex
                    v_oriented = v_idx*g_path[i].orientation
                    self.coverage[v_oriented] = self.coverage.get(v_oriented, 0) + 1
                    p_length += graph.vertices[v_idx].v_length
                    if p_length>=PARAM_b:
                        break
                    if i in {0, g_len-1}:
                        break
                    i += walk.orient
            else:
                while True:
                    if g_path[i].used==True:
                        i -= walk.orient
                        break
                    v_idx = g_path[i].vertex
                    v_oriented = v_idx*g_path[i].orientation
                    self.coverage[v_oriented] = self.coverage.get(v_oriented, 0) + 1
                    p_length += graph.vertices[v_idx].v_length
                    is_w0_before_v = ((i-w0_nr_on_path)*walk.orient) > 0
                    if is_w0_before_v:
                        distance = walk_length(CollinearWalk(g_idx, min(w0_nr_on_path,i), max(w0_nr_on_path,i), 1), graph)
                        if v_oriented not in self.shortest_walk or self.shortest_walk[v_oriented].distance>distance:
                            self.shortest_walk[v_oriented] = PathFromW0(distance, walk_nr, w0_nr_on_path, i)
                            assert (i < walk.start and i < walk.end) or (i > walk.start and i > walk.end)
                            
                    if p_length>=PARAM_b:
                        break
                    if i in {0, g_len-1}:
                        break
                    i += walk.orient # going forwards or backwards on the genome
            self.extensions[walk_nr] = CollinearWalk(g_idx, min(proximal, i), max(proximal, i), walk.orient)
            
    
    def get_carrying_path_extension(self, graph, collinear_walks, carrying_path=[], carrying_path_orientations=[]):
        '''
        Function finds 
         - index t of a vertex with field distance > 0 (reachable from w0) via a genomic walk, visited by the most extensions.
        Function returns 
         - the shortest walk r from w0 to t.
        '''
        highest_coverage = -1
        w0_to_t = None
        for v_oriented in self.shortest_walk: # for all vertices reachable from w0 within distance PARAM_b
            if self.coverage[v_oriented]>highest_coverage: # if coverage is greater than the highest one by now
                w0_to_t = self.shortest_walk[v_oriented]
                highest_coverage = self.coverage[v_oriented]
        if w0_to_t is None:
            return []
        
        r = []
        walk = collinear_walks[w0_to_t.walk_nr] 
        genome = graph.genomes[walk.genome]
        for i in range(w0_to_t.w0_nr_on_path+walk.orient, w0_to_t.t_nr_on_path+walk.orient, walk.orient):
            g_path_pos = genome.path[i]
            r.append(CarryingPathExtension(g_path_pos.vertex, g_path_pos.orientation*walk.orient))
        return r

    def update_extension(self, walk, graph, walk_idx, PARAM_b):
        '''
        Function updates self.extensions[walk_idx] based on 
        walk (collinear walk with index walk_idx), 
        graph (variation graph) and PARAM_b (maximal bubble length).
        '''
        g_idx = walk.genome
        genome = graph.genomes[g_idx]
        g_len = len(genome.path)
        proximal = walk.end + walk.orient
        if proximal>=g_len or proximal<0 or genome.path[proximal].used==True:
            if walk_idx in self.extensions:
                del self.extensions[walk_idx]
        else:
            distal = proximal
            ext_length = 0
            while True:
                g_path_pos = genome.path[distal]
                if g_path_pos.used==True:
                    distal -= walk.orient
                    break
                ext_length += graph.vertices[g_path_pos.vertex].v_length
                if ext_length>=PARAM_b:
                    break
                distal += walk.orient
                if distal<0 or distal>=g_len:
                    distal -= walk.orient
                    break
            self.extensions[walk_idx] = CollinearWalk(g_idx, min(proximal,distal), max(proximal,distal), walk.orient)