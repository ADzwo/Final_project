import math, copy
from tuples import *
from block import CollinearBlock, BlockExtensions
from block_tools import mark_vertices_as_used, scoring_function, remove_short_walks#, save_block_to_gff, save_maf

def find_collinear_block(graph, v_idx, PARAM_a, PARAM_b, PARAM_m):
    '''
    Function finds a collinear block in graph, 
    initialized with colliear walk seeds with index v_idx.
    If a block is retrieved, it is returned. Otherwise, None is returned.
    '''
    # Find the seeds --- not used occurrences of vertex with index v_idx.
    collinear_seeds, carrying_seed_orientation = graph.find_seeds(v_idx, PARAM_a=PARAM_a)
    if not collinear_seeds:
        return None
    block_dict = {} # forward and backward block in form {1: forward, -1: backward}
    for forward_backward, seeds in zip([1, -1], [collinear_seeds, [CollinearWalk(*s[:-1], -s.orient) for s in collinear_seeds]]):
        # go in one direction
        best_score = 0
        best_block = None
        new_block = CollinearBlock(v_idx, copy.copy(seeds), carrying_seed_orientation*forward_backward)
        new_score = 0
        while new_score>=0: 
            # 'large' step (find carrying path extension)
            w0_idx = new_block.carrying_path[-1]
            w0_orientation = new_block.carrying_path_orientations[-1]
            Q = BlockExtensions(new_block.collinear_walks, graph, w0_idx, w0_orientation, PARAM_b=PARAM_b) # collinear_walks, graph, w0
            r = Q.get_carrying_path_extension(graph, new_block.collinear_walks, new_block.carrying_path, new_block.carrying_path_orientations) # r - shortest walk from w0 to t (vertex, orientation), where t - index of the most frequently visited vertex;
            if not r:
                break
            for wi in r: 
                # 'small' step (take next vertex of the carrying path extension and extend collinear walks)
                new_block.carrying_path.append(wi.vertex)
                new_block.carrying_path_orientations.append(wi.orientation)
                update_collinear_walks(new_block, wi, Q, graph, PARAM_b=PARAM_b)
                new_score = scoring_function(new_block, graph, PARAM_b=PARAM_b, PARAM_m=PARAM_m)
                new_block.carrying_path_length_so_far += graph.vertices[wi.vertex].v_length
                if math.isinf(new_score):
                    break # stop
                if new_score>best_score: # keep the best state of the block
                    best_block = copy.copy(new_block)
                    best_score = new_score
        if best_score>0:
            block_dict[forward_backward] = copy.copy(best_block)
    
    if len(block_dict)>0:
        final_block = merge_blocks_and_postprocess(block_dict, graph, len(collinear_seeds), PARAM_m=PARAM_m)
        return final_block
    else:
        return None

def merge_match_carrying(len_carrying_b, matches_f=None, matches_b=None):
    '''
    Function merges match_carrying attributes of the forward and backward blocks.
    '''
    if matches_b is None:
        final_matches = matches_f
    elif matches_f is None:
        final_matches = [(-m[0], m[1]) for m in reversed(matches_b)]
    else:
        final_matches = [(-m[0], m[1]) for m in reversed(matches_b[1:])]
        final_matches += matches_f
    return [(m[0]+len_carrying_b-1, m[1]) for m in final_matches]

def merge_forward_backward_blocks(block_dict, nr_seeds):
    '''
    Function merges the forward and backward block into one.
    nr_seeds is the number of walks which initialized the block 
    (the same for both blocks).
    Function returns the merged block.
    '''
    if len(block_dict)>=3:
        raise ValueError(f'Block dict should contain at most two blocks \
                         (forward and backward) Got {len(block_dict)=}.')
    if len(block_dict)==1:
        return list(block_dict.values())[0]
    if len(block_dict)==2:
        forward_block = block_dict[1]
        backward_block = block_dict[-1]
        carrying_b = backward_block.carrying_path
        carrying_orient_b = backward_block.carrying_path_orientations
        final_block = copy.copy(forward_block)
        final_block.match_carrying = []
        # Merge common walks
        for w_idx in range(nr_seeds):
            walk_f = forward_block.collinear_walks[w_idx]
            walk_b = backward_block.collinear_walks[w_idx]
            if walk_f.genome!=walk_b.genome:
                raise ValueError(f'Genome should be the same for the seed walks of forward and backward block.\
                                   Got \n {walk_f=},\n {walk_b=}')
            final_block.collinear_walks[w_idx] = CollinearWalk(walk_f.genome, min(walk_f.start,walk_b.start), 
                                                                 max(walk_f.end,walk_b.end), walk_f.orient)
            matches_f = forward_block.match_carrying[w_idx]
            matches_b = backward_block.match_carrying[w_idx]
            final_block.match_carrying.append(merge_match_carrying(len(carrying_b), matches_f, matches_b))
            
        # Merge the remaining elements of match_carrying
        for matches in forward_block.match_carrying[nr_seeds:]:
            final_block.match_carrying.append(merge_match_carrying(len(carrying_b), matches_f=matches, matches_b=None))
        for matches in backward_block.match_carrying[nr_seeds:]:
            final_block.match_carrying.append(merge_match_carrying(len(carrying_b), matches_f=None, matches_b=matches))
            
        # Merge the remaining walks
        for walk in backward_block.collinear_walks[nr_seeds:]:
            final_block.collinear_walks.append(CollinearWalk(*walk[:-1], -walk.orient))
            
        # merge carrying paths and their orientations
        final_block.carrying_path = list(reversed(carrying_b))[:-1] + forward_block.carrying_path
        final_block.carrying_path_orientations = [-o for o in reversed(carrying_orient_b)][:-1] + forward_block.carrying_path_orientations
        final_block.scores = []
    return final_block

def merge_blocks_and_postprocess(block_dict, graph, nr_seeds, PARAM_m):
    '''
    Function merges the forward and backward blocks, marks their verties as used,
    removes walks which are too short (shorter than PARAM_m)
    and removes matches having the same first element 
    (describing the position on carrying path).
    Other arguments:
    - graph --- the variation graph,
    - nr_seeds --- the number of walks which initialized the block 
      (the same for both blocks).
    '''
    final_block = merge_forward_backward_blocks(block_dict, nr_seeds)
    mark_vertices_as_used(graph, final_block)
    remove_short_walks(final_block, graph, PARAM_m=PARAM_m)
    final_block.remove_doubled_matches()
    return final_block

def find_walk_to_extend(block, block_extensions, g_idx, g_path, o_nr_on_path, carrying_orient):
    '''
    Function finds the collinear walk to be extended 
    till the occurrence of vertex wi from genome g_idx, 
    on position o_nr_on_path of its path (g_path).
    If more then one walk is found, the one which ends further is selected.
    Other arguments:
    - block --- collinear block (Block),
    - block_extensions --- extensions of block (BlockExtensions),
    - carrying_orient --- the orientation of the vertex wi in the carrying path.
    '''
    walk_to_extend = None
    orient_on_extension = g_path[o_nr_on_path].orientation
    for e_idx, extension in block_extensions.extensions.items():
        if extension.genome==g_idx:
            walk = block.collinear_walks[e_idx]
            if walk.start<=o_nr_on_path<=walk.end:
                continue
            if extension.start<=o_nr_on_path<=extension.end and extension.orient==carrying_orient*orient_on_extension:
                if walk_to_extend is None:
                    walk_to_extend = e_idx
                # if there are multiple walks containing this occurrence of w, 
                # select the one ending further
                else:
                    if extension.orient==1:
                        if walk.end>block.collinear_walks[walk_to_extend].end:
                            walk_to_extend = e_idx
                    else:
                        if walk.start<block.collinear_walks[walk_to_extend].start:
                            walk_to_extend = e_idx
    return walk_to_extend

def update_collinear_walks(block, wi_info, block_extensions, graph, PARAM_b): # vertex_info - tuple(vertex index, orientation)
    '''
    Update collinear walks.
    For each occurrence of wi (described by wi_info --- tuple of its index and orientation),
    the function searches for a walk which can be extended till this occurrence.
    If such a walk is found, it is extended till the occurrence of wi.
    Otherwise, the occurrence of wi becomes a seed of a new collinear walk.
    Scores of all walks are updated.
    '''
    wi = graph.vertices[wi_info.vertex]
    walks_updated_score = set()
    for occurrence in wi.occurrences:
        # 1) Search for walks, whose extensions contain the occurrence of w
        g_idx = occurrence.genome
        g_path = graph.genomes[g_idx].path
        o_nr_on_path = occurrence.nr_on_path
        if g_path[o_nr_on_path].used==True:
            continue
        # 2) Find a walk to extend, if possible.
        walk_to_extend = find_walk_to_extend(block, block_extensions, g_idx, g_path, o_nr_on_path, wi_info.orientation)
        # 3a) If the path is found, extend it till w AND update the extension
        if walk_to_extend is not None:
            # if walk_to_extend==-1:
            #     continue
            walk = block.collinear_walks[walk_to_extend]
            if walk.orient==1:
                block.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, walk.start, o_nr_on_path, walk.orient)
            else:
                block.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, o_nr_on_path, walk.end, walk.orient)
            block_extensions.update_extension(block.collinear_walks[walk_to_extend], graph, walk_idx=walk_to_extend, PARAM_b=PARAM_b)
            old_score = block.scores[walk_to_extend]
            block.scores[walk_to_extend] = Score(old_score.q1, 0)
            walks_updated_score.add(walk_to_extend)

            # Update match_carrying
            block.match_carrying[walk_to_extend].append((len(block.carrying_path)-1, o_nr_on_path))
            
        # 3b) if such a walk is not found, occurrence becomes a new collinear path (provided it is not used)
        else:
            assert g_path[o_nr_on_path].vertex in block.carrying_path
            occurrence_orientation = g_path[o_nr_on_path].orientation
            orient = wi_info.orientation * occurrence_orientation
            block.collinear_walks.append(CollinearWalk(g_idx, o_nr_on_path, o_nr_on_path, orient))
            # Create extension for the new collinear walk
            block_extensions.update_extension(block.collinear_walks[-1], graph, walk_idx=len(block.collinear_walks)-1, PARAM_b=PARAM_b)
            block.scores.append(Score(block.carrying_path_length_so_far, 0))
            walks_updated_score.add(len(block.collinear_walks)-1)
            block.match_carrying.append([(len(block.carrying_path)-1, o_nr_on_path)])
            
    # Update scores
    for w_idx in range(len(block.collinear_walks)):
        if w_idx not in walks_updated_score:
            old_score = block.scores[w_idx]
            block.scores[w_idx] = Score(old_score.q1, old_score.q3+wi.v_length)


