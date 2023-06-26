import math, random, copy
from tuples import *
from config import *
from block import CollinearBlock, BlockExtensions
from block_tools import mark_vertices_as_used, scoring_function, remove_short_walks

def find_collinear_blocks(graph):
    collinear_blocks = []

    vertex_indices_order = list(range(len(graph.vertices)))
    if SORT_SEEDS=='nr_occurrences':
        vertex_indices_order = sorted(vertex_indices_order, key=lambda x: len(graph.vertices[x].occurrences), reverse=True)
    elif SORT_SEEDS=='length':
        vertex_indices_order = sorted(vertex_indices_order, key=lambda x: graph.vertices[x].v_length, reverse=True)
    else:
        random.shuffle(vertex_indices_order)
    
    for v_idx in vertex_indices_order: # select a vertex --- seed of a new CollinearBlock
        collinear_seeds, carrying_seed_orientation = graph.find_seeds(v_idx)
        if not collinear_seeds:
            continue
        
        blocks_forward_backward = {}
        for forward_backward, seeds in zip([1, -1], [collinear_seeds, [CollinearWalk(*s[:-1], -s.orient) for s in collinear_seeds]]):
            best_score = -1
            best_block = None
            new_block = CollinearBlock(v_idx, copy.copy(seeds), carrying_seed_orientation*forward_backward)
            for walk in new_block.collinear_walks:
                assert graph.genomes[walk.genome].path[walk.start].vertex in new_block.carrying_path, f'{graph.genomes[walk.genome].path[walk.start].vertex=},\n{new_block.carrying_path=}'
                assert graph.genomes[walk.genome].path[walk.end].vertex in new_block.carrying_path, f'{graph.genomes[walk.genome].path[walk.start].vertex=},\n{graph.genomes[walk.genome].path[walk.end].vertex=},\n{new_block.carrying_path=}'
            new_score = 0
            while new_score>=0:
                assert len(new_block.match_carrying) == len(new_block.collinear_walks)
                w0_idx = new_block.carrying_path[-1]
                Q = BlockExtensions(new_block.collinear_walks, graph, w0_idx) # collinear_walks, graph, w0
                r = Q.get_carrying_path_extension(graph, new_block.collinear_walks) # r - shortest walk from w0 to t (vertex, orientation), where t - index of the most frequently visited vertex;
                if not r:
                    break
                for wi in r:
                    new_block.carrying_path.append(wi.vertex)
                    new_block.carrying_path_orientations.append(wi.orientation)
                    update_collinear_walks(new_block, wi, Q, graph)
                    new_score = scoring_function(new_block, graph)
                    new_block.carrying_path_length_so_far += graph.vertices[wi.vertex].v_length
                    if math.isinf(new_score):
                        break
                    if new_score>best_score:
                        best_block = copy.copy(new_block)
                        best_score = new_score 
                for s_nr, seed in enumerate(seeds):
                    assert seed.genome==new_block.collinear_walks[s_nr].genome
                # check_walk_orientations(graph, new_block)
            if best_score>0:
                blocks_forward_backward[forward_backward] = copy.copy(best_block)
        
        if blocks_forward_backward:
            final_block = merge_blocks_and_postprocess(blocks_forward_backward, graph, len(collinear_seeds))
            # check_walk_orientations(graph, final_block)
            collinear_blocks.append(final_block)
            
    return collinear_blocks


def merge_forward_backward_blocks(block_dict, nr_seeds):
    assert len(block_dict)<3, 'Block list should contain at most two blocks (forward and backward)'
    if len(block_dict)==1:
        return list(block_dict.values())[0]
    if len(block_dict)==2:
        forward_block = block_dict[1]
        backward_block = block_dict[-1]
        carrying_b = backward_block.carrying_path
        carrying_orient_b = backward_block.carrying_path_orientations

        assert len(forward_block.match_carrying)==len(forward_block.collinear_walks)
        assert len(backward_block.match_carrying)==len(backward_block.collinear_walks)
        
        # for walk, matches in zip(forward_block.collinear_walks, forward_block.match_carrying):
        #     match_carrying_check(walk, matches)
        # for walk, matches in zip(backward_block.collinear_walks, backward_block.match_carrying):
        #     match_carrying_check(walk, matches)
        # Merge common walks
        # TO FIX - merge walks which form a bubble near seeds
        for w_idx in range(nr_seeds):
            walk_f = forward_block.collinear_walks[w_idx]
            walk_b = backward_block.collinear_walks[w_idx]
            assert walk_f.genome==walk_b.genome, f'Genome should be the same for the seed walks of forward and backward block, got \n {walk_f=},\n {walk_b=}'
            forward_block.collinear_walks[w_idx] = CollinearWalk(walk_f.genome, min(walk_f.start,walk_b.start), 
                                                                 max(walk_f.end,walk_b.end), walk_f.orient)
            matches_f = forward_block.match_carrying[w_idx]
            matches_b = backward_block.match_carrying[w_idx]
            assert -list(reversed(matches_b))[-1][0]==matches_f[0][0]
            assert list(reversed(matches_b))[-1][1]==matches_f[0][1]
            forward_block.match_carrying[w_idx] = ([(-m[0]+len(carrying_b)-1, m[1]) for m in list(reversed(matches_b))[:-1]] 
                                                 + [(m[0]+len(carrying_b)-1, m[1]) for m in matches_f])
            # match_carrying_check(forward_block.collinear_walks[w_idx], forward_block.match_carrying[w_idx])
            
        # Merge the remaining elements of match_carrying
        for w_idx, matches in enumerate(forward_block.match_carrying[nr_seeds:]):
            forward_block.match_carrying[w_idx+nr_seeds] = [(m[0]+len(carrying_b)-1, m[1]) for m in matches]
        for matches in backward_block.match_carrying[nr_seeds:]:
            forward_block.match_carrying.append([(-m[0]+len(carrying_b)-1, m[1]) for m in list(reversed(matches))])

        # Merge the remaining walks
        for walk in backward_block.collinear_walks[nr_seeds:]:
            forward_block.collinear_walks.append(CollinearWalk(*walk[:-1], -walk.orient))

        assert len(forward_block.match_carrying)==len(forward_block.collinear_walks)
        # for walk, matches in zip(forward_block.collinear_walks, forward_block.match_carrying):
        #     match_carrying_check(walk, matches)
            
        # merge carrying paths and their orientations
        assert forward_block.carrying_path[0]==carrying_b[0]
        assert forward_block.carrying_path_orientations[0]==-carrying_orient_b[0], f'{forward_block.carrying_path_orientations[0]=}, {carrying_orient_b[0]=}'
        carrying_len = len(forward_block.carrying_path) + len(carrying_orient_b) - 1
        forward_block.carrying_path = list(reversed(carrying_b))[:-1] + forward_block.carrying_path
        forward_block.carrying_path_orientations = [-o for o in reversed(carrying_orient_b)][:-1] + forward_block.carrying_path_orientations
        assert len(forward_block.carrying_path)==len(forward_block.carrying_path_orientations)==carrying_len
    return forward_block

def merge_blocks_and_postprocess(blocks_forward_backward, graph, nr_seeds):
    final_block = merge_forward_backward_blocks(blocks_forward_backward, nr_seeds)
    mark_vertices_as_used(graph, final_block)
    remove_short_walks(final_block, graph)
    final_block.remove_doubled_matches()
    for matches in final_block.match_carrying:
        for i in range(len(matches)-1):
            assert matches[i][0]<matches[i+1][0], f'{i=}, {matches[i]=}, {matches[i+1]=}'
    return final_block

def find_walk_to_extend(block, block_extensions, g_idx, g_path, o_nr_on_path, carrying_orient):
    walk_to_extend = None
    orient_on_extension = g_path[o_nr_on_path].orientation
    for e_idx, extension in block_extensions.extensions.items():
        if extension.genome==g_idx and extension.start<=o_nr_on_path<=extension.end:
            if extension.orient==carrying_orient*orient_on_extension:
                if walk_to_extend is None:
                    walk_to_extend = e_idx
                # if there are multiple walks containing this occurrence of w, 
                # select the one ending further
                else:
                    if extension.orient==1:
                        if block.collinear_walks[e_idx].end>block.collinear_walks[walk_to_extend].end:
                            walk_to_extend = e_idx
                    else:
                        if block.collinear_walks[e_idx].start<block.collinear_walks[walk_to_extend].start:
                            walk_to_extend = e_idx
    return walk_to_extend

def update_collinear_walks(block, wi_info, block_extensions, graph): # vertex_info - tuple(vertex index, orientation)
    wi = graph.vertices[wi_info.vertex]
    walks_updated_score = set()
    for occurrence in wi.occurrences:
        # 1) Search for walks, whose extensions contain the occurrence of w
        g_idx = occurrence.genome
        g_path = graph.genomes[g_idx].path
        g_len = len(g_path)
        o_nr_on_path = occurrence.nr_on_path
        if g_path[o_nr_on_path].used==True:
            continue
        # 2) Find a walk to extend, if possible.
        walk_to_extend = find_walk_to_extend(block, block_extensions, g_idx, g_path, o_nr_on_path, wi_info.orientation)
        # 3a) If the path is found, extend it till w AND update the extension
        if walk_to_extend is not None:
            walk = block.collinear_walks[walk_to_extend]
            assert g_path[walk.end].vertex in block.carrying_path
            assert g_path[walk.start].vertex in block.carrying_path
            assert g_path[o_nr_on_path].vertex in block.carrying_path
            if walk.orient==1:
                block.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, walk.start, o_nr_on_path, walk.orient)
            else:
                block.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, o_nr_on_path, walk.end, walk.orient)
            # walk_start_end_check(block.collinear_walks[walk_to_extend], g_len)
            block_extensions.update_extension(block.collinear_walks[walk_to_extend], graph, walk_idx=walk_to_extend)
            old_score = block.scores[walk_to_extend]
            block.scores[walk_to_extend] = Score(old_score.q1, 0)
            walks_updated_score.add(walk_to_extend)

            # Update match_carrying
            block.match_carrying[walk_to_extend].append((len(block.carrying_path)-1, o_nr_on_path))
            # match_carrying_check(block.collinear_walks[walk_to_extend], block.match_carrying[walk_to_extend])
            # check_walk_orientations(graph, block)

        # 3b) if such a walk is not found, occurrence becomes a new collinear path (provided it is not used)
        else:
            assert g_path[o_nr_on_path].vertex in block.carrying_path
            occurrence_orientation = g_path[o_nr_on_path].orientation
            # orient = wi_orient_on_carrying_path * occurrence_orientation
            orient = wi_info.orientation * occurrence_orientation
            block.collinear_walks.append(CollinearWalk(g_idx, o_nr_on_path, o_nr_on_path, orient))
            # walk_start_end_check(block.collinear_walks[-1], g_len)
            # Create extension for the new collinear walk
            block_extensions.update_extension(block.collinear_walks[-1], graph, walk_idx=len(block.collinear_walks)-1)
            block.scores.append(Score(block.carrying_path_length_so_far, 0))
            walks_updated_score.add(len(block.collinear_walks)-1)
            block.match_carrying.append([(len(block.carrying_path)-1, o_nr_on_path)])
            # match_carrying_check(block.collinear_walks[-1], block.match_carrying[-1])
            # check_walk_orientations(graph, block)

    # Update scores
    for w_idx in range(len(block.collinear_walks)):
        if w_idx not in walks_updated_score:
            old_score = block.scores[w_idx]
            block.scores[w_idx] = Score(old_score.q1, old_score.q3+wi.v_length)


