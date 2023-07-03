#!/usr/bin/env python
from __future__ import print_function
import poagraph
import seqgraphalignment
from collections import namedtuple
import os

complement_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
AlignmentTuple = namedtuple('AlignmentTuple', ['sequence', 'stringidxs', 'nodeidxs'])

def reverse_complement(seq):
    seq_complement = ''
    for s in seq[::-1]:
        seq_complement += complement_dict[s]
    return seq_complement

def find_edges_between(graph, start, end, walks_to_poa):
    node = graph.nodedict[start]
    walks_to_align = {w_idx for e in node.outEdges for w_idx in node.outEdges[e].labels}
    node = graph.nodedict[end]
    walks_to_align = walks_to_align.intersection({w_idx for e in node.inEdges for w_idx in node.inEdges[e].labels})
    # let's find the node IDs of the walks to align
    to_use = set(range(start+1, end))
    for w_idx in walks_to_align:
        found_start = False
        if w_idx==-1:
            continue
        for nid in walks_to_poa[w_idx]:
            if nid==end:
                break
            if nid==start:
                found_start = True
            elif found_start:
                to_use.add(nid)
        assert found_start or start==end
    return to_use                        

def carrying_to_poa_idx(block, sequences):
    carrying_to_poa = {}
    carrying_seq = ''
    poa_v_nr = 0
    for v_nr, (v_idx, v_orientation) in enumerate(zip(block.carrying_path, block.carrying_path_orientations)):
        carrying_seq += read_sequence(sequences, v_idx, v_orientation)
        carrying_to_poa[v_nr] = [poa_v_nr] # start of v
        poa_v_nr = len(carrying_seq)
        carrying_to_poa[v_nr].append(poa_v_nr-1) # end of v
    return carrying_to_poa, carrying_seq

def print_alignment(graph):
    alignments = graph.generateAlignmentStrings()
    for label, alignstring in alignments:
        print("{0:15s} {1:s}".format(str(label), alignstring))
    return alignments

def read_sequence(sequences, v_idx, v_orientation, walk_orient=None):
    if walk_orient is None:
        orientation = v_orientation
    else:
        orientation = walk_orient*v_orientation
    
    if orientation==1:
        seq = sequences[v_idx].strip()
    else:
        seq = reverse_complement(sequences[v_idx].strip())
    return seq

def align_matching_fragment(start_nid, end_nid, sequence, seq_len_so_far, old_alignment):
    if old_alignment is None:
        old_alignment = [sequence, 
                         list(range(len(sequence))), 
                         list(range(start_nid, end_nid+1))]
    else:
        old_alignment[0] += sequence
        old_alignment[1] += [sidx+seq_len_so_far for sidx in range(len(sequence))]
        old_alignment[2] += list(range(start_nid, end_nid+1))
        assert len(old_alignment[0])-1==old_alignment[1][-1], f'{len(old_alignment[0])-1=}, {old_alignment[1][-1]}'
    return old_alignment


def add_to_alignment(old_alignment, new_alignment, seq_len_so_far):
    old_alignment[0] += new_alignment.sequence
    old_alignment[1] += [None if sidx is None else sidx+seq_len_so_far for sidx in new_alignment.stringidxs]
    old_alignment[2] += new_alignment.nodeidxs
    return old_alignment

def poa_align(block, var_graph, sequences, _match, _mismatch, _gap, html=True, globalAlign=True, simple=True):
    carrying_to_poa, carrying_seq = carrying_to_poa_idx(block, sequences)
    graph = poagraph.POAGraph(carrying_seq, -1) # carrying path
    walks_to_poa = {}
    entire_sequences = []

    for w_idx in range(len(block.collinear_walks)):
        walk = block.collinear_walks[w_idx]
        matches = block.match_carrying[w_idx]

        to_start = walk.start if walk.orient==1 else walk.end # end and beginning of each walk
        to_end = walk.end+1 if walk.orient==1 else walk.start-1
        g_path = var_graph.genomes[walk.genome].path

        if len(matches)==1:
            match_carrying = matches[0][0]
            start_nid, end_nid = carrying_to_poa[match_carrying]
            entire_sequence = ''
            for i_on_g_path in range(to_start, to_end, walk.orient):
                v_idx = g_path[i_on_g_path].vertex
                v_orientation = g_path[i_on_g_path].orientation
                entire_sequence += read_sequence(sequences, v_idx, v_orientation, walk.orient)
            old_alignment = align_matching_fragment(start_nid, end_nid, entire_sequence, 0, None)
            entire_sequences.append(entire_sequence)
            old_alignment = AlignmentTuple(*old_alignment)
            walks_to_poa[w_idx] = graph.incorporateSeqAlignment(old_alignment, entire_sequence, label=w_idx)
            assert len(walks_to_poa[w_idx])>0
            continue
        
        first = True
        old_alignment = None
        i_on_match = 0
        sequence = ''
        entire_sequence = ''
        i_on_g_path = to_start
        while i_on_g_path!=to_end: # to_end does not belong to walk
            v_idx = g_path[i_on_g_path].vertex
            v_orientation = g_path[i_on_g_path].orientation
            assert i_on_match<len(matches), f'{i_on_match=}, {len(matches)=}, {matches=}, {walk=}'
            if matches[i_on_match][1]!=i_on_g_path: # there is no match on this position of walk
                sequence += read_sequence(sequences, v_idx, v_orientation, walk.orient)
                i_on_g_path += walk.orient
            elif sequence=='': # there is a match and we do not need to align anything
                match_carrying = matches[i_on_match][0]
                start_nid, end_nid = carrying_to_poa[match_carrying]
                assert start_nid<=end_nid
                assert v_orientation*walk.orient==block.carrying_path_orientations[match_carrying]
                sequence = read_sequence(sequences, v_idx, v_orientation, walk.orient)
                assert sequence==carrying_seq[start_nid:end_nid+1], f'{sequence=}, \
                    \n{sequences[v_idx].strip()=}, \n{carrying_seq[start_nid:end_nid+1]=}, \
                    \n{carrying_seq=}, \n{walk.orient*v_orientation=}, {matches[i_on_match][0]=}'

                old_alignment = align_matching_fragment(start_nid, end_nid, sequence, len(entire_sequence), old_alignment)
                i_on_match += 1
                i_on_g_path += walk.orient
                entire_sequence += sequence
                sequence = ''
                first = False
            else:
                assert first == False
                last_match = matches[i_on_match-1]
                match = matches[i_on_match]
                assert last_match[0]<=match[0], f'{last_match=}, {match=}'
                
                start_nid = carrying_to_poa[last_match[0]][1] # the last node of last_match
                end_nid = carrying_to_poa[match[0]][0] # the first node of match
                assert start_nid<=end_nid, f'{last_match[0]=}, {match[0]=}, {start_nid=}, {end_nid=}'
                assert last_match[0]<match[0], f'last_match[0] should be smaller than match[0]. Got {last_match=}, {match=}.'
                
                to_use = find_edges_between(graph, start_nid, end_nid, walks_to_poa)
                
                if to_use:
                    alignment = seqgraphalignment.SeqGraphAlignment(sequence, 
                                                                    graph, fastMethod=not simple,
                                                                    globalAlign=globalAlign,
                                                                    matchscore=_match, mismatchscore=_mismatch,
                                                                    gapscore=_gap, to_use=to_use, end=end_nid)
                else:
                    alignment = AlignmentTuple(sequence, range(len(sequence)), [None for _ in range(len(sequence))])
                old_alignment = add_to_alignment(old_alignment, alignment, len(entire_sequence))
                
                entire_sequence += sequence
                sequence = ''
            
        entire_sequences.append(entire_sequence)

        old_alignment = AlignmentTuple(*old_alignment)
        walks_to_poa[w_idx] = graph.incorporateSeqAlignment(old_alignment, entire_sequence, label=w_idx)
        assert len(walks_to_poa[w_idx])>0

    alignments = graph.generateAlignmentStrings()
    return alignments

