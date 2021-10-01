from needleman_wunsch import needleman_wunsch

if __name__ == "__main__":
    # define sequences to align
    sequenceA = 'TTGCT'
    sequenceB = 'CTTCCT'
    # define match/mismatch/gap penalties
    match = 1
    mismatch = gap = -1
    # instantiate nw class instance
    nw = needleman_wunsch(sequenceA, sequenceB, match, mismatch, gap)
    # build alignment table and get optimal score
    alignment_table, alignment_score = nw.build_table()
    # retrace path through table to find alignment
    alignment_path = nw.trace_back(alignment_table)

    print('Optimal alignment: \n', alignment_path)
    print('Optimal alignment score', alignment_score)
    print('Alignment table: ', alignment_table)
