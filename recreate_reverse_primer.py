from tnseq_tags import *

def check_hairpins_illumina(seq, illumina='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'):
    seq = illumina + seq
    trimers = []
    rc_trimers = []
    for offset in range(0,3):
        trimers.extend([seq[x:x+3] for x in range(offset,len(seq),3)])
    rc_trimers.extend(fast_rc(x) for x in trimers if len(x)==3)
    return len(set(trimers).intersection(set(rc_trimers)))


if __name__ == "__main__":
    # Set up parallelisation
    fwd_primers = []
    fwd_spacers = []
    tag_spacers = []
    tag_primers = []
    rev_spacers = []
    with open('results.txt', 'r') as fh:
        fh.readline()
        for line in fh.readlines():
            parts = line.split('\t')
            fwd_primers.append(parts[1])
            fwd_spacers.append(parts[2])
            tag_spacers.append(parts[3])
            tag_primers.append(parts[4])
            rev_spacers.append(parts[5])


    # print(fwd_primers[0])
    # print(tag_spacers[0:10])
    # print(tag_primers[0:10])
    # print(rev_spacers[0:10])q
   # print(check_hairpins2('TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'))
    #print(check_hairpins2('TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTATGAGGAGAGTAGGAGGCAATGG'))

    print(check_hairpins_illumina('', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'))
    print(check_hairpins_illumina('ATCCCTACCCTTATCCCACTTACG', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'))
    print(check_hairpins_illumina('ACTCCCTCTATTGACCTGCTATCC', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'))
    print(check_hairpins_illumina('ACGGATAGATTGAGGAGGGTAAGG', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'))
    print(check_hairpins_illumina('TTCTCCCTCATCGGTCTCTAATCC', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'))


    pool = mp.Pool(32)

    # Generate a Rev Primer
    sys.stdout.write("Generating reverse primers..\n")
    rev_primers = gen_primers(24, 1000000, pool, path='/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/chris/software/tnseq_tags')
    with open("data/reverse_primers.txt", 'w') as fo:
        for primer in rev_primers:
            fo.write(str(primer) + '\n')
    sys.stdout.write("Generated {} reverse primers.\n\n".format(len(rev_primers)))
    sys.stdout.flush()
    rev_seqs = [x[0] for x in rev_primers]
    hairpins = [check_hairpins_illumina(seq) for seq in rev_seqs]
    print(hairpins)
    min_hairpin = min(hairpins)
    rev_primers_best = [rev_primers[i] for i,v  in enumerate(hairpins) if v == min_hairpin]
    print(rev_primers_best)

    #Take the fwd and best rev primers, pair and check
    sys.stdout.write("Determining best outer primer pair..\n")
    fwd_rev_pairs = itertools.product([(fwd_primers[0], 0)], rev_primers_best)
    fwd_rev_pairs = pool.starmap(check_pair_primer3, [(x[0][0], x[1][0]) for x in fwd_rev_pairs])
    fwd_rev_pairs = sorted(list(fwd_rev_pairs), key=lambda x: x[2])
    with open("data/fwd_rev_pairs.txt", 'w') as fo:
        for pair in fwd_rev_pairs:
            fo.write(str(pair) + '\n')

    best_fwd_primer = fwd_rev_pairs[0][0]
    best_rev_primer = fwd_rev_pairs[0][1]

    # Construct and check full sequences
    sys.stdout.write("Checking {} final constructs..\n".format(min([len(tag_spacers), len(tag_primers)])))
    final_constructs = pool.starmap(check_final_construct,
                                    [(best_fwd_primer, fwd_spacers[0], x[0], x[1], rev_spacers[0], best_rev_primer)
                                     for x in zip(tag_spacers, tag_primers)])
    final_constructs = sorted(list(final_constructs), key=lambda x: x[8])
    sys.stdout.write("Finished! Generated {} final constructs.\n\n".format(len(final_constructs)))

    # Output
    with open("data/regen_reverse_results.txt", 'w') as fo:
        fo.write(
            "Full Sequence\tFwd Primer\tFwd Spacer\tTag Spacer\tTag Primer\tRev Spacer\tRev Primer\tOuter Score\tTag Score\n")
        for fc in final_constructs:
            for x in fc:
                fo.write("{}\t".format(x))
            fo.write("\n")

    sys.stdout.write("Completed in {} seconds.\n".format(time.time() - start))
