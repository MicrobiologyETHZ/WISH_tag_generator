import itertools
import multiprocessing as mp
import random
import re
import subprocess
import sys
import time
from Bio.SeqUtils import MeltingTemp

# structure of proposed molecule
# | Fwd Primer | Spacer | Tag start | Tag primer | Spacer | Rev Primer | 
# <----24bp--->|<-24bp->|<---16bp-->|<--(24bp)-->|<-8bp-->|<---24bp--->|

start = time.time()

rc_dict = {'A':'T','C':'G','G':'C','T':'A'}

def fast_rc(seq):
    seq_rc = ''.join([rc_dict[x] for x in seq])
    return(seq_rc[::-1])

def gen_seq(l):
    seq = ''.join(random.choices("ACGT",k=l))
    return(seq)

def gen_primer(l):
    seq = ''.join([''.join(random.choices("ACGT",k=l-5)),''.join(random.choices("AT",k=3)),''.join(random.choices("CG",k=2))])
    return(seq)

def check_gc(seq):
    gc = (seq.count('G')+seq.count('C'))/len(seq)
    return(0.45<=gc<=0.55)

def check_runs(seq):
    # Check single runs
    runs = []
    count = 1
    a = seq[0]
    for b in seq[1:]:
        if b==a:
            count+=1
        else:
            runs.append(count)
            count=1
        a = b
    runs.append(count)

    # Check pair runs even then odd
    count = 1
    a = seq[0:2]
    for b in [seq[x:x+2] for x in range(2,len(seq),2)]:
        if b==a:
            count+=1
        else:
            runs.append(count)
            count=1
        a = b
    runs.append(count)
    a = seq[1:3]
    for b in [seq[x:x+2] for x in range(3,len(seq),2)]:
        if b==a:
            count+=1
        else:
            runs.append(count)
            count=1
        a = b
    runs.append(count)

    return(not any(c>3 for c in runs))

def check_tm(seq):
    tm = MeltingTemp.Tm_NN(seq)
    return(59.5<=tm<=62.5)

def check_palindrome(seq):
    return(fast_rc(seq)!=seq)

def check_hairpins(seq):
    trimers = []
    rc_trimers = []
    for offset in range(0,3):
        trimers.extend([seq[x:x+3] for x in range(offset,len(seq),3)])
    rc_trimers.extend(fast_rc(x) for x in trimers if len(x)==3)
    if len(set(trimers).intersection(set(rc_trimers))) > 0:
        return(False)
    else:
        return(True)

def run_bwa(seqs, path='.'):
    # Output seqs to a fasta file for alignment
    with open("temp.fasta",'w') as fo:
        for seq in seqs:
            fo.write(">{}\n{}\n".format(seq,seq))

    # Align seqs fasta file against a target file with BWA, allowing 2 errors
    sys.stdout.write("Aligning primers..\n")
    aln = subprocess.check_output(f"bwa aln -t 32 -n 2 {path}/db/db.fasta temp.fasta 2> /dev/null | bwa samse {path}/db/db.fasta - temp.fasta 2> /dev/null",shell=True,executable="/bin/bash").decode()
    aln = aln.strip().split("\n")
    aln = [x.strip().split("\t") for x in aln if x[0]!='@']
    aln = {x[0]:x for x in aln}

    return(aln)

def check_misprime(seq,aln):
    return(aln[seq][1] == '4')

def check_primer_primer3(seq):
    argument = "SEQUENCE_ID=primer_pair\\nPRIMER_TASK=check_primers\\nSEQUENCE_PRIMER={}\\nPRIMER_MIN_TM=59.5\\nPRIMER_MAX_TM=62.5\\nPRIMER_OPT_TM=61\\nPRIMER_OPT_SIZE=24\\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(seq)
    output = subprocess.check_output("primer3_core <(printf \"{}\")".format(argument),shell=True,executable="/bin/bash").decode()
    primer_penalty = 100
    if not re.findall('PRIMER_LEFT_NUM_RETURNED=0',output):
        lines = output.split("\n")[1:-2]
        param = {}
        for line in lines:
            parts = line.split("=")
            param[parts[0]] = parts[1]
        primer_penalty = float(param['PRIMER_LEFT_0_PENALTY'])
    return((seq,primer_penalty))

def check_pair_primer3(fwd,rev):
    argument = "SEQUENCE_ID=primer_pair\\nPRIMER_TASK=check_primers\\nSEQUENCE_PRIMER={}\\nSEQUENCE_PRIMER_REVCOMP={}\\nPRIMER_MIN_TM=59.5\\nPRIMER_MAX_TM=62.5\\nPRIMER_OPT_TM=61\\nPRIMER_OPT_SIZE=24\\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(fwd,rev)
    output = subprocess.check_output("primer3_core <(printf \"{}\")".format(argument),shell=True,executable="/bin/bash").decode()
    pair_penalty = 100
    if not re.findall('PRIMER_PAIR_NUM_RETURNED=0',output):
        lines = output.split("\n")[1:-2]
        param = {}
        for line in lines:
            parts = line.split("=")
            param[parts[0]] = parts[1]
        pair_penalty = float(param['PRIMER_PAIR_0_PENALTY'])
    return((fwd,rev,pair_penalty))

def check_illumina_primers(fwd, rev,
                           illumina_fwd="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
                           illumina_rev="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"):
    pcr1_fwd = illumina_fwd + fwd
    pcr1_rev = illumina_rev + rev
    return check_pair_primer3(pcr1_fwd, pcr1_rev)


def check_final_construct(fwd_primer, fwd_spacer, tag_spacer, tag_primer, rev_spacer, rev_primer):
    construct = ''.join([fwd_primer, fwd_spacer, tag_spacer, fast_rc(tag_primer), rev_spacer, fast_rc(rev_primer)])

    outer_argument = "SEQUENCE_ID=primer_pair\\nPRIMER_TASK=check_primers\\nSEQUENCE_PRIMER={}\\nSEQUENCE_PRIMER_REVCOMP={}\\nSEQUENCE_TEMPLATE={}\\nPRIMER_MIN_TM=59.5\\nPRIMER_MAX_TM=62.5\\nPRIMER_OPT_TM=61\\nPRIMER_OPT_SIZE=24\\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(fwd_primer,rev_primer,construct)
    outer_output = subprocess.check_output("primer3_core <(printf \"{}\")".format(outer_argument),shell=True,executable="/bin/bash").decode()
    outer_penalty = 100
    if not re.findall('PRIMER_PAIR_NUM_RETURNED=0',outer_output):
        lines = outer_output.split("\n")[1:-2]
        param = {}
        for line in lines:
            parts = line.split("=")
            param[parts[0]] = parts[1]
        outer_penalty = float(param['PRIMER_PAIR_0_PENALTY'])

    inner_argument = "SEQUENCE_ID=primer_pair\\nPRIMER_TASK=check_primers\\nSEQUENCE_PRIMER={}\\nSEQUENCE_PRIMER_REVCOMP={}\\nSEQUENCE_TEMPLATE={}\\nPRIMER_MIN_TM=59.5\\nPRIMER_MAX_TM=62.5\\nPRIMER_OPT_TM=61\\nPRIMER_OPT_SIZE=24\\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(fwd_primer,tag_primer,construct)
    inner_output = subprocess.check_output("primer3_core <(printf \"{}\")".format(inner_argument),shell=True,executable="/bin/bash").decode()
    inner_penalty = 100
    if not re.findall('PRIMER_PAIR_NUM_RETURNED=0',inner_output):
        lines = inner_output.split("\n")[1:-2]
        param = {}
        for line in lines:
            parts = line.split("=")
            param[parts[0]] = parts[1]
        inner_penalty = float(param['PRIMER_PAIR_0_PENALTY'])
    
    return((construct, fwd_primer, fwd_spacer, tag_spacer, tag_primer, rev_spacer, rev_primer, outer_penalty, inner_penalty))

def gen_primers(l,n,pool, path='.'):
    seqs = map(lambda x: gen_primer(l), range(n))
    seqs = filter(check_gc, seqs)
    seqs = filter(check_palindrome, seqs)
    seqs = filter(check_runs, seqs)
    seqs = filter(check_hairpins, seqs)
    seqs = list(dict.fromkeys(seqs))
    sys.stdout.write("{} of {} random {}mers passed initial checks.\n".format(len(seqs), n, l))
    aln = run_bwa(seqs, path)
    seqs = filter(lambda x: check_misprime(x, aln), seqs)
    seqs = pool.map(check_primer_primer3, seqs)
    seqs = sorted(list(seqs), key=lambda x:x[1])
    return(seqs)

def gen_spacers(l,n):
    seqs = map(lambda x: gen_seq(l), range(n))
    seqs = filter(check_gc, seqs)
    seqs = filter(check_palindrome, seqs)
    seqs = list(dict.fromkeys(seqs))
    sys.stdout.write("{} of {} random {}mers passed initial checks.\n".format(len(seqs), n, l))
    return(list(seqs))



if __name__ == "__main__":
    # Set up parallelisation
    pool = mp.Pool(32)

    # Generate a Fwd Primer
    sys.stdout.write("Generating forward primers..\n")
    fwd_primers = gen_primers(24,1000000,pool)
    sys.stdout.write("Generated {} forward primers.\n\n".format(len(fwd_primers)))
    sys.stdout.flush()

    # Generate a spacer
    sys.stdout.write("Generating forward spacers..\n")
    fwd_spacers = gen_spacers(24,1000)
    sys.stdout.write("Generated {} forward spacers.\n\n".format(len(fwd_spacers)))
    sys.stdout.flush()

    # Generate first part of tag
    sys.stdout.write("Generating tag spacers..\n")
    tag_spacers = gen_spacers(16,100000)
    sys.stdout.write("Generated {} tag spacers.\n\n".format(len(tag_spacers)))
    sys.stdout.flush()

    # Generate tag primers
    sys.stdout.write("Generating tag primers..\n")
    tag_primers = gen_primers(24,1000000,pool)
    sys.stdout.write("Generated {} tag primers.\n\n".format(len(tag_primers)))
    sys.stdout.flush()

    # Generate a spacer
    sys.stdout.write("Generating reverse spacers..\n")
    rev_spacers = gen_spacers(8,1000)
    sys.stdout.write("Generated {} reverse spacers.\n\n".format(len(rev_spacers)))
    sys.stdout.flush()

    # Generate a Rev Primer
    sys.stdout.write("Generating reverse primers..\n")
    rev_primers = gen_primers(24,1000000,pool)
    sys.stdout.write("Generated {} reverse primers.\n\n".format(len(rev_primers)))
    sys.stdout.flush()

    # Take the best 10 fwd and rev primers, pair and check
    sys.stdout.write("Determining best outer primer pair..\n")
    fwd_rev_pairs = itertools.product(fwd_primers[0:10],rev_primers[0:10])
    fwd_rev_pairs = pool.starmap(check_pair_primer3,[(x[0][0],x[1][0]) for x in fwd_rev_pairs])
    fwd_rev_pairs = sorted(list(fwd_rev_pairs), key=lambda x:x[2])

    best_fwd_primer = fwd_rev_pairs[0][0]
    best_rev_primer = fwd_rev_pairs[0][1]

    # Construct and check full sequences
    sys.stdout.write("Checking {} final constructs..\n".format(min([len(tag_spacers),len(tag_primers)])))
    final_constructs = pool.starmap(check_final_construct,[(best_fwd_primer,fwd_spacers[0],x[0],x[1][0],rev_spacers[0],best_rev_primer) for x in zip(tag_spacers,tag_primers)])
    final_constructs = sorted(list(final_constructs), key=lambda x:x[8])
    sys.stdout.write("Finished! Generated {} final constructs.\n\n".format(len(final_constructs)))

    # Output
    with open("results.txt", 'w') as fo:
        fo.write("Full Sequence\tFwd Primer\tFwd Spacer\tTag Spacer\tTag Primer\tRev Spacer\tRev Primer\tOuter Score\tTag Score\n")
        for fc in final_constructs:
            for x in fc:
                fo.write("{}\t".format(x))
            fo.write("\n")

    sys.stdout.write("Completed in {} seconds.\n".format(time.time()-start))
