import subprocess
import re


rc_dict = {'A': 'T','C':'G','G':'C','T':'A'}

def fast_rc(seq):
    if not seq:
        return ''
    seq_rc = ''.join([rc_dict[x] for x in seq])
    return(seq_rc[::-1])



def check_primer_primer3(seq):
    argument = "SEQUENCE_ID=primer_pair\\n" \
               "PRIMER_TASK=check_primers\\n" \
               "SEQUENCE_PRIMER={}\\n" \
               "PRIMER_MIN_TM=59.5\\n" \
               "PRIMER_MAX_TM=80.5\\n" \
               "PRIMER_OPT_SIZE=24\\n" \
               "PRIMER_OPT_TM=71\\n" \
               "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" \
               "/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(seq)

    output = subprocess.check_output("primer3_core <(printf \"{}\")".format(argument),shell=True,executable="/bin/bash").decode()
    primer_penalty = 100
    if not re.findall('PRIMER_LEFT_NUM_RETURNED=0',output):
        lines = output.split("\n")[1:-2]
        param = {}
        for line in lines:
            parts = line.split("=")
            param[parts[0]] = parts[1]
        primer_penalty = float(param['PRIMER_LEFT_0_PENALTY'])
    return ((seq,primer_penalty))


def check_pair_primer3(fwd,rev):
    argument = "SEQUENCE_ID=primer_pair\\n" \
               "PRIMER_TASK=check_primers\\n" \
               "SEQUENCE_PRIMER={}\\n" \
               "SEQUENCE_PRIMER_REVCOMP={}\\n" \
               "PRIMER_MIN_TM=59.5\\n" \
               "PRIMER_MAX_TM=62.5\\n" \
               "PRIMER_OPT_TM=61\\n" \
               "PRIMER_OPT_SIZE=24\\n" \
               "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" \
               "/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(fwd,rev)

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


def check_final_construct(fwd_primer, fwd_spacer, tag_spacer, tag_primer, rev_spacer, rev_primer,
                          ill_fwd='', ill_rvr=''):

    construct = ''.join([ill_fwd, fwd_primer, fwd_spacer, tag_spacer, fast_rc(tag_primer), rev_spacer, fast_rc(rev_primer), fast_rc(ill_rvr)])

    outer_argument = "SEQUENCE_ID=primer_pair\\n" \
                     "PRIMER_TASK=check_primers\\n" \
                     "SEQUENCE_PRIMER={}\\n" \
                     "SEQUENCE_PRIMER_REVCOMP={}\\n" \
                     "SEQUENCE_TEMPLATE={}\\n" \
                     "PRIMER_MIN_TM=59.5\\n" \
                     "PRIMER_MAX_TM=62.5\\n" \
                     "PRIMER_OPT_TM=61\\n" \
                     "PRIMER_OPT_SIZE=24\\n" \
                     "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" \
                     "/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(
        fwd_primer, rev_primer, construct)
    outer_output = subprocess.check_output("primer3_core <(printf \"{}\")".format(outer_argument), shell=True,
                                           executable="/bin/bash").decode()
    outer_penalty = 100
    if not re.findall('PRIMER_PAIR_NUM_RETURNED=0', outer_output):
        lines = outer_output.split("\n")[1:-2]
        param = {}
        for line in lines:
            parts = line.split("=")
            param[parts[0]] = parts[1]
        outer_penalty = float(param['PRIMER_PAIR_0_PENALTY'])

    inner_argument = "SEQUENCE_ID=primer_pair\\n" \
                     "PRIMER_TASK=check_primers\\n" \
                     "SEQUENCE_PRIMER={}\\n" \
                     "SEQUENCE_PRIMER_REVCOMP={}\\n" \
                     "SEQUENCE_TEMPLATE={}\\n" \
                     "PRIMER_MIN_TM=59.5\\n" \
                     "PRIMER_MAX_TM=62.5\\n" \
                     "PRIMER_OPT_TM=61\\n" \
                     "PRIMER_OPT_SIZE=24\\n" \
                     "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" \
                     "/nfs/modules/modules/software/Primer3/2.4.0-foss-2018b/primer3-2.4.0/src/primer3_config/\\n=".format(
        fwd_primer, tag_primer, construct)
    inner_output = subprocess.check_output("primer3_core <(printf \"{}\")".format(inner_argument), shell=True,
                                           executable="/bin/bash").decode()
    inner_penalty = 100
    if not re.findall('PRIMER_PAIR_NUM_RETURNED=0', inner_output):
        lines = inner_output.split("\n")[1:-2]
        param = {}
        for line in lines:
            parts = line.split("=")
            param[parts[0]] = parts[1]
        inner_penalty = float(param['PRIMER_PAIR_0_PENALTY'])

    return (
    (construct, fwd_primer, fwd_spacer, tag_spacer, tag_primer, rev_spacer, rev_primer, outer_penalty, inner_penalty))

if __name__ == "__main__":
    fwd ='TATGAGGAGAGTAGGAGGCAATGG'
    fwd_sapcer = 'GGAGGTTCACAATGTGGGAGGTCA'
    tag_spacer = 'AGACCAGTACATGACG'
    tag_primer = ['ATACAGGAGTGGCAGAGAGATACC', 'TCTCTTCTGGGTATGACGGTATCC', 'ACGGCATCCCATTATTCCTTTACC']
    rev_spacer = 'CTATTGCC'
    rev = 'ATCCCTACCCTTATCCCACTTACG'
    ill_fwd = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
    ill_rev = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    result = check_primer_primer3(ill_rev)
    print(result)
    # for t in tag_primer:
    #     result = check_final_construct(fwd, fwd, tag_spacer, t, rev_spacer, rev, ill_fwd = 'TCGTCGGCAGCGTCAGATGTG',
    #                                    ill_rvr='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG')
    #     print(result[0], result[-2:])

