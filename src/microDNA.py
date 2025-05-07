import pysam
import argparse
import find_start_end
import get_seq
import sw

MAX_SIZE_THRESHOLD = 1000            # max length of circular microDNA
MIN_SIZE_THRESHOLD = 100             # min length of circular microDNA
MIN_SCORE_THRESHOLD = 60
NO_MICROH_PENALTY = 0.2

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_file",
                        type=str,
                        default="data/SRR413984.sorted.NC_000001.10.bam",
                        help="BAM file to read")
    parser.add_argument("--fasta_file",
                        type=str,
                        default="data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna",
                        help="FASTA (.fna) file to read")

    return parser.parse_args()

def get_alignment_score(ref_genome, start_pos, end_pos, s_s, s_e):
    # Get the sequence from the reference genome
    seq = ref_genome.fetch("NC_000001.10", start_pos, end_pos)

    start_seq = seq[:10].upper()
    end_seq = seq[-10:].upper()

    H = sw.sw_fill_matrix(start_seq, end_seq, -2, -1, 1)
    align_A, align_B, score = sw.sw_traceback(H,start_seq, end_seq, -2, -1, 1)

    micro_h = align_A

    if micro_h == start_seq[:len(micro_h)]:

        c_begin = micro_h + s_e
        c_begin = c_begin.upper()
        c_end = s_s + micro_h
        c_end = c_end.upper()

        start_seq = seq[:len(c_begin)].upper()

        H = sw.sw_fill_matrix(start_seq, c_begin, -2, -1, 1)
        align_A, align_B, score = sw.sw_traceback(H,start_seq, c_begin, -2, -1, 1)

        score1 = score/len(start_seq)

        end_seq = seq[-1*len(c_end):].upper()

        H = sw.sw_fill_matrix(end_seq, c_end, -2, -1, 1)
        align_A, align_B, score = sw.sw_traceback(H, end_seq, c_end, -2, -1, 1)

        score2 = score/len(end_seq)

        return (score1 + score2)/2

    else:   # no microhomology present
        c_begin = s_e.upper()
        c_end = s_s.upper()

        start_seq = seq[:len(c_begin)].upper()

        H = sw.sw_fill_matrix(start_seq, c_begin, -2, -1, 1)
        align_A, align_B, score = sw.sw_traceback(H,start_seq, c_begin, -2, -1, 1)

        score1 = score/len(start_seq)         

        end_seq = seq[-1*len(c_end):].upper()
        
        H = sw.sw_fill_matrix(end_seq, c_end, -2, -1, 1)
        align_A, align_B, score = sw.sw_traceback(H, end_seq, c_end, -2, -1, 1)

        score2 = score/len(end_seq)      

        return (score1 + score2)/2 - NO_MICROH_PENALTY     # take off points because no microhomology


def main():

    args = get_args()

    junctions = find_start_end.find_junctions(args.bam_file)

    ref_genome = pysam.FastaFile(args.fasta_file)

    circles = []
    scoredCircles = []

    
    # find max/min number of junction tags for normalization
    maxNumTags = junctions[1][3]
    minNumTags = junctions[1][3]

    for i, curr_junc in enumerate(junctions):
        if curr_junc[3] > maxNumTags:
            maxNumTags = curr_junc[3]
        if curr_junc[3] < minNumTags:
            minNumTags = curr_junc[3]

        if curr_junc[0] == '+':
            k = 0
            if i + 1 >= len(junctions):
                break
            
            next_junc = junctions[i+1]

            start_pos = curr_junc[1]

            while next_junc[1]-start_pos <= MAX_SIZE_THRESHOLD and (i+k) < (len(junctions)):
                if next_junc[0] == '-':
                    end_pos = next_junc[1]
                    score = get_alignment_score(ref_genome, start_pos, end_pos, curr_junc[2], next_junc[2])
                    if end_pos-start_pos >= MIN_SIZE_THRESHOLD:
                        circles.append([start_pos, end_pos, score, curr_junc[3]+next_junc[3]])
                        
                        
                k += 1
                if i+k < len(junctions):
                    next_junc = junctions[i+1+k]


    # final scoring of circles
    for circle in circles:
        # find normalized junction score (take avg of start/end junction tags)
        junction_score = ( (circle[3]/2) - minNumTags) / (maxNumTags - minNumTags)

        # add criteria together with weights
        norm_score = 80*circle[2] + 20*junction_score

        if norm_score > MIN_SCORE_THRESHOLD:
            scoredCircles.append([circle[0], circle[1], norm_score])



    ## PRINT SCORED CIRCLES ABOVE THRESHOLD##
    print('{:25s} {:12s} {:12s}'.format("POSITION", "LENGTH", "SCORE"))
    for circle in scoredCircles:
        pos = str(circle[0]) + "-" + str(circle[1])
        circ_length = circle[1]-circle[0]
        cl_str = str(circ_length) + " bp"
        score = circle[2]
        print('{:25s} {:12s} {:5.2f}'.format(pos, cl_str, score))

    print('\nNumber of circles:', len(scoredCircles))

    # for circle in scoredCircles:
    #     print(circle, "Length (bp):", (circle[1]-circle[0]))
    # print('Number of circles:', len(scoredCircles))

if __name__ == "__main__":
    main()
