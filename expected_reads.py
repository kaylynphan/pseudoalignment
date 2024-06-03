from input import parse_reads, parse_transcriptome

def verify_expected_reads():
    read_data = parse_reads('reads.fasta')
    transcriptome_data = parse_transcriptome('chr11_transcriptome.fasta')

    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    sequence_map = {}

    for item in transcriptome_data:
        isoform_name = item['isoform_name']
        fwd_strand = item['sequence']
        rev_strand = str([complement[ch] for ch in fwd_strand[::-1]])
        
        sequence_map[isoform_name] = (fwd_strand, rev_strand)


    for item in read_data:
        read_id = item['read_id']
        isoform_name = item['isoform_name']
        mate1_start = int(item['mate1_start'])
        mate1_end = int(item['mate1_end'])
        mate2_start = int(item['mate2_start'])
        mate2_end = int(item['mate2_end'])
        sequence = item['sequence']

        print(f'Read ID: {read_id}')
        print("Expected Read:")

        (fwd_strand, rev_strand) = sequence_map[isoform_name]
        print("Mate 1:")
        print(fwd_strand[mate1_start: mate1_end])

        print("Mate 2:")
        print(rev_strand[mate2_start: mate2_end])

        print("Actual Read Sequence:")
        print(sequence)


if __name__ == '__main__':
    verify_expected_reads()
