from input import parse_reads, parse_transcriptome
import argparse

def verify_expected_reads(reads_file, transcriptome_file):
    read_data = parse_reads(reads_file)
    transcriptome_data = parse_transcriptome(transcriptome_file)

    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    sequence_map = {}

    for item in transcriptome_data[0:min(10, len(transcriptome_data))]:
        isoform_name = item['isoform_name']
        fwd_strand = item['sequence']
        sequence_map[isoform_name] = fwd_strand

    for item in read_data[0:min(10, len(read_data))]:
        read_id = item['read_id']
        isoform_name = item['isoform_name']
        mate1_start = int(item['mate1_start'])
        mate1_end = int(item['mate1_end'])
        mate2_start = int(item['mate2_start'])
        mate2_end = int(item['mate2_end'])
        sequence = item['sequence']

        print(f'Read ID: {read_id}')
        print(f'Matched to: {isoform_name}')
        print("Expected Read:")

        fwd_strand = sequence_map[isoform_name]
        print("Mate 1:")
        print(fwd_strand[mate1_start-1: mate1_end])

        print("Mate 2:")
        print(''.join([complement[ch] for ch in fwd_strand[mate2_end-1:mate2_start-3:-1]]))

        print("Actual Read Sequence:")
        print(''.join(sequence))
        
        print("----------------------------------------------")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Adding optional argument
    parser.add_argument("--reads", dest='read_file', type=str)
    parser.add_argument("--transcriptome", dest='transcriptome_file', type=str)

    args = parser.parse_args()
    verify_expected_reads(args.read_file, args.transcriptome_file)
