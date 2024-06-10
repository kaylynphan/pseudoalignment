import re

# Define a function to parse the fasta-like file
def parse_reads(file_path):
    with open(file_path, 'r') as file:
        data = file.read().strip().split('>')
    
    parsed_data = []

    for entry in data:
        if not entry:
            continue
        lines = entry.split('\n')
        header = lines[0]
        sequence = ''.join(lines[1:])

        # print(f'header: {header}')

        # Use regular expressions to parse the header
        match = re.match(r'(\S+)/(\S+);mate1:(\d+)-(\d+);mate2:(\d+)-(\d+)', header)
        if match:
            read_id = match.group(1)
            isoform_name = match.group(2)
            mate1_start = match.group(3)
            mate1_end = match.group(4)
            mate2_start = match.group(5)
            mate2_end = match.group(6)

            parsed_data.append({
                'read_id': read_id,
                'isoform_name': isoform_name,
                'mate1_start': mate1_start,
                'mate1_end': mate1_end,
                'mate2_start': mate2_start,
                'mate2_end': mate2_end,
                'sequence': sequence
            })
    
    return parsed_data

def parse_transcriptome(file_path):
    with open(file_path, 'r') as file:
        parsed_data = []
        current_entry = {}

        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_entry:
                    parsed_data.append(current_entry)
                current_entry = {'isoform_name': line[1:], 'sequence': ''}
            else:
                current_entry['sequence'] += line

        # Add the last entry
        if current_entry:
            parsed_data.append(current_entry)

    return parsed_data
