from input import parse_reads, parse_transcriptome
import argparse

class Node:
  def __init__(self, kmer):
    self.kmer = kmer
    self.equivalence_class = set()
    self.next = set()

  def __repr__(self):
    return f'kmer: {self.kmer} | equivalence class: {self.equivalence_class} | next: {[n.kmer for n in self.next]}'

  def add_isoform(self, isoform_name):
    self.equivalence_class.add(isoform_name)
    

# skipping not implemented yet
class DeBruijnGraph:
  def __init__(self, k):
    self.k = k

    # Mapping of (isoform_name, kallisto index) to node
    # e.g. ('ENST000003679', 0) --> Node. This node is the first occurring in the De Bruijn graph whose equivalence class contains ENST000003679

    self.contig_mapping = {}
    self.contig_lengths = {}

    # Mapping of kmer to node e.g. AGC --> Node
    self.kmer_to_node_map = {}

  def get_node(self, kmer):
    if kmer in self.kmer_to_node_map:
      return self.kmer_to_node_map[kmer]
    else:
      new_node = Node(kmer)
      self.kmer_to_node_map[kmer] = new_node
      return new_node

  def add_transcript(self, sequence, isoform_name):
    if len(sequence) < self.k:
      ValueError(f"Transcript length is less than k")
    
    first_kmer = sequence[0:self.k]
    # find location of first kmer
    prev = self.get_node(first_kmer)
    prev.add_isoform(isoform_name)
    # self.contig_mapping[(isoform_name, 0)] = prev
    # self.contig_lengths[isoform_name] = 1
    
    for i in range(1, len(sequence) - self.k + 1):
      kmer = sequence[i:i+self.k]
      # print(f"kmer: {kmer}")
      node = self.get_node(kmer)
      node.add_isoform(isoform_name)
      # self.contig_mapping[(isoform_name, i)] = node
      # self.contig_lengths[isoform_name] += 1
      prev.next.add(node)
      prev = node


  def print_nodes(self):
    print("TBG Nodes:")
    for node in self.kmer_to_node_map.values():
      print(node)

  def classify_read(self, read_sequence):
    print(f'classifying read {read_sequence}')
    if len(read_sequence) < self.k:
      ValueError(f"Read length is less than k")
    
    read_equivalence_class = set()

    # find head
    first_kmer = read_sequence[0:self.k]
    # print(f"first kmer: {first_kmer}")
    if first_kmer not in self.kmer_to_node_map:
      print(f'returning empty set because kmer {first_kmer} is not in self.kmer_to_node_map')
      return set() # empty set
    else:
      # initialize equivalence class set
      first_node = self.kmer_to_node_map[first_kmer]
      read_equivalence_class = first_node.equivalence_class

    for i in range(1, len(read_sequence) - self.k + 1):
      kmer = read_sequence[i:i+self.k]
      # print(f"kmer: {kmer}")
      if kmer not in self.kmer_to_node_map:
        # kmer does not exist in any transcript
        print(f'returning empty set because kmer {kmer} not found')
        return set() # empty set
      else:
        node = self.kmer_to_node_map[kmer]
        # start a dfs from start node. When a cycle is detected, stop the branch
        # print(f'intersecting with node {node}')
        read_equivalence_class = read_equivalence_class.intersection(node.equivalence_class)
    
        # # perform kallisto 'skipping'
        # children = set()
        # for eq in curr.equivalence_class:
        #   contig_length = self.contig_lengths[eq]
        #   last_node_of_contig = self.contig_mapping[(eq, contig_length - 1)]
        #   node_following_last_node_of_contig = last_node_of_contig.next
        #   children.add(node_following_last_node_of_contig)
              
    return read_equivalence_class


def run_kallisto(k, reads_file, transcriptome_file, test):
  read_data = parse_reads(reads_file)
  # print(read_data)
  transcriptome_data = parse_transcriptome(transcriptome_file)

  tbg = DeBruijnGraph(k)

  if test: # truncate for ease of reading
    transcriptome_data = transcriptome_data[0:min(10, len(transcriptome_data))]
    read_data = read_data[0:min(10, len(read_data))]

  for transcript in transcriptome_data:
    sequence = transcript['sequence']
    isoform_name = transcript['isoform_name']

    print(f'Adding transcript for isoform {isoform_name}')

    tbg.add_transcript(sequence, isoform_name)

  tbg.print_nodes()

  for item in read_data:
    read_id = item['read_id']
    read_sequence = item['sequence']
    isoforms = tbg.classify_read(read_sequence)
    print(f"{read_id} matched to isoforms {isoforms}")

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  # Adding optional argument
  parser.add_argument("--k", dest='k', type=int)
  parser.add_argument("--reads", dest='read_file', type=str)
  parser.add_argument("--transcriptome", dest='transcriptome_file', type=str)
  parser.add_argument("--test", action='store_true', default=False,
        help="Use sabre to get SWAP upper bound")

  args = parser.parse_args()

  run_kallisto(args.k, args.read_file, args.transcriptome_file, args.test)
