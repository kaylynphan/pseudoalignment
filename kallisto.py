from input import parse_reads, parse_transcriptome

class Node:
  def __init__(self, kmer):
    self.kmer = kmer
    self.equivalence_class = set()
    self.next = set()

  def add_isoform(self, isoform_name):
    self.equivalence_class.add(isoform_name)
    
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
    self.contig_mapping[(isoform_name, 0)] = prev
    self.contig_lengths[isoform_name] = 1
    
    for i in range(1, len(sequence) - self.k - 1):
      kmer = sequence[i:i+self.k]
      node = self.get_node(kmer)

      node.add_isoform(isoform_name)
      prev.next.add(node)
      prev = node

  def classify_read(self, read_sequence):
    if len(read_sequence) < self.k:
      ValueError(f"Read length is less than k")
    
    read_equivalence_class = set()
    # find head
    for i in range(0, len(read_sequence) - self.k - 1):
      kmer = read_sequence[i:i+self.k]
      if kmer not in self.kmer_to_node_map:
        # kmer does not exist in any transcript
        return read_equivalence_class # empty set
      else:
        start = self.kmer_to_node_map[kmer]
        # start a dfs from start node. When a cycle is detected, stop the branch
        stack = [start]
        seen_kmers = set(start.kmer)
        while stack:
          curr = stack.pop()

          # perform kallisto 'skipping'
          children = set()
          for eq in curr.equivalence_class:
            contig_length = self.contig_lengths[eq]
            last_node_of_contig = self.contig_mapping[(eq, contig_length - 1)]
            node_following_last_node_of_contig = last_node_of_contig.next
            children.add(node_following_last_node_of_contig)
          
          # add to stack, include cycle prevention
          for c in children:
            if c.kmer not in seen_kmers:
              stack.append(c)
              seen_kmers.add(c.kmer)

          read_equivalence_class = read_equivalence_class.intersection(curr.equivalence_class)
          
        return read_equivalence_class


def run_kallisto(k):
  read_data = parse_reads('reads.fasta')
  transcriptome_data = parse_transcriptome('chr11_transcriptome.fasta')

  tbg = DeBruijnGraph(k)

  for transcript in transcriptome_data:
    sequence = transcript['sequence']
    isoform_name = transcript['isoform_name']

    tbg.add_transcript(sequence, isoform_name)

  for item in read_data:
    read_sequence = item['sequence']
    isoforms = tbg.classify_read(read_sequence)
