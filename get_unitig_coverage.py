def get_unitig_coverage(graph_infile, coverage_infile, kmer, coverage_outfile):
    from Bio import SeqIO
    from Bio.Seq import Seq
    unitig_cov_dict = {}
    kmer_coverage_dict = {}
    with open(graph_infile, "r") as g:
        for line in g:
            if line[0] == "S":
                line = line.strip("\n")
                items = line.split("\t")
                id = items[1]
                seq = items[2]
                unitig_cov_dict[seq] = 0

    records = list(SeqIO.parse(coverage_infile, "fasta"))
    for rec in records:
        coverage = int(rec.id)
        kmer_seq = str(rec.seq)
        kmer_coverage_dict[kmer_seq] = coverage

    #iterate over untig kmers, calculate average coverage
    for unitig in unitig_cov_dict.keys():
        num_kmers = len(unitig) - kmer + 1
        total_coverage = 0

        kmer_cov = 0
        for i in range(0, num_kmers):
            kmer_seq = unitig[i:i + kmer]
            if kmer_seq in kmer_coverage_dict:
                kmer_cov = kmer_coverage_dict[kmer_seq]
            else:
                kmer_seq = str(Seq(kmer_seq).reverse_complement())
                kmer_cov = kmer_coverage_dict[kmer_seq]

            total_coverage += kmer_cov

        average_coverage = total_coverage / num_kmers
        unitig_cov_dict[unitig] = average_coverage

    with open(coverage_outfile, "w") as o:
        for unitig, coverage in unitig_cov_dict.items():
            o.write(unitig + "\t" + str(coverage) + "\n")
