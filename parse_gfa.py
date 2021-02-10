def parse_gfa(gfa_file):
	nodes = {}
	edges = []
	with open(gfa_file, "r") as f:
		for line in f:
			line = line.strip()
			parsed_line = line.split("\t")
			if parsed_line[0] == "S":
				id = int(parsed_line[1])
				seq = parsed_line[2]
				nodes[seq] = id
			elif parsed_line[0] == "L":
				if parsed_line[2] == "+":
					source_dir = "F"
				else:
					source_dir = "R"
				if parsed_line[4] == "+":
					sink_dir = "F"
				else:
					sink_dir = "R"
				edges.append([parsed_line[1], parsed_line[3], source_dir + sink_dir])
	return(nodes, edges)

def parse_unitig_coverage(coverage_file):
	unitig_coverage = {}
	with open(coverage_file, "r") as f:
		for line in f:
			line = line.strip()
			parsed_line = line.split("\t")
			unitig_coverage[parsed_line[0]] = parsed_line[1]
	return(unitig_coverage)


def generate_output(gfa_file, coverage_file, outpref):
	nodes, edges = parse_gfa(gfa_file)
	unitig_coverage = parse_unitig_coverage(coverage_file)

	nodes_output = outpref + ".nodes"
	edges_output = outpref + ".edges"
	cov_output = outpref + ".coverage"

	with open(outpref + ".nodes", "w") as n:
		for node, id in nodes.items():
			n.write(str(id) + "\t" + str(node) + "\n")

	with open(outpref + ".edges", "w") as e:
		for edge in edges:
			e.write(str(edge[0]) + "\t" + str(edge[1]) + "\t" + str(edge[2]) + "\n")

	with open(outpref + ".coverage", "w") as c:
		for seq, cov in unitig_coverage.items():
			c.write(str(nodes[seq]) + "\t" + str(cov) + "\n")


if __name__ == "__main__":
	import sys
	args = sys.argv
	gfa_file = args[1]
	coverage_file = args[2]
	outpref = args[3]
	generate_output(gfa_file, coverage_file, outpref)
