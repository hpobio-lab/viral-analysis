import sys
from collections import defaultdict


if __name__ == "__main__":

    header = "#ref\tname\tgene\tstart\tend\tstrand\tphase\tfeature_type"
    print(header)
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            if not line.startswith("#"):
                tokens = line.strip().split("\t")
                ref = tokens[0]
                feature_type = tokens[2]
                start = tokens[3]
                end = tokens[4]
                strand = tokens[6]
                phase = tokens[7]
                name = ""
                gene = ""
                attr = tokens[8]
                attr_d = defaultdict(str)
                for i in attr.split(";"):
                    i_splits = i.split("=")
                    if len(i_splits) > 1:
                        attr_d[i_splits[0]] = i_splits[1]
                    else:
                        attr_d[i_splits[0]] = i_splits[0]
                if "Name" in attr_d:
                    name = attr_d["Name"]
                if "gene" in attr_d:
                    gene = attr_d["gene"]
                print("\t".join([ref, name, gene, start, end, strand, phase, feature_type]))
