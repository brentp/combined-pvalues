from toolshed import reader
import sys
import re

header = open(sys.argv[1]).readline().rstrip()
print header + "\tnearest\tgene"
# calling bedtools
for r in reader("|closestBed -a %s -b data/arabidopsis_thaliana.gff -d -t all \
                 | groupBy -g 1,2,3,4 -c 7,14,13 -o collapse,min,collapse" \
                 % sys.argv[1], header=False):
    names = ",".join(set(re.findall("(?:ID|Name|Parent)=(\w+)", r[-1])))
    r[-1] = names

    # try to get more specific types.
    ftypes = r[5].split(",")
    if len(ftypes) > 2 and "gene" in ftypes:
        ftypes = [f for f in ftypes if f not in ("gene", "mRNA", "protein")]
        if len(ftypes): r[5] = ",".join(ftypes)

    dist = r[-2]
    if dist != '0':
        r[5] = dist
    print "\t".join(r[:6] + [r[-1]])

