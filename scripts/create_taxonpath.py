import argparse
import json

parser = argparse.ArgumentParser(description="Convert nodes.dmp to taxonpath.json")
parser.add_argument("-i", dest="input", help="Input nodes.dmp file name", required=True)
parser.add_argument("-o", dest="output", help="Output taxonpath.json file name", required=True)
args = parser.parse_args()

DISPLAY_RANKS = ['k', 'p', 'c', 'o', 'f', 'g', 's']
RANK_FULL_NAMES = {"k": "superkingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}

class TaxonEntry(object):
    def __init__(self, taxonId):
        self.taxonId = taxonId
        self.parent = ""
        
    def setRank(self, rank):
        self.rank = rank
        
    def setParent(self, parent):
        self.parent = parent
        
    def getTaxonId(self):
        return self.taxonId
        
    def getRank(self):
        return self.rank
        
    def getParent(self):
        return self.parent

    def getAllParents(self):
        ret = []
        if self.parent:
            ret.append(self.parent)
            ret.extend(self.parent.getAllParents())
        return ret
        
    def getTaxonPath(self):
        map = {}
        for x in self.getAllParents():
            map[x.getRank()] = x.getTaxonId()
        
        map[self.rank] = self.getTaxonId()
        
        ret = {'rank': self.getRank()}
        for r in DISPLAY_RANKS:
            ret[r] = map.get(RANK_FULL_NAMES[r], "")
            
        return ret

taxonEntryMap = {}
with open(args.input) as file:
    while line := file.readline():
        cols = line.split('|')
        taxonId = cols[0].strip()
        parentId = cols[1].strip()
        rank = cols[2].strip()

        if taxonId not in taxonEntryMap:
            taxonEntryMap[taxonId] = TaxonEntry(taxonId)

        if taxonId != parentId:
            if parentId not in taxonEntryMap:
                taxonEntryMap[parentId] = TaxonEntry(parentId)

            taxonEntryMap[taxonId].setParent(taxonEntryMap[parentId])

        taxonEntryMap[taxonId].setRank(rank)

data = {}
for tid in taxonEntryMap:
    rank = taxonEntryMap[tid].getRank()
    if rank in RANK_FULL_NAMES.values():
        data[tid] = taxonEntryMap[tid].getTaxonPath()

with open(args.output, "w") as j:
    json.dump(data, j, indent=3)
