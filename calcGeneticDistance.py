#! usr/bin/python


import dendropy
from dendropy.calculate import treemeasure
import os.path
import cPickle as pickle


path = '/home/jaggu/jagannath/operons'
fname = 'bactTree.mod.1.newick'
treeF = os.path.join(path,fname)

tns = dendropy.TaxonNamespace()
tree = dendropy.Tree.get_from_path(treeF,'newick',taxon_namespace=tns)
pdm = treemeasure.PatristicDistanceMatrix(tree)
print "Calculation done"

db = dict()

for t1 in tns:
	for t2 in tns:
		db[(t1,t2)]=pdm(t1,t2)

		
pklF = os.path.join(path,'all.dist.pkl')
pickle.dump(db, open(pklF,'wb'))
print "Pickled"
 
