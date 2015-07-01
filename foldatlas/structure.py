# Contains code for analysing RNA structure predictions.
# @author Matthew Norris <matthew.norris@jic.ac.uk>

from sklearn import decomposition

data = [
	[0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1],
	[0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1],
	[0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1],
	[0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1],
	[0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0],
	[0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1],
	[0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1],
	[1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0],
	[1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
	[1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1]
]

def do_pca(data):
	pca = decomposition.PCA(n_components=2)
	pca.fit(data)
	results = pca.transform(data)

	print(results)

do_pca(data);





