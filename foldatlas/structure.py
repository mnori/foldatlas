# Contains code for analysing RNA structure predictions.
# @author Matthew Norris <matthew.norris@jic.ac.uk>

from sklearn import decomposition

data = [
	[0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1],
	[0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
	[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
	[0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
	[1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
	[1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1],
	[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1],
	[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
	[1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1],
	[1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1]
]

def do_pca(data):

	# hmm seems to work. BUT how do you know which entries are which?
	# Simple - results always listed in the order that they were added.

	pca = decomposition.PCA(n_components=2)
	pca.fit(data)
	results = pca.transform(data)

	print(results)

do_pca(data);

# from sklearn import datasets
# iris = datasets.load_iris()
# do_pca(iris.data);
# print(iris.target)





