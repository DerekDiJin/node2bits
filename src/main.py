import sys
import datetime
from pathlib import Path
import numpy as np, scipy as sp, networkx as nx
import math, time, os, sys, random
from collections import deque
import pickle
import itertools
from util import *

import matplotlib
import matplotlib.pyplot as plt

import scipy.sparse as sps
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import svds, eigs
import sparsesvd
import scipy.spatial.distance as distance

# import sklearn
# from sklearn.model_selection import KFold
# from sklearn.linear_model import LogisticRegression
# from sklearn.naive_bayes import MultinomialNB
# from sklearn.model_selection import train_test_split
# from sklearn.model_selection import validation_curve
# from sklearn.multiclass import OneVsRestClassifier
# from sklearn.svm import SVC
# from sklearn import svm
# from sklearn import preprocessing
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF, DictionaryLearning
from sklearn.manifold import TSNE
from sklearn.preprocessing import normalize

from sklearn.feature_extraction import FeatureHasher
import collections
from collections import defaultdict


def aggregate_walks(walks, walks_num, nodes_to_explore, dist_scope):

	print '[Aggregating walks]'

	result_out = [ {i:[] for i in range(dist_scope)} for _ in nodes_to_explore ]
	# result_in = [ {i:[] for i in range(dist_scope)} for _ in nodes_to_explore ]
	'''
		Return (2 dicts): 
			node - {1:[nodes], 2:[nodes], 3:[nodes]}
		Parameters: 
			dist_scope: dist scopres to explore
	'''


	N = len( nodes_to_explore )

	counter = 0
	for walk in walks:

		if counter % 50000 == 0:
			print '[Current walk] ' + str(counter)
		counter += 1

		if len(walk) <= dist_scope:

			head = walk[0]

			for dist in range(1, len(walk)):
				result_out[head][dist-1].append( walk[dist] )


		else:
			for idx in range(len(walk) - dist_scope):
				
				walk_seq = walk[idx: idx+dist_scope+1]
				head = walk_seq[0]
				tail = walk_seq[-1]

				for dist in range(1, dist_scope+1):

					# print 'walk_seq:'
					# print walk_seq

					# print 'head: ' + str(head)
					# print 'tail: ' + str(tail)

					result_out[head][dist-1].append( walk_seq[dist] )
					# result_in[tail][dist].append( walk_seq[dist_scope-dist-1] )

	return result_out#, result_in

def get_sketch_b(K, P, idx):

	file_name = 'sketch_' + str(idx) + '.tsv'
	result = np.genfromtxt(file_name, delimiter='\t', dtype=int)

	return result


def get_sketch(K, P, idx):

	# result = sps.random(N, P, density=0.5, format='coo')
	# result = np.zeros((K, P))

	# for k in range(K):
	# 	tmp = np.random.choice([-1, 1], size=(P))
	# 	tmp_map = collections.Counter(tmp)

	# 	while (tmp_map[1] - tmp_map[-1]) > P * 0.05:
	# 		tmp = np.random.choice([-1, 1], size=(P))
	# 		tmp_map = collections.Counter(tmp)

	# 	result[k,:] = tmp

	file_name = 'sketch_' + str(idx) + '.tsv'
	tmp = np.random.choice([-1, 1], size=(K, P))
	np.savetxt(file_name, tmp, delimiter='\t', fmt='%i')
	return tmp

def bit_to_int(binary_list):
	# print binary_list
	return int(''.join([str(ele) for ele in binary_list]), 2)

def list_to_ints(l):
	'''
	Converts a list of values to str of ints
	'''
	# return ''.join(str(int(ele)) for ele in l)
	return [int(ele) for ele in l]

def combine_list(dict_in):
	result = []

	for ele in dict_in:
		result += dict_in[ele]
	return result


def hash_func(PSs, Ks, Ts):
	print '[Hasing cosine sketch begins]'

	b = 4
	tables = [defaultdict(list) for _ in range(len(PSs))]
	f_rep = None

	for idx in range(len(PSs)):

		PS = PSs[idx]
		K = Ks[idx]
		T = Ts[idx]
		B = K/b

		table = tables[idx]
		indices_back = [ b*i for i in range(B) ]
		indices = random.sample(indices_back, B/2)	# OR-construction: randomly select B-T bands
		print 'indices:'
		print indices

		sps.save_npz('./PS.npz', sps.coo_matrix(PS))

		N, P = PS.shape
		# projection = np.zeros((K, N), dtype=int)
		# result = np.zeros((N, K), dtype=int)

		sketches = get_sketch(K, P, idx)
		print ':('
		print sketches.shape
		projection = np.matmul(PS, sketches.T).astype(int)
		result = (np.sign(projection) + 1) / 2

		print result

		if idx == 0:
			f_rep = result
		else:
			f_rep = np.concatenate((f_rep, result), axis=1)

		##################################################################

		for node in range(N):
			row = result[node, :]
			signature = tuple(  [ bit_to_int(row[index:index + b]) for index in indices]  ) # sign format: tuple
			table[signature].append(node)




	
	return tables, f_rep
	

	# l = result.shape[1]
	# b = l/16

	# for i in range(N):
	# 	
	# 	
	# 	# print 'signature:'
	# 	# print signature
	# 	



def emb_svd(feature_matrix, k):
	temp = scipy.sparse.csc_matrix(feature_matrix)
	U,s,V = sparsesvd.sparsesvd(temp, k)
	# S = sps.diag(s,dtype=float)
	#####################

	# U,s,V = svds(feature_matrix, k)
	S = np.diag(s)
	
	# print '----'
	# print U
	# print S
	emb = np.dot(U.T, (S ** 0.5))
	# g_sum = np.dot((S**0.5), V)

	return emb


def get_init_features(graph, base_features, nodes_to_explore):

	init_feature_matrix = np.zeros((len(nodes_to_explore), len(base_features)))
	adj = graph.adj_matrix

	if "outdegree" in base_features:
		init_feature_matrix[:,0] = (adj.sum(axis=0).transpose() +  adj.sum(axis=1)).ravel()

	if "indegree" in base_features:
		init_feature_matrix[:,1] = adj.sum(axis=0).transpose().ravel()

	if "degree" in base_features:
		init_feature_matrix[:,2] = adj.sum(axis=1).ravel()

	return init_feature_matrix

def get_feature_n_buckets(feature_matrix, num_buckets, bucket_max_value):

	result_sum = 0
	result_ind = []
	result_cum = []
	N, cur_P = feature_matrix.shape

	if num_buckets is not None:
		for i in range(cur_P):
			result_cum.append(result_sum)
			n_buckets = max(bucket_max_value, int(math.log(max(max(feature_matrix[:,i]), 1), num_buckets) + 1))
			# print max(feature_matrix[:,i])
			result_sum += n_buckets
			result_ind.append(n_buckets)
			
			
	else:
		for i in range(cur_P):
			result_cum.append(result_sum)
			n_buckets = max(bucket_max_value, int( max(feature_matrix[:,i]) ) + 1)
			result_sum += n_buckets
			result_ind.append(n_buckets)

	return result_sum, result_ind, result_cum



def feature_binning(graph, init_feature_matrix, nodes_to_explore, S, i):

	feature_wid_sum, feature_wid_ind, feature_wid_cum = get_feature_n_buckets(init_feature_matrix, graph.num_buckets, graph.bucket_max_value)
	feature_matrix_seq = np.zeros([graph.num_nodes, feature_wid_sum * len(graph.cat_dict.keys())])

	N, P = init_feature_matrix.shape

	id_cat_dict = graph.id_cat_dict

	

	for node in nodes_to_explore:
		if node % 50000 == 0:
			print "[Generate combined feature vetor] node: " + str(node)

		# combined_feature_sequence = get_combined_feature_sequence(graph, node, 
			# input_dense_matrix = init_feature_matrix, feature_wid_ind = feature_wid_ind, S_list=S[node][i])#------< resume here.
		
		S_list = S[node][i]
		cur_neighbors = set(S_list)
		cur_neighbor_dict = dict([x,S_list.count(x)/float(len(S_list))] for x in set(S_list))
		
		for neighbor in cur_neighbors:

			cat_idx = id_cat_dict[neighbor]


			for p in range(P):

				feature_idx = feature_wid_cum[p]
				node_feature = init_feature_matrix[neighbor, p]

				if (graph.num_buckets is not None) and (node_feature != 0):
					bucket_index = max( int(math.log(node_feature, graph.num_buckets)), 0 )
				else:
					bucket_index = int(node_feature)

				bucket_index = min( bucket_index, (feature_wid_ind[p]-1) )


				global_idx = cat_idx * feature_wid_sum + feature_idx + bucket_index
				feature_matrix_seq[node, global_idx] += cur_neighbor_dict[neighbor]


	# np.set_printoptions(threshold=np.nan)

	# print '========='
	# print feature_matrix_seq

	return feature_matrix_seq




def construct_prox_structure(graph, nodes_to_explore, base_features, S_out, dist_scope):

	init_feature_matrix = get_init_features(graph, base_features, nodes_to_explore)

	feature_matrices = []

	for i in range(dist_scope):

		feature_matrices.append( feature_binning(graph, init_feature_matrix, nodes_to_explore, S_out, i) )
	
	
	# init_feature_matrix = init_feature_matrix_seq

	# U = emb_svd(init_feature_matrix, feature_rank)

	# kmeans = KMeans(n_clusters=num_roles, random_state=0).fit(U)
	# role_mapping = kmeans.labels_
	# result = dict(zip(nodes_to_explore, role_mapping))

	return feature_matrices

def parse_weighted_temporal(input_file_path, delimiter):

	check_eq = True
	num_nodes = 0
	num_edges = 0
	adj_matrix_global = None
	edge_time_dict = None
	time_edge_dict = None
	start_time = 0
	end_time = 0
	# node_roadmap = None


	raw = np.genfromtxt(input_file_path, dtype=int, delimiter=delimiter)
	print raw
	ROW, COL = raw.shape
	num_edges = ROW

	if COL == 3:
		print '[input_file does not contain timestamps. Processing as static graphs]'

		srcs = raw[:,0]
		dsts = raw[:,1]
		weis = raw[:,2]

		max_id = int(max(max(srcs), max(dsts)))
		num_nodes = max_id + 1
		print '[max_node_id] ' + str(max_id)
		print '[num_nodes] ' + str(num_nodes)

		if max(srcs) != max(dsts):
			srcs = np.append(srcs, max(max(srcs), max(dsts)))
			dsts = np.append(dsts, max(max(srcs), max(dsts)))
			weis = np.append(weis, 0)
			check_eq = False

		# srcs_tmp = srcs

		# srcs = np.append(srcs, dsts)
		# dsts = np.append(dsts, srcs_tmp)
		# weis = np.append(weis, weis)

		adj_matrix_global = sps.lil_matrix( sps.csc_matrix((weis, (srcs, dsts))))

	elif COL == 4:
		print '[input_file contains timestamps. Processing as dynamic graphs]'

		edge_time_dict = defaultdict(list)
		time_edge_dict = defaultdict(list)

		srcs = raw[:,0]
		dsts = raw[:,1]
		weis = raw[:,2]
		times = raw[:,3]

		start_time = min(times)
		end_time = max(times)

		max_id = int(max(max(srcs), max(dsts)))
		num_nodes = max_id + 1
		print '[max_node_id] ' + str(max_id)
		print '[num_nodes] ' + str(num_nodes)

		if max(srcs) != max(dsts):
			srcs = np.append(srcs, max(max(srcs), max(dsts)))
			dsts = np.append(dsts, max(max(srcs), max(dsts)))
			weis = np.append(weis, 0)
			check_eq = False

		adj_matrix_global = sps.lil_matrix( sps.csc_matrix((weis, (srcs, dsts))))

		fIn = open(input_file_path, 'r')
		lines = fIn.readlines()
		for line in lines:
			parts = line.strip('\r\n').split(delimiter)
			src = int(parts[0])
			dst = int(parts[1])
			wei = float(parts[2])
			timestamp = int(parts[3])

			edge = (src, dst, wei)
			edge_time_dict[edge].append(timestamp)
			time_edge_dict[timestamp].append(edge)


			# edge_r = (dst, src, wei)
			# edge_time_dict[edge_r].append(timestamp)
			# time_edge_dict[timestamp].append(edge_r)
			# tup = (dst, wei, timestamp)
			# node_roadmap[src].append(tup)

		fIn.close()


	else:
		sys.exit('[input_file format error. Please make sure the input file with the format <src, dst, wei> or <src, dst, wei, timestamps>')

	return check_eq, num_nodes, num_edges, adj_matrix_global, edge_time_dict, time_edge_dict, start_time, end_time


def svd_emb(PSs, Ks):


	for idx in range(len(Ks)):
		Ki = Ks[idx]
		PSi = PSs[idx]
		emb_temp = feature_2_embedding(PSi, Ki)

		if idx == 0:
			emb = emb_temp
		else:
			emb = np.concatenate((emb, emb_temp), axis=1)

	return emb


def feature_2_embedding(feature_matrix = None, k = 17):

	####### TODO #######
	temp = sps.csc_matrix(feature_matrix)
	U,s,V = sparsesvd.sparsesvd(temp, k)
	# S = sps.diag(s,dtype=float)
	#####################

	# U,s,V = svds(feature_matrix, k)
	S = np.diag(s)
	
	# print '----'
	# print U
	# print S
	emb = np.dot(U.T, (S ** 0.5))
	# g_sum = np.dot((S**0.5), V)

	return emb


def construct_cat(input_gt_path, delimiter):
	
	####################################################
	# Input: per line, 1) cat-id_init, id_end or 2) cat-id
	#
	# Return: 1) dict: cat-ids and 2) id-cat
	####################################################

	result = defaultdict(set)
	id_cat_dict = dict()

	fIn = open(input_gt_path, 'r')
	lines = fIn.readlines()
	for line in lines:

		parts = line.strip('\r\n').split(delimiter)
		if len(parts) == 3:
			cat = parts[0]
			node_id_start = parts[1]
			node_id_end = parts[2]

			for i in range( int(node_id_start), int(node_id_end)+1 ):
				result[ int(cat) ].add( i )
				id_cat_dict[i] = int(cat)

		elif len(parts) == 2:
			cat = parts[0]
			node_id = parts[1]

			result[int(cat)].add( int(node_id) )
			id_cat_dict[int(node_id)] = int(cat)



		else:
			sys.exit('Cat file format not supported')

	fIn.close()
	return result, id_cat_dict


if __name__ == '__main__':
	
	if len(sys.argv) != 1:
		sys.exit('usage: stats_edges.py <input_file_path>')


	# assume the graph is undirected, unweighted
	weighted = True
	directed = True

	# input_file_path = sys.argv[1]
	cur_file_path = str(Path().resolve().parent)

	# input_file_path = cur_file_path + '/link_prediction_graphs/yahoo-msg-heter/yahoo-msg-heter_wei_sample_0.6.tsv'
	# input_gt_path = cur_file_path + '/link_prediction_graphs/yahoo-msg-heter/yahoo-msg-heter_cat.tsv'
	# input_file_path = cur_file_path + '/link_prediction_graphs/citseer/citseer_wei_splitted_base_sample_0.6.tsv'
	# input_gt_path = cur_file_path + '/link_prediction_graphs/citseer/citseer_wei_splitted_cat.tsv'

	# input_file_path = cur_file_path + '/static_graphs/citseer_wei_splitted_base.tsv'
	# input_gt_path = cur_file_path + '/static_graphs/citseer_wei_splitted_cat.tsv'

	# input_file_path = cur_file_path + '/static_graphs/yahoo-msg-heter_wei_splitted.tsv'
	# input_gt_path = cur_file_path + '/static_graphs/yahoo-msg-heter_wei_splitted_cat.tsv'
	# input_file_path = cur_file_path + '/real_graphs/wiki-talk-temporal_temp_splitted.tsv'
	# input_gt_path = cur_file_path + '/real_graphs/wiki-talk-temporal_temp_splitted_cat.tsv'
	# input_file_path = cur_file_path + '/real_graphs/soc-sign-bitcoinotc_wei_temp_splitted.tsv'
	# input_gt_path = cur_file_path + '/real_graphs/soc-sign-bitcoinotc_wei_temp_splitted_cat.tsv'
	input_file_path = cur_file_path + '/real_graphs/soc-sign-bitcoinalpha_wei_temp_splitted_base.tsv'
	input_gt_path = cur_file_path + '/real_graphs/soc-sign-bitcoinalpha_wei_temp_splitted_cat.tsv'
	# input_file_path = cur_file_path + '/real_graphs/digg_wei_temp_undir_splitted.tsv'
	# input_gt_path = cur_file_path + '/real_graphs/digg_wei_temp_undir_splitted_cat.tsv'
	# input_file_path = cur_file_path + '/toy_graphs/yahoo-msg-heter_wei_splitted.tsv'
	# input_gt_path = cur_file_path + '/toy_graphs/yahoo-msg-heter_wei_splitted_cat.tsv'
	# input_gt_path = cur_file_path + '/toy_graphs/yahoo-msg-heter_cat.tsv'
	# input_file_path = cur_file_path + '/toy_graphs/test1_splitted.tsv'
	# input_file_path = cur_file_path + '/toy_graphs/test1_temporal.tsv'
	# input_gt_path = cur_file_path + '/toy_graphs/test1_splitted_cat.tsv'
	# input_file_path = cur_file_path + '/toy_graphs/test1.tsv'
	# input_gt_path = cur_file_path + '/toy_graphs/test1_cat.tsv'

	##############################
	# for exp
	##############################
	# input_file_path = cur_file_path + '/real_graphs/yahoo-msg-heter/yahoo-msg-heter_wei.tsv'
	# input_gt_path = cur_file_path + '/real_graphs/yahoo-msg-heter/yahoo-msg-heter_cat.tsv'
	# input_file_path = cur_file_path + '/exp/ER_100_500.tsv'
	# input_gt_path = cur_file_path + '/exp/ER_100_500_cat.tsv'
	# input_file_path = cur_file_path + '/exp/ER_1000_5000.tsv'
	# input_gt_path = cur_file_path + '/exp/ER_1000_5000_cat.tsv'
	# input_file_path = cur_file_path + '/exp/ER_10000_50000.tsv'
	# input_gt_path = cur_file_path + '/exp/ER_10000_50000_cat.tsv'
	# input_file_path = cur_file_path + '/exp/ER_100000_500000.tsv'
	# input_gt_path = cur_file_path + '/exp/ER_100000_500000_cat.tsv'
	# input_file_path = cur_file_path + '/exp/ER_1000000_5000000.tsv'
	# input_gt_path = cur_file_path + '/exp/ER_1000000_5000000_cat.tsv'
	# input_file_path = cur_file_path + '/exp/ER_10000000_50000000.tsv'
	# input_gt_path = cur_file_path + '/exp/ER_10000000_50000000_cat.tsv'
	##############################

	delimiter = " "
	if ".csv" in input_file_path:
		delimiter = ","
	elif ".tsv" in input_file_path:
		delimiter = "\t"
	else:
		sys.exit('Format not supported.')

	''' adj_matrix: lil format adj matrix
		edge_time_dict: (src, dst, wei) - time_1, time_2, ...
	'''
	check_eq, num_nodes, num_edges, adj_matrix, edge_time_dict, time_edge_dict, start_time, end_time = parse_weighted_temporal(input_file_path, delimiter)

	############################################################################################################
	# Setup
	############################################################################################################

	walks_num = 10
	walk_length = 20
	nodes_to_explore = range(num_nodes)
	base_features = ['degree', 'indegree', 'outdegree']

	K = 128
	num_buckets = 4	#2
	bucket_max_value = 30

	Ks = [K*1/3, K/3, K/3]
	Ts = [4, 4, 4]
	# Ks = [K/5, K/5, K/5, K/5, K/5]

	# Ks = [K/4, K/4, K/4, K/4]

	# Ks = [K/3, K/3, K/3]

	# Ks = [K/2, K/2]

	# Ks = [K]
	Ts = [4 for _ in range(len(Ks))]
	dist_scope = len(Ks)

	init_mod = 'early'
	walk_mod = 'early'

	graph_mod_external = 'static' # <Set this to 'static' only when we want to perform static analysis on dynamic graphs todo: remove this if necessary>

	############################################################################################################


	graph_mod = 'static' if edge_time_dict is None else 'dynamic'
	print '[Graph mode detected] ' + graph_mod
	graph_mod = graph_mod_external if graph_mod != graph_mod_external else graph_mod
	print '[Graph mode set as] ' + graph_mod


	SM = Static_Methods(adj_matrix = adj_matrix, nodes_to_explore = nodes_to_explore)
	DM = Dynamic_Methods(adj_matrix = adj_matrix, nodes_to_explore = nodes_to_explore, edge_time_dict = edge_time_dict)
	neighbor_list_static = SM.construct_neighbor_list()
	neighbor_list_dynamic = DM.construct_neighbor_list()


	CAT_DICT, ID_CAT_DICT = construct_cat(input_gt_path, delimiter)


	G = Graph(adj_matrix = adj_matrix, edge_time_dict = edge_time_dict, num_nodes = num_nodes, num_edges = num_edges, weighted = weighted, directed = directed, check_eq = check_eq, 
		start_time = start_time, end_time = end_time, num_buckets = num_buckets, bucket_max_value = bucket_max_value, cat_dict = CAT_DICT, id_cat_dict = ID_CAT_DICT, 
		neighbor_list_static = neighbor_list_static, neighbor_list_dynamic = neighbor_list_dynamic, time_edge_dict = time_edge_dict, dist_scope = dist_scope)

	walks, S_out = G.simulate_walks(walks_num, walk_length, nodes_to_explore, walk_mod, graph_mod, init_mod)
	
	fOut = open('walks.txt', 'w')
	for walk in walks:
		fOut.write(str(walk) + '\n')

	fOut.close()

	# print S_out
	# S_out = aggregate_walks(walks, walks_num, nodes_to_explore, dist_scope)
	sps.save_npz('./S.npz', sps.coo_matrix(S_out))

	PSs = construct_prox_structure(G, nodes_to_explore, base_features, S_out, dist_scope)
	# print '====='
	# print PSs

	# Ks = [K*1/4, K/4, K/4, K/4]
	# Ts = [4, 4, 4, 4]

	
	# Ks = [K*2/4, K*2/4]
	# Ts = [4, 4]

	tables, rep = hash_func(PSs, Ks, Ts)

	ifSVD = False

	if ifSVD:
		emb = svd_emb(PSs, Ks)
		np.savetxt('rep_svd.tsv', emb.round(6), fmt='%.6f', delimiter = '\t')

	print '-------------'
	np.savetxt('rep.tsv', rep, fmt='%i', delimiter = '\t')
	sps.save_npz('./rep.npz', sps.coo_matrix(rep))


	
	# table_1s = tables[0]
	# table_2s = tables[1]
	# print len(table_1s)
	# print len(table_2s)
	

	fOut = open('hashtable_1s.tsv', 'w')
	for key in tables[0]:
		key_trans = ''.join([str(ele) for ele in key])
		fOut.write(str(key_trans) + delimiter + str(tables[0][key]) + '\n')
	fOut.close()
	fOut = open('hashtable_2s.tsv', 'w')
	for key in tables[1]:
		key_trans = ''.join([str(ele) for ele in key])
		fOut.write(str(key_trans) + delimiter + str(tables[1][key]) + '\n')
	fOut.close()
	fOut = open('hashtable_3s.tsv', 'w')
	for key in tables[2]:
		key_trans = ''.join([str(ele) for ele in key])
		fOut.write(str(key_trans) + delimiter + str(tables[2][key]) + '\n')
	fOut.close()

	print rep[1]















