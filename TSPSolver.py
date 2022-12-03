#!/usr/bin/python3

from which_pyqt import PYQT_VER
import itertools
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def greedy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		count = 0
		start_time = time.time()
		distMatrix = np.full((ncities, ncities), float('inf'), dtype=float)

		for x in range(ncities):
			for y in range(ncities):
				distMatrix[x][y] = cities[x].costTo(cities[y])
		
		while time.time() - start_time < time_allowance:
			reducedMatrix = distMatrix.copy()
			# start at a random city
			city = random.randint(0, ncities-1)
			reducedMatrix[:, city] = float('inf')
			# start a new path
			path = [city]

			for i in range(ncities):
				# pick shortest path from current city
				m = np.min(reducedMatrix[city])
				# get index number of next city with the shortest path
				nextCity = 0
				for num in reducedMatrix[city]:
					if num == m:
						break
					nextCity += 1
				reducedMatrix[:, nextCity] = float('inf')
				# update path
				path.append(nextCity)
				city = nextCity

			count += 1

			# a full path has been found
			if len(cities) == ncities:
				bssf = TSPSolution([cities[x] for x in path])
				# a cost has been found
				if bssf.cost != float('inf'):
					break

		end_time = time.time()
		results['cost'] = bssf.cost
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def branchAndBound( self, time_allowance=60.0 ):
		class Node:
			def __init__(self, _lowBound, _costMatrix, _path, _city):
				self.lowBound = _lowBound
				self.costMatrix = _costMatrix
				self.path = _path
				self.city = _city

			def reduceMatrix(self):
				for row in range(len(self.costMatrix)):
					minCost = np.min(self.costMatrix[row])
					if minCost != float('inf') and minCost != 0:
						self.lowBound += minCost
						self.costMatrix[row] -= minCost
				for col in range(len(self.costMatrix)):
					minCost = np.min(self.costMatrix[:, col])
					if minCost != float('inf') and minCost != 0:
						self.lowBound += minCost
						self.costMatrix[:, col] -= minCost

			def getNewSubNode(self, dest):
				tempLowBound = self.lowBound + self.costMatrix[self.city][dest]

				tempCostMatrix = self.costMatrix.copy()
				tempCostMatrix[self.city] = float('inf')
				tempCostMatrix[:, dest] = float('inf')

				tempPath = self.path.copy()
				tempPath.append(self.city)

				temp = Node(tempLowBound, tempCostMatrix, tempPath, dest)
				temp.reduceMatrix()

				return temp

			def __gt__(self, other):
				if other.lowBound / len(other.path) < self.lowBound / len(self.path):
					return True
				else:
					return False

		time_allowance *= 10
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		count = 0
		bssf = self.defaultRandomTour()['soln']
		for i in range(9):
			minn = self.defaultRandomTour()['soln']
			if minn.cost < bssf.cost:
				bssf = minn
		nstates = 0
		npruned = 0
		queue = []
		maxqueue = 1
		start_time = time.time()

		costMatrix = np.full((ncities, ncities), float('inf'), dtype=float)

		# fill up the cost matrix for the root
		for x in range(ncities):
			for y in range(ncities):
				costMatrix[x][y] = cities[x].costTo(cities[y])

		root = Node(0, costMatrix, [], 0)
		nstates += 1

		# reduce matrix
		root.reduceMatrix()

		# create priority queue
		heapq.heappush(queue, root)

		while len(queue) != 0 and time.time() - start_time < time_allowance:
			maxqueue = max(maxqueue, len(queue))
			subNode = heapq.heappop(queue)

			for dest in range(ncities):
				if dest != subNode.city:
					nstates += 1
					# check to see if next cities are fruitful
					if subNode.costMatrix[subNode.city][dest] + subNode.lowBound <= bssf.cost:
						# a fruitful path, find new node from sub node
						nextNode = subNode.getNewSubNode(dest)
						if nextNode.lowBound < bssf.cost:
							if len(nextNode.path) == ncities:
								newSolution = TSPSolution([cities[i] for i in nextNode.path])
								if nextNode.lowBound < bssf.cost:
									count += 1
									bssf = newSolution
							else:
								heapq.heappush(queue, nextNode)
								maxqueue = max(maxqueue, len(queue))
						else:
							npruned += 1
					else:
						npruned += 1

		end_time = time.time()

		for node in queue:
			if node.lowBound > bssf.cost:
				npruned += 1

		results['cost'] = bssf.cost
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = maxqueue
		results['total'] = nstates
		results['pruned'] = npruned
		return results



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''

	def fancy(self, time_allowance=60.0):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)

		distMatrix = np.full((ncities, ncities), float('inf'), dtype=float)

		for x in range(ncities):
			for y in range(ncities):
				distMatrix[x][y] = cities[x].costTo(cities[y])

		# held-karp implementation
		# used this: https://github.com/CarlEkerot/held-karp/blob/master/held-karp.py
		costs = {}

		for i in range(1, len(distMatrix)):
			costs[(1 << i, i)] = (distMatrix[0][i], 0)

		start_time = time.time()
		count = 0
		for subSize in range(2, len(distMatrix)):
			if time.time() - start_time > time_allowance:
				break
			for subset in itertools.combinations(range(1, len(distMatrix)), subSize):
				bits = 0
				count += 1
				for bit in subset:
					bits |= 1 << bit
				for k in subset:
					prev = bits & ~(1 << k)

					res = []
					for m in subset:
						if m == 0 or m == k:
							continue
						res.append((costs[(prev, m)][0] + distMatrix[m][k], m))
					costs[(bits, k)] = min(res)
		bits = (2 ** len(distMatrix) - 1) - 1

		# Calculate optimal cost
		res = []
		for k in range(1, len(distMatrix)):
			res.append((costs[(bits, k)][0] + distMatrix[k][0], k))
		opt, parent = min(res)

		end_time = time.time()
		results['cost'] = opt
		results['time'] = end_time - start_time
		results['count'] = 1
		results['total'] = count
		results['soln'] = None
		return results




