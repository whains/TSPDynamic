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

		while True:
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
		pass



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
		start_time = time.time()

		distMatrix = np.full((ncities, ncities), float('inf'), dtype=float)

		for x in range(ncities):
			for y in range(ncities):
				distMatrix[x][y] = cities[x].costTo(cities[y])

		min_distance = self.tsp_distance(distMatrix, ncities)

		end_time = time.time()
		results['cost'] = min_distance
		results['time'] = end_time - start_time
		results['count'] = 0
		results['soln'] = None
		return results

	# I copied this from a YouTube video and really have no clue how it works. Needs more edits and to make more sense lol
	def tsp_distance(self, distMatrix, ncities):
		cost = {(1 << dest, dest): distMatrix[0][dest] for dest in range(1, ncities)}  # (subset, endpoint)
		for size in range(2, ncities):
			for sub in itertools.combinations(range(1, ncities), size):
				sub_i = sum([1 << x for x in sub])
				for dest in sub:
					cost[sub_i, dest] = min([cost[sub_i & ~ (1 << dest), side] + distMatrix[side][dest] for side in sub if side != dest])
		return min([cost[(1 << ncities) - 2, final] + distMatrix[final][0] for final in range(1, ncities)])


