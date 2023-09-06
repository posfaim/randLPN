#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools


#=============Network generators

def LinearPhysNetModel(nodes, d=0, L=10, node_link_int=False):
    '''
    Generates an instance of a random linear physical network. For details see https://arxiv.org/abs/2211.13265
    
    Parameters
    ----------
    nodes : ndarray
        A numpy array with shape (N,3) containing the coordinates of the center of each node. 
    d : float
        Diameter of links and nodes.
    L : int
        Maximum number of links added to the network. The number of physical links that are possible to add might be smaller than L.
    node_link_int : bool
        If True, include node-link volume exclusion into the model.
    
    Returns
    -------
    vs : ndarray
        Array of integers containing the index of the start point of each physical link.
    ws : ndarray
        Array of integers containing the index of the end point of each physical link.
    '''
    
    N = len(nodes)  #number of nodes
    
    #randomize order of link candidates
    candidate_order = list(itertools.combinations(range(len(nodes)),2))
    np.random.shuffle(candidate_order)
    
    vs=np.zeros(L,dtype=int)    #start point of successfully added links
    ws=np.zeros(L,dtype=int)    #end point of successfully added links
    added_links = np.zeros((L,6))
    l=0 #number of successfully added links

    for v,w in candidate_order:
        #get array representing new link
        new_link = np.concatenate([nodes[v],nodes[w]])
        
        # test for node-link overlap
        node_link_overlap=False
        if node_link_int:
            #get distances from all points
            point_dists = dist3D_Point_to_Segment(nodes,new_link)
            # check if link axis is sufficiently far away from point centers
            if (point_dists<=d).sum()>2: #this is 2 because endpoints overlap with the link
                node_link_overlap = True
        
        if not node_link_overlap:
            # get distances from all added links
            if l>0:
                link_dists = np.apply_along_axis(lambda S: dist3D_Segment_to_Segment(S,new_link),1,added_links[:l])
                
                # links can overlap if they share endpoints
                # check if endpoint is the same as added link
                endpoint_overlap = (vs[:l]==v) | (vs[:l]==w) | (ws[:l]==v) | (ws[:l]==w)
                

                # add link if no overlap
                link_link_overlap = False
                if ((link_dists>d) | endpoint_overlap).all():
                    link_link_overlap = True
                    
            if l==0 or link_link_overlap:
                vs[l]=v
                ws[l]=w
                added_links[l] = new_link
                l+=1

        if l==L:
            break
            
    return vs[:l],ws[:l]


def PlaceRandPoints(N,d, max_num_trials=10000):
	'''
	Randomly places points in the unit cube, such that no points are closer than d to each other
	
	Parameters
    ----------
	N : int
	    Number of points to be placed.
	d : float
	    Minimum distance between points.
	max_num_trials: int
	    Maximum number of times we attempt to place a point.
	    
	Returns
    -------
	rand_points : ndarray or None
	    Array with shape (N,3) containing the coordinates of the points or None if we failed to place all nodes in a given number of trials.
	'''
	
	rand_points =  np.zeros((N,3))
	rand_points[0] =  np.random.uniform(size=3) #place first point
	for i in range(1,N):
		count = 0 #number of trials
		
		new_point = np.random.uniform(size=3) # random location for new point
		while (np.linalg.norm(rand_points-new_point,axis=1)<d).any(): #test for distance
			new_point = np.random.uniform(size=3) #select new random location
			count +=1
			if count==max_num_trials:
			    #Trials exceed threshold, give up
				print("Oh boy, this is taking too long...",file=sys.stderr)
				return None
				
		#store point coordinate
		rand_points[i]=new_point

	return rand_points


def dist3D_Point_to_Segment(P,S):
    '''
    Distance between a point or a collection of points and a segment.
    Based on http://geomalgorithms.com/a07-_distance.html
	
	Parameters
    ----------
	P : ndarray
        Coordinates of a point as a numpy array with shape (3,) or (3,N) for N points
	S : ndarray
	    Coordinates of endpoints of a segment as a numpy array with shape (6,), i.e., [x1,y1,z1,x2,y2,z2]
	    
	Returns
    -------
	
	dist : ndarray or float
	    Distance(s) between point(s) and segment.
	'''
   
    A = S[:3]
    B = S[3:]
    AB = A-B
    e  = AB/np.linalg.norm(AB)
    # let Q be the point on AB line that is perpendicular to P
    QA = ((P-A)*e).sum(axis=1) #negative if outside AB
    QB = ((B-P)*e).sum(axis=1)

    QX = np.max(np.stack((QA,QB,np.zeros(len(P)))),axis=0)
    
    PQ = np.apply_along_axis(np.linalg.norm,1,np.cross(P-A,e))
    
    dist = np.hypot(PQ,QX)
    
    return dist


def dist3D_Segment_to_Segment(S1,S2):
    '''
    Distance between two segments.
    Based on http://geomalgorithms.com/a07-_distance.html
	
	Parameters
    ----------
	S1 : ndarray
	    Coordinates of endpoints of a segment as a numpy array with shape (6,), i.e., [x1,y1,z1,x2,y2,z2]
	S2 : ndarray
	    Coordinates of endpoints of a segment as a numpy array with shape (6,), i.e., [x1,y1,z1,x2,y2,z2]
	
	Returns
    -------
	
	dist : ndarray or float
	    Distance(s) between point(s) and segment.
	'''
	
    '''
    S1: endpoints of segment as a single numpy array with shape (6,), i.e., [x1,y1,z1,x2,y2,z2]
    S2: endpoints of segment as a single numpy array with shape (6,), i.e., [x1,y1,z1,x2,y2,z2]
    returns the distance between the two segments
    
    Based on http://geomalgorithms.com/a07-_distance.html
    '''
    
    SMALL_NUM=1e-6
    
    S1P0 = S1[:3]
    S1P1 = S1[3:]
    
    S2P0 = S2[:3]
    S2P1 = S2[3:]
    
    u = S1P1 - S1P0
    v = S2P1 - S2P0
    w = S1P0 - S2P0
    a = np.dot(u,u)         # always >= 0
    b = np.dot(u,v)
    c = np.dot(v,v)         # always >= 0
    d = np.dot(u,w)
    e = np.dot(v,w)
    D = a*c - b*b;        # always >= 0
    sD = D
    tD = D

    # compute the line parameters of the two closest points
    if D < SMALL_NUM: #{ // the lines are almost parallel
        sN = 0.0     #    // force using point P0 on segment S1
        sD = 1.0     #    // to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:                 # get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if sN < 0.0:       # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD: # sc > 1  => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0:            # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a
    elif tN > tD:      # tc > 1  => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if -d + b < 0.0:
            sN = 0
        elif (-d + b) > a:
            sN = sD
        else:
            sN = (-d +  b)
            sD = a
        
    # finally do the division to get sc and tc
    if np.abs(sN) < SMALL_NUM:
        sc = 0.
    else:
        sc= sN/sD
    if np.abs(tN) < SMALL_NUM:
        tc = 0.
    else:
        tc = tN/tD
    
    # get the difference of the two closest points
    
    dP = w + (sc * u) - (tc * v)  # =  S1(sc) - S2(tc)

    dist = np.linalg.norm(dP)

    return dist    # return the closest distance


def DrawLinearPhysNet(nodes,vs,ws):
	'''
	Draws a linear physical network using matplotlib.
	Colors of the links correspond to the position of its half point (red = x, green = y, blue = z)
	
	Input
	points:		A numpy array with shape (N,3) containing the coordinates of the points
	vs, ws:		Two integer arrays containing the indices of the endpoints of the links
	'''
	
	'''
    Draws a linear physical network using matplotlib.
	Colors of the links correspond to the position of its half point (red = x, green = y, blue = z)
	
	Parameters
    ----------
    nodes : ndarray
        A numpy array with shape (N,3) containing the coordinates of the center of each node.
	vs : ndarray
        Array of integers containing the index of the start point of each physical link.
    ws : ndarray
        Array of integers containing the index of the end point of each physical link.
	
	Returns
    -------
	ax : matplotlib.axes
	    Axes of plot.
	'''
	
	fig = plt.figure(figsize=(10,10))
	ax = fig.add_subplot(111, projection='3d')

	for v,w in zip(vs,ws):
	   color = (nodes[v]+nodes[w])/2.
	   ax.plot(*(nodes[[v,w]].T),color=color)

	ax.set_xlim([0,1])
	ax.set_ylim([0,1])
	ax.set_zlim([0,1])
	ax.axis('off')
	plt.subplots_adjust(left=0., right=1, top=1, bottom=0)

	return ax


def main():
    #Example
    
    #Generate network
    L = 100
    N = 20
    d =.05
    
    nodes = PlaceRandPoints(N, d)
    vs, ws = LinearPhysNetModel(nodes, d, L, node_link_int=True)
    
    DrawLinearPhysNet(nodes,vs,ws)
    plt.show()
    
    '''
    #Create igraph object representing the abstract network
    import igraph as ig
    net = ig.Graph(N)
    net.add_edges(list(zip(vs,ws)))

    #Set node size and color for plotting
    net.vs['size']=10*np.log(np.array(net.degree())+1)+2
    net.vs['color'] = ['#%02x%02x%02x' % tuple(int(x*256) for x in rp) for rp in nodes]

    #Plot abstract network
    ig.plot(net)
    '''
    
    return
    
    
if __name__ == '__main__':
    main()
    
