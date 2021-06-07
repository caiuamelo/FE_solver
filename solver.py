__author__ = 'vskfv_000'

from scipy import *
import numpy

class Solver:
    nnd = 3 # Number of nodes:
    nel = 3 # Number of elements:
    nne = 2 # Number of nodes per element:
    nodof = 2 # Number of degrees of freedom per node
    eldof = nne*nodof # Number of degrees of freedom per element

    #Nodes coordinates X and Y
    geom = array([[0., 0.],       # X and Y coord. node 1
                 [4000., 0.],     # and Y coord. node 2
                 [4000.,6000.]])  # X and X coord. node 3

    #Element connectivity
    connec = array([[1, 2],   # 1st and 2nd node of element 1
                    [2, 3],   # 1st and 2nd node of element 2
                    [1, 3]])  # 1st and 2nd node of element 3


    # Geometrical properties
    # prop(1,1) = E; prop(1,2)= Area
    prop = array([[200000., 2300.],  # E and A of element 1
                  [200000., 2300.],  # E and A of element 2
                  [200000., 2300.]]) # E and A of element 3

    # Boundary conditions
    nf = ones((nnd, nodof)) # Initialize the matrix nf to 1
    nf[0][0] = 0 # Prescribed nodal freedom of node 1
    nf[0][1] = 0
    nf[1][1] = 0 # Prescribed nodal freedom of node 3

    # Degrees of freedom
    n = 0

    # loading
    load = array([[0., 0.],
                  [0., 0.],
                  [12000., 0.]])  # forces in X and Y directions at node 3

    # The global stiffness matrix
    KK = array([])

    # The global force vector
    F = array([])

    # The displacement vector
    node_disp = array([])

    def __init__(self):
        n = 0
        for i in arange(0, self.nnd):
            for j in arange(0, self.nodof):
                if self.nf[i][j] != 0:
                    n += 1
                    self.nf[i][j] = n
        self.n = n

        self.KK = zeros((self.n,self.n))  # Initialize global stiffness matrix to zero
        self.F = zeros((self.n,1)) # Initialize global force vector
        self.node_disp = zeros((self.nnd,self.nodof))


    def truss_g(self, i):  # This method forms the steering vector for element i

        # retrieve the nodes of element i
        node_1 = self.connec[i][0]
        node_2 = self.connec[i][1]

        # Form the steering vector from elements degrees of freedom

        g = array([[self.nf[node_1-1][0]],
                   [self.nf[node_1-1][1]],
                   [self.nf[node_2-1][0]],
                   [self.nf[node_2-1][1]]])

        return g

    def truss_kl(self, i):  # This method forms the element stiffness matrix in local coordinates

        # retrieve the nodes of element i
        node_1 = self.connec[i][0]
        node_2 = self.connec[i][1]

        # Retrieve the x and y coordinates of nodes 1 and 2
        x1 = self.geom[node_1-1][0]
        y1 = self.geom[node_1-1][1]
        x2 = self.geom[node_2-1][0]
        y2 = self.geom[node_2-1][1]

        # Evaluate length of element i
        L = sqrt((x2-x1)**2 + (y2-y1)**2)

        # Retrieve section properties of element i
        E = self.prop[i][0]
        A = self.prop[i][1]

        # Calculate element stiffness matrix in its local coordinates

        kl = array([[E*A/L, 0., -E*A/L, 0.],
                    [0., 0., 0., 0.],
                    [-E*A/L, 0., E*A/L, 0.],
                    [0., 0., 0., 0.]])

        return kl


    def truss_c(self, i):  # This method forms the transformation between local and global coordinates

        node_1 = self.connec[i][0]
        node_2 = self.connec[i][1]

        # Retrieve the x and y coordinates of nodes 1 and 2
        x1 = self.geom[node_1-1][0]
        y1 = self.geom[node_1-1][1]

        x2 = self.geom[node_2-1][0]
        y2 = self.geom[node_2-1][1]

        # Evaluate the angle that the member makes with the global axis X
        if (x2-x1) == 0:
            if (y2 > y1):
                theta = numpy.pi/2
            else:
                theta = -numpy.pi/2

        else:
            theta=math.atan((y2-y1)/(x2-x1))

        # Construct the transformation matrix
        #
        C = array([[numpy.cos(theta), -numpy.sin(theta), 0., 0. ],
                   [numpy.sin(theta), numpy.cos(theta), 0., 0. ],
                   [0., 0., numpy.cos(theta), -numpy.sin(theta)],
                   [0., 0., numpy.sin(theta), numpy.cos(theta) ]])


        return C
        
    def form_kk(self, g=array([]), kg=array([])):  # This method assembles the global stiffness matrix

        for i in arange(0,self.eldof):
            if g[i] != 0:
                for j in arange(0, self.eldof):
                    if g[j] != 0:
                       self.KK[int(g[i]-1)][int(g[j]-1)] += kg[i][j]


    def form_truss_f(self):  # This function forms the global force vector

        force = zeros((self.n,1))
        for i in arange(0, self.nnd):
            for j in arange(0, self.nodof):
                if self.nf[i][j] != 0:
                    force[self.nf[i][j]-1][0] = self.load[i][j]

        self.F = force


    def print_results(self): # This method prints the results of the solved problem (Displacements and forces)
        pass

    def run(self):  # Run the model

        for i in arange(0, self.nel):

            kl = self.truss_kl(i)  # Form element matrix in local xy
            C = self.truss_c(i)  # Form transformation matrix

            kg = C.dot(kl.dot((C.transpose())))  # Transform the element matrix from local to global coordinates

            g = self.truss_g(i)  # Retrieve the element steering vector
            self.form_kk(g, kg)  # assemble global stiffness  matrix
            self.form_truss_f()  # assemble global force vector

        # End of assembly
        delta = numpy.linalg.solve(self.KK,self.F)  # solve for unknown displacements

        print delta

        # Extract nodal displacements

        for i in arange(0, self.nnd):
            for j in arange(0, self.nodof):
                if self.nf[i][j]!= 0:
                    self.node_disp[i][j] = delta[int(self.nf[i][j])-1]

        # Calculate the forces acting on each element
        # in local coordinates, and store them in the
        # vector force().

        edg = zeros((self.eldof,1))

        for i in arange(0,self.nel):

            kl = self.truss_kl(i) # Form element matrix in local xy
            C = self.truss_c(i)   # Form transformation matrix
            kg = C.dot(kl.dot(C.transpose()))  # Transform the element matrix from
            # local to global coordinates
            g = self.truss_g(i)  # Retrieve the element steering vector
            for j in arange(0,self.eldof):
                if g[j] == 0:
                    edg[j] = 0. # displacement = 0. for restrained freedom
                else:
                    edg[j] = delta[int(g[j]-1)]

        fg = kg.dot(edg) # Element force vector in global XY
        fl = (C.transpose()).dot(fg) # Element force vector in local xy



a = Solver()
a.run()
#print a.node_disp[2][1]
