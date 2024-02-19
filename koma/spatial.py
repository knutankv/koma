'''
MODULE FOR SPATIAL MAPPING AND VISUALIZATION
'''
import pyvista as pv
import pyvistaqt as pvqt
import numpy as np
import sys

from copy import deepcopy

def nodes_from_matrix(node_matrix):
    return [Node(*row) for row in node_matrix]

def elements_from_matrix(element_matrix):
    return [Element(el[1:], label=el[0]) for el in element_matrix]

class Element:
    '''
    Element class for creation of elements in model. Elements can have 2, 3 or 4 nodes, which 
    dictate what type of element it representes (line, triangle, rectangle).

    Arguments
    ------------
    nodes : `Node`
        list of `Node` objects to connect with element
    label : None, optional
        requested label given to element

    '''
    def __init__(self, nodes, label=None):
        self.nodes = nodes
        self.label = label

    def get_cog(self, deformed=False):
        if deformed:
            return np.mean([node.xyz for node in self.nodes], axis=0)
        else:
            return np.mean([node.xyz0 for node in self.nodes], axis=0)

    @property
    def num_points(self):
        return len(self.nodes)

    @property
    def nodelabels(self):
        return np.array([n.label for n in self.nodes]).astype(int)
    
    @property
    def padded_conn(self):
        return np.hstack([self.num_points, self.nodelabels])

    @property
    def el_type(self):
        if self.num_points == 4:
            return 'rectangle'
        elif self.num_points == 3:
            return 'triangle'
        elif self.num_points == 2:
            return 'line'
    
    # CORE METHODS
    def __str__(self):
        if self.label is not None:
            return f'Element ({self.el_type}) {self.label}'
        else:
            return f'Element ({self.el_type}) <{self.__hash__}>'  
        
    def __repr__(self):
        return self.__str__()
    
    def __eq__(self, other):
        if isinstance(other, Element):
            return self.label == other.label
        elif isinstance(other, int):
            return self.label == other
        
class Node:
    '''
    Node class for creation of nodes to use in model. Nodes can be defined by two or three coordinates.
    If the latter, z is assumed 0. Label is required.

    Arguments
    ------------
    label : int
        requested label given to node
    x : float
        x-coordinate of node
    y : float
        y-coordinate of node
    z : 0.0, optional
        float to specify z-coordinate of node

    '''

    def __init__(self, label, x, y, z=0.0):
        self.x0 = x
        self.y0 = y
        self.z0 = z

        self.x = self.x0
        self.y = self.y0
        self.z = self.z0


        self._label = label
        
    @property
    def label(self):
        return self._label
    @label.setter
    def label(self, val):
        self._label = val
    
    @property
    def xyz(self):
        return np.array([self.x, self.y, self.z])
    
    @property
    def xyz0(self):
        return np.array([self.x0, self.y0, self.z0])

    # CORE METHODS
    def __str__(self):
        if self.label is not None:
            return f'Node {self.label}'
        else:
            return f'Node <{self.__hash__}>'
        
    def __repr__(self):
        return self.__str__()
    
    def __eq__(self, other):
        if isinstance(other, Node):
            return self.label == other.label
        elif isinstance(other, int):
            return self.label == other

def rel3(master):
    return [Rel(master, 0), Rel(master, 1), Rel(master, 2)]

class Rel:
    def __init__(self, master=None, dof=None, fun=None):
        '''
        Class to generate relative constraints used for more complex `dofmap` creation.
        These are placed in relevant entries of the model `dofmap`.

        Parameters
        -------------
        master : None, int, optional
            node label of master node (not input if fun is given)
        dof : None, int, optional
            dof index of master node dof (not input if fun is given)
        fun : Non, fun, optional
            function to create more complex relationships

        Examples
        -------------
        `Rel(5, 1)` instructs to assign displacements in the chosen model DOF <--- node 5, dof index 1       
        `Rel(fun=lambda n: n(1,0)*0.5 - n(2,1)*0.2) instructs to assign displacements in the chosen model DOF from sum:
                                                    node 1 dof index 0 (scaled 0.5) - node 2 dof index 0.2 (scaled 0.2)

    
        '''

        self.dof = dof          # master dof
        self.master = master    # master node label
        self._fun = fun


    @property
    def fun(self):
        if self._fun is None:
            def f(n):
                return n(self.master, self.dof)*1.0
            return f
        else:
            return self._fun
        
    @fun.setter
    def fun(self, val):
        self._fun = val

        
class Model:
    '''
    Model class for creation of models for visualizing mode shapes.

    Arguments
    ------------
    nodes : `Node` 
        list of nodes
    elements : `Line` or `Triangle`
        list of elements (either lines or triangle patches)
    dofmap : dict, optional
        dictionary to map DOFs of model to DOFs of displacements
    sensors : dict, optional
    u : None, optional
        displacements (with dimensions corresponding to index values provided in dofmap)

    Example
    ------------
    Let's create a beam with an assumed sensor measuring vertical accelerations at midspan:
    __________v__________
    ^                   ^
    1         2         3    < --- Node labels

    The beam is 10 meters long, giving us these three nodes:

        nodes = [Node(1, 0, 0, 0), Node(2, 5, 0, 0), Node(3, 10, 0, 0)]

    Our elements are simply connecting node 1 to 2 and node 2 to 3:   
        
        elements = [Element([1,2])]
    
    Our dofmap would then be given as follows (assuming one DOF at midspan, along the z-direction):

        dofmap = {2: [0, 0, 1]}

    Let's go ahead and create our model object:
        
        model = Model(nodes, elements, dofmap=dofmap, sensors={'acc_z': 2})

    The last (optional) input, the `sensors` dictionary simply tells the model at which nodes our sensors are placed, for convenient
    plotting later.

    A more comprehensive example is provided in the Examples folder on GitHub.

    '''
    def __init__(self, nodes, elements, dofmap=None, sensors={}, u=None):
        self.nodes = nodes
        self.elements = elements
        self.dofmap = dofmap
        self.sensors = sensors
        self.u = u
    
        # Ensure node objects and not labels
        for el in self.elements:
            el.nodes = [self.get_node(n) for n in el.nodes] 

    @property 
    def u(self):
        return self._u
    
    @u.setter
    def u(self, val):
        ufull = self.expand_field(val)
        self.deform(ufull)
        self._u = val
    
    def deform(self, ufull):
        '''
        Deform system based on input full field displacements, ufull.
        '''

        for node in self.nodes:
            dxyz = ufull[self.get_node_ixs(node)]
            node.x, node.y, node.z = dxyz[0]+node.x0, dxyz[1]+node.y0, dxyz[2]+node.z0

    def expand_field(self, u):
        '''
        Exapend field of subsystem u to full field based on given dofmap.

        Arguments
        ----------
        u : float
            displacements to assign to system
            with dimensions corresponding to index values provided in dofmap
        
        Returns
        ---------
        ufull : float
            displacements corresponding to full field
        '''
        ufull = np.zeros(len(self.nodes)*3)

        if u is not None:
            # First, assign all direct (if not given - stays zero)
            for node in self.dofmap:
                rels = self.dofmap[node]
                for dof_ix, rel in enumerate(rels):
                    if type(rel) == int:
                        global_ix = self.get_node_ixs(node)[dof_ix]
                        if rel is not None:
                            ufull[global_ix] = u[rel]

            # Second, assign all relative
            n = self.grab_fun(ufull)
            for node in self.dofmap:
                rels = self.dofmap[node]
                for dof_ix, rel in enumerate(rels):
                    if type(rel) == Rel:
                        global_ix = self.get_node_ixs(node)[dof_ix]
                        ufull[global_ix] = rel.fun(n)

        return ufull


    def grab_fun(self, ufull):
        def n(node, dof):
            return ufull[self.get_node_ixs(node)[dof]]
        
        return n

    def get_elements(self, elements):
        if elements is None:
            return self.elements
        else:
            return [self.get_element(el) for el in elements]
    
    def get_element(self, element):
        return self.elements[self.get_element_index(element)]
    
    def get_element_index(self, element):
        return self.elements.index(element)

    def get_nodes(self, nodes):
        if nodes is None:
            return self.nodes
        else:
            return [self.get_node(node) for node in nodes]
    
    def get_node(self, node):
        return self.nodes[self.get_node_index(node)]
    
    def get_node_index(self, node):
        # Accepts both label and object
        return self.nodes.index(node)
    
    def get_node_ixs(self, node):
        return self.get_node_index(node)*3 + np.array([0,1,2])
    
    def get_node_indices(self, element):
        return np.array([self.get_node_index(node) for node in element.nodes])

    # Establish inputs for pyvista
    def get_points(self, nodes=None, deformed=False, flattened=True):
        nodes = self.get_nodes(nodes)
        
        if deformed:
            out = [n.xyz for n in self.nodes]
        else:
            out = [n.xyz0 for n in self.nodes]
        
        if flattened:
            return np.hstack(out)
        else:
            return np.vstack(out)
    

    def get_lines(self):
        lines = np.array([])
        
        for el in self.get_line_elements():
            lines = np.hstack([lines, el.num_points, self.get_node_indices(el)]).astype(int)
                     
        return lines
    
    def get_faces(self):
        faces = np.array([])
        
        for el in self.get_face_elements():
            faces = np.hstack([faces, el.num_points, self.get_node_indices(el)]).astype(int)
                               
        return faces

    # Grab elements based on type
    def filter_elements(self, allowed_types=['line', 'rectangle', 'triangle']):
        return [el for el in self.elements if el.el_type in allowed_types]
    
    def get_element_cogs(self, elements=None, deformed=False, flattened=False):
        elements = self.get_elements(elements)
        cogs = [el.get_cog(deformed=deformed) for el in elements]

        if flattened:
            return np.hstack(cogs)
        else:
            return np.vstack(cogs)
    
    def get_face_elements(self):
        return self.filter_elements(allowed_types=['rectangle', 'triangle'])
        
    def get_line_elements(self):
            return self.filter_elements(allowed_types=['line'])
        
    @property
    def n_faces(self):
        return len(self.get_elements(filter_elements=['rectangle', 'triangle']))

    def plot(self, pl=None, show=True, plot_lines=True, plot_nodes=True, 
                node_labels=False, element_labels=False,
                plot_faces=True, canvas_settings={}, 
                node_settings={}, line_settings={}, 
                face_settings={}, nodelabel_settings={}, sensorlabel_settings={},
                elementlabel_settings={}, view=None, sensor_labels=False,
                deformed=True, 
                node_label_fun=None, element_label_fun=None, perspective_cam=False, background_plotter=True):


        '''
        Plot model (either undeformed or deformed).

        Arguments
        ----------
        pl : `pyvista.Plotter()` object, optional
            plotter object - one will be created if not input
        show : True, optional
            whether or not to show plot at end - output is pl-object, which can be shown later
        plot_lines : True, optional
            whether or not to plot the line elements
        plot_nodes : True, optional
            whether or not to plot the nodes
        node_labels : True, optional
            whether or not to show node labels
            if input is a list of node labels, only these are shown
        element_labels :True, optional
            whether or not to show element labels
            if input is a list of element labels, only these are shown
        plot_faces : True, optional
            whether or not to plot face elements (triangles and rectangles)
        canvas_settings : dict, optional
            dictionary with additional settings for canvas
        node_settings : dict, optional
            dictionary with additional settings for undeformed nodes / points
            given to `pyvista.Plotter().add_points`
        line_settings : dict, optional
            dictionary with additional settings for undeformed line elements
            given to `pyvista.Plotter().add_mesh`
        face_settings : dict, optional
            dictionary with additional settings for undeformed patch elements
            given to `pyvista.Plotter().add_mesh`
        nodelabel_settings : dict, optional
            dictionary with additional settings for node labels
            given to `pyvista.Plotter().add_point_labels`
        sensorlabel_settings : dict, optional
            dictionary with additional settings for sensor labels 
            (=nodes with defined additional label)
            given to `pyvista.Plotter().add_point_labels`
        elementlabel_settings : dict, optional
            dictionary with additional settings for element labels
            given to `pyvista.Plotter().add_point_labels`
        view : {None, 'xy', 'yx', 'xz', 'top', 'side', 'front'}, optional
            string defining view, if None standard isometric perspective plot is given
        sensor_labels : False, optional
            whether or not to show sensor labels
        plot_states : ['deformed', 'undeformed'], optional
            which states of the system to plot
        label_priority : str, optional 
        node_label_fun : fun, optional
            function to define what property from each node is given in the node label strings
            if None, fun = lambda node: node.label
        element_label_fun : fun, optional
            function to define what property from each element is given in the element label strings
            if None, fun = lambda el: el.label
        perspective_cam : False, optional
            whether or not to use perspective projection (parallel projection otherwise)
        background_plotter : True, optional
            whether or not to plot in background 
            (meaning that Python processes in terminal can run while plot is open) 
            if True, pyvistaqt is used to generate plot (different interface)
     
        Returns
        ----------
        pl : `pyvista.Plotter()` object
        
        Example
        ---------
        See separate notebook on GitHub for example.

        '''

        if node_label_fun is None:
            node_label_fun = lambda n: str(int(n.label))
        if element_label_fun is None:
            element_label_fun = lambda e: str(int(e.label))

        # Element label settings
        elementlabel_settings = dict(always_visible=True, 
                                      text_color='blue', 
                                      shape_color='white', 
                                      shape_opacity=0.2) | elementlabel_settings

        # Node and sensor label settings
        nodelabel_settings = dict(always_visible=True, shape_opacity=0.2)
        sensorlabel_settings = nodelabel_settings | sensorlabel_settings

        # Node plotting settings
        node_settings = dict(
                            render_points_as_spheres=True,
                            color='black',
                            lighting=True,
                            point_size=6) | node_settings  

        # Line plotting settings
        line_settings = dict(render_lines_as_tubes=True,
                                style='wireframe',
                                lighting=True,
                                line_width=4) | line_settings
        
        canvas_settings = dict(background_color='white') | canvas_settings

        # Face plotting settings
        face_settings = dict(show_edges=True, color='lightgray') | face_settings

        if plot_nodes is not False:
            if plot_nodes is True:
                nodes_to_plot = self.nodes*1
            else:
                nodes_to_plot = [node for node in self.nodes if node in plot_nodes]

        if pl is None:
            if background_plotter:
                pl = pvqt.BackgroundPlotter()
            else:
                pl = pv.Plotter()

            for key in canvas_settings:
                setattr(pl, key, canvas_settings[key])

            if view is not None:
                if view in ['xy', 'top']:
                    pl.view_xy()
                if view in ['yz', 'front']:
                    pl.view_yz()
                if view in ['xz', 'side']:
                    pl.view_xz()

        if node_labels is not False:
            if node_labels is not True:
                labeled_nodes = [node for node in self.nodes if node in node_labels]
            else:
                labeled_nodes = self.nodes

            lbl = [node_label_fun(node) for node in labeled_nodes]
            pl.add_point_labels(self.get_points(nodes=labeled_nodes, deformed=deformed, flattened=False), lbl, **nodelabel_settings)
      

        if element_labels is not False:
            if element_labels is not True:
                labeled_els = [el for el in self.elements if el in element_labels]
            else:
                labeled_els = self.elements

            lbl = [node_label_fun(el) for el in labeled_els]
            pl.add_point_labels(self.get_element_cogs(elements=labeled_els, deformed=deformed, flattened=False), lbl, **elementlabel_settings)

                    
        if sensor_labels:
            for sensor in self.sensors:
                node = self.get_node(self.sensors[sensor])
                if deformed: 
                   pl.add_point_labels(np.array([node.xyz]), [sensor], **sensorlabel_settings)     
                else:
                   pl.add_point_labels(np.array([node.xyz0]), [sensor], **sensorlabel_settings)     

        points = pv.pyvista_ndarray(self.get_points(deformed=deformed, flattened=False))

        if plot_lines:
            lines = self.get_lines()
            self.line_mesh = pv.PolyData(points, 
                        lines=lines)        
            pl.add_mesh( self.line_mesh, **line_settings)

        if plot_faces:
            faces = self.get_faces()
            self.face_mesh = pv.PolyData(points,
                        faces=faces)
        
            pl.add_mesh(self.face_mesh, **face_settings)

        if plot_nodes:
            node_points = pv.pyvista_ndarray(self.get_points(nodes=nodes_to_plot, 
                                                        deformed=deformed, 
                                                        flattened=False))
            pts = pl.add_points(node_points, **node_settings)
            self.point_mesh = pts.mapper.dataset
        
        pl.camera.SetParallelProjection(not perspective_cam)

        if show:
            pl.show()

        return pl



    def plot_all(self, pl=None, show=True, plot_lines=True, plot_nodes=True, 
             node_labels=False, element_labels=False,
             plot_faces=True, canvas_settings={}, 
             node_settings={}, deformed_node_settings={},
             line_settings={}, deformed_line_settings={}, 
             face_settings={}, deformed_face_settings={}, 
             nodelabel_settings={}, sensorlabel_settings={},
             elementlabel_settings={}, view=None, sensor_labels=False,
             plot_states=['deformed','undeformed'], label_priority='deformed',
             node_label_fun=None, element_label_fun=None, perspective_cam=False, background_plotter=True):


        '''
        Plot model (undeformed or deformed or both).

        Arguments
        ----------
        pl : `pyvista.Plotter()` object, optional
            plotter object - one will be created if not input
        show : True, optional
            whether or not to show plot at end - output is pl-object, which can be shown later
        plot_lines : True, optional
            whether or not to plot the line elements
        plot_nodes : True, optional
            whether or not to plot the nodes
        node_labels : True, optional
            whether or not to show node labels
            if input is a list of node labels, only these are shown
        element_labels :True, optional
            whether or not to show element labels
            if input is a list of element labels, only these are shown
        plot_faces : True, optional
            whether or not to plot face elements (triangles and rectangles)
        canvas_settings : dict, optional
            dictionary with additional settings for canvas
        node_settings : dict, optional
            dictionary with additional settings for undeformed nodes / points
            given to `pyvista.Plotter().add_points`
        deformed_node_settings : dict, optional
            dictionary with additional settings for deformed nodes / points
            given to `pyvista.Plotter().add_points`
        line_settings : dict, optional
            dictionary with additional settings for undeformed line elements
            given to `pyvista.Plotter().add_mesh`
        deformed_line_settings : dict, optional
            dictionary with additional settings for deformed line elements
            given to `pyvista.Plotter().add_mesh`
        face_settings : dict, optional
            dictionary with additional settings for undeformed patch elements
            given to `pyvista.Plotter().add_mesh`
        deformed_face_settings : dict, optional
            dictionary with additional settings for deformed patch elements
            given to `pyvista.Plotter().add_mesh`
        nodelabel_settings : dict, optional
            dictionary with additional settings for node labels
            given to `pyvista.Plotter().add_point_labels`
        sensorlabel_settings : dict, optional
            dictionary with additional settings for sensor labels 
            (=nodes with defined additional label)
            given to `pyvista.Plotter().add_point_labels`
        elementlabel_settings : dict, optional
            dictionary with additional settings for element labels
            given to `pyvista.Plotter().add_point_labels`
        view : {None, 'xy', 'yx', 'xz', 'top', 'side', 'front'}, optional
            string defining view, if None standard isometric perspective plot is given
        sensor_labels : False, optional
            whether or not to show sensor labels
        plot_states : ['deformed', 'undeformed'], optional
            which states of the system to plot
        label_priority : str, optional 
        node_label_fun : fun, optional
            function to define what property from each node is given in the node label strings
            if None, fun = lambda node: node.label
        element_label_fun : fun, optional
            function to define what property from each element is given in the element label strings
            if None, fun = lambda el: el.label
        perspective_cam : False, optional
            whether or not to use perspective projection (parallel projection otherwise)
        background_plotter : True, optional
            whether or not to plot in background 
            (meaning that Python processes in terminal can run while plot is open) 
            if True, pyvistaqt is used to generate plot (different interface)
     
        Returns
        ----------
        pl : `pyvista.Plotter()` object
        
        Example
        ---------
        See separate notebook on GitHub for example.

        '''

        if node_label_fun is None:
            node_label_fun = lambda n: str(int(n.label))
        if element_label_fun is None:
            element_label_fun = lambda e: str(int(e.label))

        canvas_settings = dict(background_color='white') | canvas_settings

        if 'deformed' in plot_states:
            add_transparency = dict(opacity=0.4)
        else:
            add_transparency = {}

        # Element label settings
        elementlabel_settings = dict(always_visible=True, 
                                      text_color='blue', 
                                      shape_color='white', 
                                      shape_opacity=0.2) | elementlabel_settings

        # Node and sensor label settings
        nodelabel_settings = dict(always_visible=True, shape_opacity=0.2)
        sensorlabel_settings = nodelabel_settings | sensorlabel_settings

        # Node plotting settings
        node_settings = dict(
                            render_points_as_spheres=True,
                            color='black',
                            lighting=True,
                            point_size=6) | add_transparency | node_settings  
        deformed_node_settings = node_settings | deformed_node_settings
        
        # Line plotting settings
        line_settings = dict(render_lines_as_tubes=True,
                                style='wireframe',
                                lighting=True,
                                line_width=4) | add_transparency | line_settings
        
        deformed_line_settings = line_settings | dict(color='#ee8899') | deformed_line_settings

        # Face plotting settings
        face_settings = dict() | add_transparency | face_settings
        deformed_face_settings = face_settings | dict(color='#ee8899') | deformed_face_settings


        if plot_nodes is not False:
            if plot_nodes is True:
                nodes_to_plot = self.nodes*1
            else:
                nodes_to_plot = [node for node in self.nodes if node in plot_nodes]

        if pl is None:
            if background_plotter:
                pl = pvqt.BackgroundPlotter()
            else:
                pl = pv.Plotter()

        if node_labels is not False:
            if node_labels is not True:
                labeled_nodes = [node for node in self.nodes if node in node_labels]
            else:
                labeled_nodes = self.nodes

            lbl = [node_label_fun(node) for node in labeled_nodes]
            
            if 'deformed' in plot_states and (label_priority!='undeformed' or 'undeformed' not in plot_states):
                pl.add_point_labels(self.get_points(nodes=labeled_nodes, deformed=False, flattened=False), lbl, **nodelabel_settings)
            if 'undeformed' in plot_states and (label_priority!='deformed' or 'deformed' not in plot_states):
                pl.add_point_labels(self.get_points(nodes=labeled_nodes, deformed=True, flattened=False), lbl, **nodelabel_settings)            

        if element_labels is not False:
            if element_labels is not True:
                labeled_els = [el for el in self.elements if el in element_labels]
            else:
                labeled_els = self.elements

            lbl = [node_label_fun(el) for el in labeled_els]
            
            if 'deformed' in plot_states and (label_priority!='undeformed' or 'undeformed' not in plot_states):
                pl.add_point_labels(self.get_element_cogs(elements=labeled_els, deformed=False, flattened=False), lbl, **elementlabel_settings)
            if 'undeformed' in plot_states and (label_priority!='deformed' or 'deformed' not in plot_states):
                pl.add_point_labels(self.get_element_cogs(elements=labeled_els, deformed=True, flattened=False), lbl, **elementlabel_settings)            
                    
        if sensor_labels:
            for sensor in self.sensors:
                node = self.get_node(self.sensors[sensor])
                if 'undeformed' in plot_states and (label_priority!='deformed' or 'deformed' not in plot_states):
                    pl.add_point_labels(np.array([node.xyz0]), [sensor], **sensorlabel_settings)     
                if 'deformed' in plot_states and (label_priority!='undeformed' or 'undeformed' not in plot_states):
                    pl.add_point_labels(np.array([node.xyz]), [sensor], **sensorlabel_settings)     
                    

        if plot_lines:
            if 'undeformed' in plot_states:
                lines = self.get_lines()
                mesh = pv.PolyData(self.get_points(), 
                            lines=lines)        
                pl.add_mesh(mesh, **line_settings)
            if 'deformed' in plot_states:
                lines = self.get_lines()
                mesh = pv.PolyData(self.get_points(deformed=True), 
                            lines=lines)        
                pl.add_mesh(mesh, **deformed_line_settings)

        if plot_faces:
            if 'undeformed' in plot_states:
                faces = self.get_faces()
                mesh = pv.PolyData(self.get_points(),
                            faces=faces)
            
                pl.add_mesh(mesh, **face_settings)
            if 'deformed' in plot_states: 
                faces = self.get_faces()
                mesh = pv.PolyData(self.get_points(deformed=True),
                            faces=faces)
            
                pl.add_mesh(mesh, **deformed_face_settings)

        if view is not None:
            if view in ['xy', 'top']:
                pl.view_xy()
            if view in ['yz', 'front']:
                pl.view_yz()
            if view in ['xz', 'side']:
                pl.view_xz()

        if plot_nodes:
            if 'undeformed' in plot_states:
                pl.add_points(self.get_points(nodes=nodes_to_plot, deformed=False, flattened=False), **node_settings)
            if 'deformed' in plot_states:
                pl.add_points(self.get_points(nodes=nodes_to_plot, deformed=True, flattened=False), **deformed_node_settings)

        for key in canvas_settings:
            setattr(pl, key, canvas_settings[key])

        pl.camera.SetParallelProjection(not perspective_cam)
   
        if show:
            pl.show()

        return pl

    
    def animate_mode(self, phi, filename=None, pl=None, add_undeformed=False, undeformed_settings={'opacity':0.2},
                     cycles=1, f=1.0, fps=60, **kwargs):
        
        '''
        Plot model (either undeformed or deformed).

        Arguments
        ----------
        phi : float, complex
            1d array with considered mode shape
            complex values supported
        filename : None, str, optional
            path to filename, use either `mp4` or `gif` file types (additional packages may be required)
            if `None`, interactive mode is assumed
        pl : `pyvista.Plotter()` object, optional
            plotter object to base animation on - one will be created if not input
        add_undeformed : False, optional
            whether or not to add undeformed (reference) structure in animation
            only applicable if `pl=None`
        undeformed_settings : dict, optional
            settings dict used for both lines, edges and nodes of undeformed structure
            keywords must be valid for all these object types (therefore limited flexibility)
            only applicable if `pl=None` and `add_undeformed=True`
        cycles : 1, optional
            number of cycles for stored file (not used for interactive mode)
        f : 1.0, optional
            frequency of rotation for stored file (not used for interactive mode)
        fps : 60, optional
            frames per second used for stored file (not used for interactive mode)
        **kwargs
            arguments supported by `spatial.plot` are passed to that method
        
        Example
        ---------
        See separate notebook on GitHub for example.

        '''
        
        frames_per_cycle = int(np.round(fps/f))

        # Save video
        if filename is not None:

            if pl is None:
                if add_undeformed:
                    pl = self.plot(deformed=False, background_plotter=False, show=False, 
                                        line_settings=undeformed_settings, node_settings=undeformed_settings,
                                        face_settings=undeformed_settings, **kwargs)
                else:
                    pl = pv.Plotter()
                    pl.background_color='white'


            if filename.split('.')[1].lower()=='gif':
                pl.open_gif(filename, fps=fps)
            else:
                pl.open_movie(filename)

            pl = self.plot(pl=pl, deformed=True, **kwargs)
            
            frames = cycles*frames_per_cycle
            dtheta = 2*np.pi * f/fps
            theta = 0.0

            for frame in range(frames):
                theta = theta+dtheta
                self.u = np.real(phi * np.exp(theta*1j))
                pts = pv.pyvista_ndarray(self.get_points(deformed=True, flattened=False))
                self.face_mesh.points = pts
                self.line_mesh.points = pts
                self.point_mesh.points = pts
                pl.update()
                pl.write_frame()

            pl.close()

        else:   # interactive animation           
            # Initial deformed plot

            def update_shape():
                self.step += 1
                self.u = np.real(phi * np.exp(2*np.pi*1j*(self.step/frames_per_cycle)))
                pts = pv.pyvista_ndarray(self.get_points(deformed=True, flattened=False))
                self.face_mesh.points = pts
                self.line_mesh.points = pts
                self.point_mesh.points = pts
                pl.update()

            if pl is None:
                if add_undeformed:
                    pl = self.plot(deformed=False, background_plotter=True, show=False, 
                                        line_settings=undeformed_settings, node_settings=undeformed_settings,
                                        face_settings=undeformed_settings, **kwargs)
                else:
                    pl = pvqt.BackgroundPlotter()
                    pl.background_color='white'

            pl = self.plot(pl=pl, deformed=True, **kwargs)

            self.step = 0
            pl.add_callback(update_shape, interval=0)  
            pl.show()
            pl.app.exec_()
        
            sys.exit()


    