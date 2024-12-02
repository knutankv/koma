'''
MODULE FOR SPATIAL MAPPING AND VISUALIZATION
'''
import pyvista as pv
import pyvistaqt as pvqt
import numpy as np
from koma.modal import maxreal_vector
import sys
from scipy.interpolate import griddata

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

        # Initialize deformed
        self.x = self.x0*1
        self.y = self.y0*1
        self.z = self.z0*1
        self._label = int(label)
        
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
    
    @property
    def u(self):
        return self.xyz-self.xyz0

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
    x_ixs : None, optional
        indices of input model DOFs that correspond to global x
        used to construct interpolation field for automatic extension of mode shapes
        when None (standard value), 0::3 is used
    y_ixs : None, optional
        indices of input model DOFs that correspond to global x
        used to construct interpolation field for automatic extension of mode shapes
        when None (standard value), 1::3 is used
    z_ixs : None, optional
        indices of input model DOFs that correspond to global x
        used to construct interpolation field for automatic extension of mode shapes
        when None (standard value), 2::3 is used
    undefined_dofs : {'zero', 'linear', 'quadratic'}, optional
        how to treat undefined nodes (either set to zero or interpolate between DOFs present)  
        interpolation (values other than 'zero') is experimental
    interpolation_axes : [0,1,2], optional
        axes used for interpolation

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
    def __init__(self, nodes, elements, dofmap=None, sensors={}, 
                 u=None, x_ixs=None, y_ixs=None, z_ixs=None, 
                 undefined_dofs='zero', interpolation_axes=[0,1,2]):
        
        self.nodes = nodes
        self.elements = elements
        self.dofmap = dofmap
        self.sensors = sensors
        self.u = u

        if x_ixs is None:
            self.x_ixs = np.arange(0, len(self.nodes)*3, 3)
        else:
            self.x_ixs = x_ixs

        if y_ixs is None:
            self.y_ixs = np.arange(1, len(self.nodes)*3, 3)
        else:
            self.y_ixs = y_ixs

        if z_ixs is None:
            self.z_ixs = np.arange(2, len(self.nodes)*3, 3)
        else:
            self.z_ixs = z_ixs

        self.interpolation_axes = np.array(interpolation_axes)

        # Ensure node objects and not labels
        for el in self.elements:
            el.nodes = [self.get_node(n) for n in el.nodes] 
        
        self.undefined_dofs = undefined_dofs

    @property 
    def u(self):
        return self._u
    
    @u.setter
    def u(self, val):
        self.ufull = self.expand_field(val)
        self.deform(np.real(self.ufull))
        self._u = val

    def rotate_phase(self, angle):
        self.ufull = self.ufull*np.exp(1j*angle)
        self.deform(np.real(self.ufull))
    
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
        if u is None:
            return np.zeros(len(self.nodes)*3)

        ufull = np.zeros(len(self.nodes)*3).astype(complex)*np.nan

        if u is not None:
            # First, assign all direct (if not given - stays zero)
            for node in self.dofmap:
                rels = self.dofmap[node]
                for dof_ix, rel in enumerate(rels):
                    if type(rel) in [int, float, np.int32]:
                        global_ix = self.get_node_ixs(node)[dof_ix]
                        if rel is None:
                            ufull[global_ix] = np.nan
                        else:
                            ufull[global_ix] = u[rel]
                        

            # Second, assign all relative
            n = self.grab_fun(ufull)
            for node in self.dofmap:
                rels = self.dofmap[node]
                for dof_ix, rel in enumerate(rels):
                    if type(rel) == Rel:
                        global_ix = self.get_node_ixs(node)[dof_ix]
                        ufull[global_ix] = rel.fun(n)

        
        if self.undefined_dofs != 'zero':
            ufull[self.x_ixs] = self.fill_from_interpolation(self.x_ixs, 
                                                                ufull[self.x_ixs],
                                                                method=self.undefined_dofs)
            
            ufull[self.y_ixs]  = self.fill_from_interpolation(self.y_ixs, 
                                                                ufull[self.y_ixs], 
                                                                method=self.undefined_dofs)
            
            ufull[self.z_ixs]  = self.fill_from_interpolation(self.z_ixs, 
                                                                ufull[self.z_ixs],
                                                                method=self.undefined_dofs)

        ufull[np.isnan(ufull)] = 0.0

        return ufull

    def fill_from_interpolation(self, ixs, u, method='linear'):
        source_nodes = self.get_nodes_with_dofs(ixs[~np.isnan(u)])
        target_nodes = self.get_nodes_with_dofs(ixs[np.isnan(u)])

        if len(source_nodes)>1:
            source_xyz = np.vstack([node.xyz0 for node in source_nodes])[:, self.interpolation_axes]
            source_u = np.array(u[~np.isnan(u)])
            target_xyz = np.vstack([node.xyz0 for node in target_nodes])[:, self.interpolation_axes]
            target_u = griddata(source_xyz, source_u, target_xyz, method=method, fill_value=0).flatten()
    
            u[np.isnan(u)] = target_u
            
        return u


    def get_nodes_with_dofs(self, ixs):
        '''
        Establish list of nodes conatining given DOF indices (all nodes have 3 DOFs).
        '''

        return [self.nodes[int(ix)] for ix in np.floor(ixs/3)]

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
            out = [n.xyz for n in nodes]
        else:
            out = [n.xyz0 for n in nodes]
        
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
    
    @property
    def el_types(self):
        return set([el.el_type for el in self.elements])

    def plot(self, pl=None, show=True, plot_lines=True, plot_nodes=True, 
                plot_sensor_nodes=True, 
                node_labels=False, element_labels=False,
                plot_faces=True, canvas_settings={}, 
                node_settings={}, line_settings={}, 
                face_settings={}, nodelabel_settings={}, sensorlabel_settings={},
                sensor_node_settings={'color':'red', 'point_size': 8}, elementlabel_settings={}, view=None, sensor_labels=False,
                deformed=True, node_label_fun=None, element_label_fun=None, perspective_cam=False, background_plotter=True):


        '''
        Plot model (either undeformed or deformed).

        Arguments
        ----------
        pl : `pyvista.Plotter()` object, optional
            plotter object - one will be created if not input
        show : True, optional
            whether or not to show plot at end - output is pl-object, which can be shown later
        plot_lines : True, optional
            whether or not to plot the line elements;
            only applicable if line elements are present
        plot_nodes : True, optional
            whether or not to plot the nodes
        node_labels : True, optional
            whether or not to show node labels
            if input is a list of node labels, only these are shown
        element_labels :True, optional
            whether or not to show element labels
            if input is a list of element labels, only these are shown
        plot_faces : True, optional
            whether or not to plot face elements (triangles and rectangles);
            only applicable if face elements are present
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
        plot_sensor_nodes : True, optional
            whether or not to plot sensors
        sensor_node_settings : dict, optional
            overriding settings dict on sensor nodes
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

        # Check what type of elements are present
        el_types = self.el_types
        plot_faces = plot_faces & (('triangle' in el_types) or ('rectangle' in el_types))
        plot_lines = plot_lines & ('line' in el_types)

        # Element label settings
        elementlabel_settings = dict(always_visible=True, 
                                      text_color='blue', 
                                      shape_color='white', 
                                      shape_opacity=0.2) | elementlabel_settings

        # Node and sensor label settings
        nodelabel_settings = dict(always_visible=True, shape_opacity=0.2, text_color='black') | nodelabel_settings
        sensorlabel_settings = nodelabel_settings | sensorlabel_settings

        # Node plotting settings
        node_settings = dict(
                            render_points_as_spheres=True,
                            color='black',
                            lighting=True,
                            point_size=6) | node_settings 
         
        sensor_node_settings = dict(node_settings) | sensor_node_settings

        # Line plotting settings
        line_settings = dict(render_lines_as_tubes=True,
                                style='wireframe',
                                lighting=True,
                                line_width=4) | line_settings
        
        canvas_settings = dict(background_color='white') | canvas_settings

        # Face plotting settings
        face_settings = dict(show_edges=True, color='lightgray') | face_settings

        if plot_nodes is True:
            nodes_to_plot = self.nodes*1
        elif plot_nodes is False and sensor_node_settings != node_settings and hasattr(self, 'sensors') and len(self.sensors)>0:
            nodes_to_plot = self.get_nodes(list(self.sensors.values()))
        elif plot_nodes is False:
            nodes_to_plot = []
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
            pl.add_mesh(self.line_mesh, **line_settings)

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

        # Sensor nodes
        if plot_sensor_nodes and len(self.sensors)>0:
            sensor_node_points = pv.pyvista_ndarray(self.get_points(nodes=self.sensors.values(), 
                                                        deformed=deformed, 
                                                        flattened=False))
            sensor_pts = pl.add_points(sensor_node_points, **sensor_node_settings)
            self.sensor_point_mesh = sensor_pts.mapper.dataset
        
        pl.camera.SetParallelProjection(not perspective_cam)

        if show:
            pl.show()

        return pl


    
    def animate_mode(self, phi, filename=None, pl=None, add_undeformed=False, 
                     undeformed_settings={'opacity':0.2},
                     node_settings={}, face_settings={}, line_settings={},
                     cycles=1, f=1.0, fps=60, ensure_exact_cycle=True, **kwargs):
        
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
        node_settings : dict, optional
            dictionary with additional settings for undeformed nodes / points
            given to `pyvista.Plotter().add_points`
        line_settings : dict, optional
            dictionary with additional settings for undeformed line elements
            given to `pyvista.Plotter().add_mesh`
        face_settings : dict, optional
            dictionary with additional settings for undeformed patch elements
            given to `pyvista.Plotter().add_mesh`
        cycles : 1, optional
            number of cycles for stored file (not used for interactive mode)
        f : 1.0, optional
            frequency of rotation for stored file (not used for interactive mode)
        fps : 60, optional
            frames per second used for stored file (not used for interactive mode)
        ensure_exact_cycle : True, optional
            ensure that fps is such that 0 and 0+T both are exactly sampled/shown
        **kwargs
            arguments supported by `spatial.plot` are passed to that method
        
        Example
        ---------
        See separate notebook on GitHub for example.

        '''
        
        # Ensure repeating pattern
        if ensure_exact_cycle:
            fps = np.ceil(fps/f) * f

        frames_per_cycle = int(np.round(fps/f))
        self.u = maxreal_vector(phi)
        self.dt = 1/fps
        self.dangle = self.dt * f * np.pi * 2.0

        # Update function
        def update_shape():
            self.rotate_phase(self.dangle)
            pts = pv.pyvista_ndarray(self.get_points(deformed=True, flattened=False))
            
            if hasattr(self, 'face_mesh'):
                 self.face_mesh.points = pts
           
            if hasattr(self, 'line_mesh'):
                self.line_mesh.points = pts

            if hasattr(self, 'point_mesh'):
                self.point_mesh.points = pts

            if hasattr(self, 'sensor_point_mesh'):
                pts_sensors = pv.pyvista_ndarray(self.get_points(nodes=self.sensors.values(), 
                                                                 deformed=True, flattened=False))
                self.sensor_point_mesh.points = pts_sensors

            pl.update()

        # Save video
        if filename is not None:
            if pl is None:
                if add_undeformed:
                    pl = self.plot(deformed=False, background_plotter=False, show=False, 
                                        line_settings=line_settings | undeformed_settings, node_settings=node_settings | undeformed_settings,
                                        face_settings=face_settings | undeformed_settings)
                else:
                    pl = pv.Plotter()
                    pl.background_color='white'

            if filename.split('.')[-1].lower()=='gif':
                pl.open_gif(filename, fps=fps)
            else:
                pl.open_movie(filename)

            pl = self.plot(pl=pl, show=False, deformed=True, line_settings=line_settings, 
                                node_settings=node_settings,
                                face_settings=face_settings, **kwargs)
            
            pl.show(interactive_update=True)
            frames = cycles*frames_per_cycle
            
            for frame in range(frames):
                update_shape()
                pl.write_frame()
            
            pl.close()

        else:   # Interactive animation           

            if pl is None:
                if add_undeformed:
                    pl = self.plot(deformed=False, background_plotter=True, show=False, 
                                        line_settings=undeformed_settings, node_settings=undeformed_settings,
                                        face_settings=undeformed_settings, **kwargs)
                else:
                    pl = pvqt.BackgroundPlotter()
                    pl.background_color='white'

            pl = self.plot(pl=pl, deformed=True, line_settings=line_settings, 
                                node_settings=node_settings,
                                face_settings=face_settings, **kwargs)
            self.dt = 1/fps
            self.dangle = self.dt * f * np.pi * 2.0
            pl.add_callback(update_shape, interval=int(np.ceil(self.dt*1000)))  
            pl.show()
            pl.app.exec_()
        
            sys.exit()

    

def save_model(model, path):
    import dill
    with open(path, 'wb') as f:
        dill.dump(model, f, -1)

def load_model(path):
    import dill
    with open(path, 'rb') as f:
        model = dill.load(f)
    return model