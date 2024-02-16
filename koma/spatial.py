'''
(BASIC) MODULE FOR SPATIAL MAPPING AND VISUALIZATION
'''
import pyvista as pv
import pyvistaqt as pvqt
import numpy as np

def nodes_from_matrix(node_matrix):
    return [Node(*row) for row in node_matrix]

def elements_from_matrix(element_matrix):
    return [Element(el[1:], label=el[0]) for el in element_matrix]

class Element:
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
    def __init__(self, label, x, y, z=0):
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

    Arguments
    ------------
    nodes : `Node` 
        list of nodes
    elements : `Line` or `Triangle`
        list of elements (either lines or triangle patches)
    dofmap : dict optional
        doc

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
        for node in self.nodes:
            dxyz = ufull[self.get_node_ixs(node)]
            node.x, node.y, node.z = dxyz[0]+node.x0, dxyz[1]+node.y0, dxyz[2]+node.z0

    def expand_field(self, u):
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
             node_settings={}, deformed_node_settings={},
             line_settings={}, deformed_line_settings={}, 
             face_settings={}, deformed_face_settings={}, 
             nodelabel_settings={}, sensorlabel_settings={},
             elementlabel_settings={}, view=None, sensor_labels=False,
             plot_states=['deformed','undeformed'], label_priority='deformed',
             node_label_fun=None, element_label_fun=None, perspective_cam=True, background_plotter=True):

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

    
    def animate_mode(self, phi, frequency=1.0, n_cycles=1, 
                    frames_per_cycle=30, show=True,  **kwargs):
        omega = frequency*2*np.pi
        pl = self.plot(show=False, use_pvqt=True, **kwargs)

        duration = n_cycles*frames_per_cycle
        def update_shape(step):
            self.u = np.real(phi * np.exp(1j*omega*(step/duration)))
            print(self.u)
            pl.render()

        timer = pvqt.plotting.QTimer()
        timer.timeout.connect(update_shape)
        timer.setInterval(300)
        timer.setDuration(1000)
        timer.start()